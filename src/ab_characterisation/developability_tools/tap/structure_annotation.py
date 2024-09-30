import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio import PDB
from Bio.PDB.NeighborSearch import NeighborSearch

from ab_characterisation.developability_tools.tap.definitions import (
    acceptors, anchor_residues, donors, imgt_cdr_definitions,
    normalised_hydrophobicities, residue_charges)


class PSAError(Exception):
    """Raised when something has gone awry when running psa to calculate surface areas."""


@dataclass
class AnnotatedResidue(PDB.Residue.Residue):  # type: ignore
    salt_bridge_partner: Optional[tuple] = None
    neighbours: dict = field(default_factory=dict)
    cdr_number: int = field(init=False)
    relative_surface_area: float = field(init=False)
    hydrophobicity: float = field(init=False)
    charge: float = field(init=False)
    in_cdr_vicinity: bool = field(init=False)

    def __eq__(self, other):  # type: ignore
        """Direct copy from Biopython Entity"""
        if isinstance(other, type(self)):
            if self.parent is None:
                return self.id == other.id
            return self.full_id[1:] == other.full_id[1:]
        return NotImplemented

    def __hash__(self) -> int:
        """Direct copy from Biopython Entity"""
        return hash(self.full_id)

    @property
    def is_cdr(self) -> bool:
        """Whether the residue is part of a CDR (IMGT definition) or not."""
        return self.cdr_number > 0

    @property
    def is_anchor(self) -> bool:
        """
        Whether the residue an anchor residue to a CDR (IMGT definition) or not.
        Anchor residues are defined as the two residues on each side of the CDR.
        """
        return self.id[1] in anchor_residues

    @property
    def is_surface(self) -> bool:
        """
        Uses the relative surface area as calculated by psa to determine whether the residue is on the surface or not.
        Surface residues have a relative sidechain surface area of 7.5 or above.
        """
        return self.relative_surface_area >= 7.5

    @property
    def res_number(self) -> str:
        """Returns a formatted string containing the residue number and insertion code, if present."""
        return f"{self.id[1]}{self.id[2]}".strip()

    @property
    def is_donor(self) -> bool:
        """Whether the residue is a salt bridge donor residue type."""
        return self.resname in donors

    @property
    def is_acceptor(self) -> bool:
        """Whether the residue is a salt bridge acceptor residue type."""
        return self.resname in acceptors


@dataclass
class StructureAnnotator:
    """
    Class containing methods for structural annotation of properties required by TAP":
        - relative surface area
        - minimum distances between neighbouring residues
        - CDR vicinity (surface residues within 4 A of CDRs/anchors)
        - salt bridges (donor/acceptor atoms within 3.2 A)
        - hydrophobicity
        - charge
    """

    neighbour_cutoff: float = 7.5
    salt_bridge_cutoff: float = 3.2
    vicinity_cutoff: float = 4.0
    psa_path: Path = field(init=False)
    cdr_lookup_dict: dict[tuple[str, int], int] = field(init=False)

    def __post_init__(self) -> None:
        lookup_dict = {}
        for chain in "HL":
            for cdr, residue_range in imgt_cdr_definitions.items():
                for res in residue_range:
                    lookup_dict[(chain, res)] = cdr
        self.cdr_lookup_dict = lookup_dict

        psa_version = "psa_mac" if sys.platform == "darwin" else "psa"
        self.psa_path = (
            Path(__file__).resolve().parent / f"psa_executables/{psa_version}"
        )

    @staticmethod
    def _convert_residues(structure: PDB.Structure.Structure) -> None:
        """Converts normal Biopython residues to our annotated version with extra properties/methods."""
        for res in structure[0].get_residues():
            res.__class__ = AnnotatedResidue
            res.neighbours = {}
        return

    def _run_psa(self, structure_path: str) -> list[str]:
        """Runs the psa executable on the .pdb file to get surface accessibility information."""
        if self.psa_path.exists() is False:
            raise PSAError("psa executable was not found.")

        result, error = subprocess.Popen(
            [str(self.psa_path), "-t", structure_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).communicate()
        if not result:
            raise PSAError(error.decode())

        psa_output = result.decode().split("\n")
        return psa_output

    def _annotate_sasa(
        self, structure: PDB.Structure.Structure, structure_path: str
    ) -> None:
        """Runs psa and extracts surface accessibility information from its output."""
        # Run psa and get the relevant lines from the output
        psa_output = self._run_psa(structure_path)
        residue_lines = [line for line in psa_output if line.startswith("ACCESS")]

        # Check that the number of residues in the psa output is the same as the number of residues in our structure
        all_residues = list(structure[0].get_residues())
        if len(all_residues) != len(residue_lines):
            raise PSAError("PSA output contained the wrong number of residues.")

        # Iterate through residues and annotate the structure
        for res, psa_line in zip(all_residues, residue_lines):
            # Check we are on the correct residue with the correct type
            psa_residue_number = psa_line[6:12].strip()
            if res.res_number != psa_residue_number:
                raise PSAError(
                    f"Residue number mismatch: {res.res_number} != {psa_residue_number}"
                )
            psa_residue_type = psa_line[14:17]
            if res.resname != psa_residue_type:
                raise PSAError(
                    f"Expected type {res.resname} for residue {res.parent}{res.res_number}; got {psa_residue_type}"
                )

            res.relative_surface_area = float(psa_line[61:67])
        return

    @staticmethod
    def _get_minimum_distance(res1: AnnotatedResidue, res2: AnnotatedResidue) -> float:
        """Calculates the minimum distance between the heavy atoms of two residues."""
        min_dist = 100.0
        for atom1 in res1.get_atoms():
            for atom2 in res2.get_atoms():
                dist = atom1 - atom2
                if dist < min_dist:
                    min_dist = dist
        return min_dist

    def _get_neighbours(self, structure: PDB.Structure.Structure) -> None:
        """For each residue in the structure, gets a list of neighbouring residues and their minimum distance"""
        # Quickly get list of residue pairs that are less than 7.5A apart
        all_heavy_atoms = [
            atom for atom in structure[0].get_atoms() if atom.element != "H"
        ]
        residue_pairs = NeighborSearch(atom_list=all_heavy_atoms).search_all(
            self.neighbour_cutoff, level="R"
        )

        for res1, res2 in residue_pairs:
            if res1 == res2:
                continue
            min_dist = self._get_minimum_distance(res1, res2)
            res1.neighbours[res2.get_full_id()] = min_dist
            res2.neighbours[res1.get_full_id()] = min_dist
        return

    def _cdr_lookup(self, chain_id: str, residue_number: int) -> int:
        """
        Returns the number of the CDR a residue is part of from its residue number.
        Returns zero if the residue is not part of a CDR.
        """
        return self.cdr_lookup_dict.get((chain_id, residue_number), 0)

    def _annotate_cdrs(self, structure: PDB.Structure.Structure) -> None:
        """Annotates residues with their CDR number (0 if not in a CDR)"""
        for res in structure[0].get_residues():
            chain_id = res.parent.id
            res_number = res.id[1]
            res.cdr_number = self._cdr_lookup(chain_id, res_number)

    def _annotate_cdr_vicinity(self, structure: PDB.Structure.Structure) -> None:
        """Finds and annotates which residues are on the surface and less than 4A away from the CDRs/anchors."""
        surface_cdrs_and_anchors = []
        for res in structure[0].get_residues():
            if res.is_surface and (res.is_cdr or res.is_anchor):
                res.in_cdr_vicinity = True
                surface_cdrs_and_anchors.append(res)
            else:
                res.in_cdr_vicinity = False

        for res in surface_cdrs_and_anchors:
            res.in_cdr_vicinity = True
            for neighbour_id, distance in res.neighbours.items():
                res2 = structure[0][neighbour_id[2]][neighbour_id[3]]
                if distance < self.vicinity_cutoff and res2.is_surface:
                    res2.in_cdr_vicinity = True
        return

    def _annotate_salt_bridges(self, structure: PDB.Structure.Structure) -> None:
        """Identifies which residues form salt bridges based on distance."""
        donor_atoms = [
            atom
            for atom in structure[0].get_atoms()
            if atom.id in donors.get(atom.parent.resname, [])
        ]
        acceptor_atoms = [
            atom
            for atom in structure[0].get_atoms()
            if atom.id in acceptors.get(atom.parent.resname, [])
        ]
        residue_pairs = NeighborSearch(
            atom_list=donor_atoms + acceptor_atoms
        ).search_all(self.salt_bridge_cutoff, level="R")

        for res1, res2 in residue_pairs:
            # Ignore if already part of a salt bridge
            if res1.salt_bridge_partner or res2.salt_bridge_partner:
                continue
            if res1.is_surface and res2.is_surface:
                if (res1.is_donor and res2.is_acceptor) or (
                    res2.is_donor and res1.is_acceptor
                ):
                    res1.salt_bridge_partner = res2.get_full_id()
                    res2.salt_bridge_partner = res1.get_full_id()
        return

    def _annotate_hydrophobicity(self, structure: PDB.Structure.Structure) -> None:
        """
        Annotates residues with their normalised (between 1 and 2) hydrophobicity values.
        If the residue forms part of a salt bridge, it is assigned the hydrophobicity of glycine.
        """
        for res in structure[0].get_residues():
            if res.salt_bridge_partner:
                res.hydrophobicity = normalised_hydrophobicities["GLY"]
            else:
                res.hydrophobicity = normalised_hydrophobicities[res.resname]
        return

    def _annotate_charge(self, structure: PDB.Structure.Structure) -> None:
        """
        Annotates residues with their charges.
        If the residue forms part of a salt bridge, it is assigned a charge of zero.
        """
        for res in structure[0].get_residues():
            if res.salt_bridge_partner:
                res.charge = 0
            else:
                res.charge = residue_charges[res.resname]
        return

    def load_and_annotate_structure(
        self, structure_path: str
    ) -> PDB.Structure.Structure:
        """
        Loads an antibody structure from the provided file, and annotates the residues for later use in TAP metric
        calculations.
        Assumes the structure is already IMGT-numbered!!

        Args:
            structure_path: the path to the structure that is to be annotated.

        Returns:
            the annotated structure (Biopython Structure entity, with residues converted to an AnnotatedResidue type).
        """
        structure = PDB.PDBParser(QUIET=True).get_structure(
            "input_structure", structure_path
        )
        self._convert_residues(structure)
        self._annotate_sasa(structure, structure_path)
        self._get_neighbours(structure)
        self._annotate_cdrs(structure)
        self._annotate_cdr_vicinity(structure)
        self._annotate_salt_bridges(structure)
        self._annotate_hydrophobicity(structure)
        self._annotate_charge(structure)
        return structure
