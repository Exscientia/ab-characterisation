import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional

from loguru import logger

from ab_characterisation.developability_tools.sequence_liabilities.definitions import \
    custom_regions
from ab_characterisation.utils.anarci_utils import Accept


@dataclass
class Position:
    chain: str
    number: int
    ins_code: str

    def to_string(self) -> str:
        if self.ins_code != " ":
            return f"{self.chain}{self.number}{self.ins_code}"
        return f"{self.chain}{self.number}"


@dataclass
class SequenceLiability:
    liability_type: str
    motif: str
    positions: list[Position]

    @property
    def positions_string(self) -> str:
        return "-".join([pos.to_string() for pos in self.positions])


@dataclass
class BaseScannerDataclassMixin:
    name: str
    description: str


class BaseScanner(ABC, BaseScannerDataclassMixin):
    name: str
    description: str

    @abstractmethod
    def scan(
        self,
        numbering_dict: dict[str, list[tuple[tuple[int, str], str]]],
        quiet: bool = False,
    ) -> List[SequenceLiability]:
        """
        Scans an input sequence for liabilities.

        Args:
            numbering_dict: a dictionary of ANARCI numberings
                       e.g. {"H": [((1, ' '), 'E'), ((2, ' '), 'L'), ((3, ' '), 'K'), ...],
                             "L": [((1, ' '), 'D'), ((2, ' '), 'V'), ((3, ' '), 'L'), ...]}

        Returns:
            a list of identified liabilities.
        """


@dataclass
class RegexScanner(BaseScanner):
    regions: list[str]
    regex_search_string: str
    ignored_positions: Optional[list[tuple[int, str]]] = None

    def __post_init__(self) -> None:
        self.regex_pattern = re.compile(self.regex_search_string)

    def _get_acceptor(self, chain: str) -> Accept:
        acceptor = Accept(numbering_scheme="imgt", definition="imgt")
        for region in self.regions:
            if region in custom_regions:
                acceptor.add_positions(custom_regions[region][chain], chain)
            else:
                acceptor.add_regions([region])
        return acceptor

    def scan(
        self,
        numbering_dict: dict[str, list[tuple[tuple[int, str], str]]],
        quiet: bool = False,
    ) -> List[SequenceLiability]:
        identified = []
        for chain, numbering in numbering_dict.items():
            acceptor = self._get_acceptor(chain)

            sequence = "".join([res[1] for res in numbering if res[1] != "-"])
            numbers = [res[0] for res in numbering if res[1] != "-"]

            for match in self.regex_pattern.finditer(sequence):
                start, end = match.start(), match.end()
                identified_positions = numbers[start:end]

                # Check if any of the residues identified should be ignored; skip if so
                if self.ignored_positions:
                    if set(identified_positions).intersection(self.ignored_positions):
                        continue

                # Check if the first residue of the set identified belongs to a region of interest
                if acceptor.accept(identified_positions[0], chain):
                    identified.append(
                        SequenceLiability(
                            liability_type=self.name,
                            motif=match.group(),
                            positions=[
                                Position(chain=chain, number=pos[0], ins_code=pos[1])
                                for pos in identified_positions
                            ],
                        )
                    )

        if not quiet:
            color = "red" if identified else "green"
            logger.opt(colors=True).info(
                f"<b><{color}>{self.name}:</{color}></b> identified <b><{color}>{len(identified)}</{color}></b> liabilities"
            )

        return identified


class NTerminalGlutamateScanner(BaseScanner):
    # This does not look for a consecutive pattern like the other liabilities
    # Checks for E residues at the start of each chain instead
    def scan(
        self,
        numbering_dict: dict[str, list[tuple[tuple[int, str], str]]],
        quiet: bool = False,
    ) -> List[SequenceLiability]:
        if "H" not in numbering_dict or "L" not in numbering_dict:
            if not quiet:
                logger.opt(colors=True).warning(
                    f"<y><b>{self.name}</b></y>: both H and L chain sequences are required for this check; skipping"
                )
            return []

        heavy_dict: dict[tuple[int, str], str] = dict(numbering_dict["H"])
        light_dict: dict[tuple[int, str], str] = dict(numbering_dict["L"])

        identified = []
        if (
            heavy_dict.get((1, " "), None) == "E"
            and light_dict.get((1, " "), None) == "E"
        ):
            identified = [
                SequenceLiability(
                    liability_type=self.name,
                    motif="EE",
                    positions=[
                        Position(chain="H", number=1, ins_code=" "),
                        Position(chain="L", number=1, ins_code=" "),
                    ],
                )
            ]

        if not quiet:
            color = "red" if identified else "green"
            logger.opt(colors=True).info(
                f"<b><{color}>{self.name}:</{color}></b> identified <b><{color}>{len(identified)}</{color}></b> liabilities"
            )

        return identified
