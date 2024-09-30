from dataclasses import dataclass, field
from typing import Optional, Union

from ab_characterisation.utils.anarci_region_definition_utils import (_index_to_imgt_state,
                                                      _regions)

_reg_one2three = {
    "1": "fw%s1",
    "2": "cdr%s1",
    "3": "fw%s2",
    "4": "cdr%s2",
    "5": "fw%s3",
    "6": "cdr%s3",
    "7": "fw%s4",
}


@dataclass
class Accept:  # pylint: disable=R0902
    """
    Class that taken ANARCI numbering and classifies each position according to antibody region
    """

    _defined_regions: list[str] = field(
        init=False,
        repr=False,
        default_factory=lambda: [
            "fwh1",
            "fwh2",
            "fwh3",
            "fwh4",
            "fwl1",
            "fwl2",
            "fwl3",
            "fwl4",
            "cdrh1",
            "cdrh2",
            "cdrh3",
            "cdrl1",
            "cdrl2",
            "cdrl3",
        ],
    )
    numbering_scheme: str = field(default="imgt")
    definition: str = field(default="imgt")
    not_defined: bool = field(default=False)
    positions: dict[str, set[tuple[int, str]]] = field(init=False)
    exclude: dict[str, set[tuple[int, str]]] = field(init=False)
    regions: set[str] = field(init=False, default_factory=set)

    def __post_init__(self) -> None:

        self._macro_regions = {
            "hframework": {"fwh1", "fwh2", "fwh3", "fwh4"},
            "hcdrs": {"cdrh1", "cdrh2", "cdrh3"},
            "lframework": {"fwl1", "fwl2", "fwl3", "fwl4"},
            "lcdrs": {"cdrl1", "cdrl2", "cdrl3"},
        }
        self._macro_regions.update(
            {
                "framework": self._macro_regions["hframework"]
                | self._macro_regions["lframework"],
                "cdrs": self._macro_regions["hcdrs"] | self._macro_regions["lcdrs"],
                "vh": self._macro_regions["hcdrs"] | self._macro_regions["hframework"],
                "vl": self._macro_regions["lcdrs"] | self._macro_regions["lframework"],
            }
        )

        self._macro_regions.update(
            {"fv": self._macro_regions["vh"] | self._macro_regions["vl"]}
        )

        self.positions = {"H": set(), "L": set()}
        self.exclude = {"H": set(), "L": set()}

    def set_regions(self, regions: Union[list, str, None] = None) -> None:
        """
        Set the regions to be used. Will clear anything added using add regions.
        """
        if not regions:
            raise AssertionError(
                f"Need to specify a list of regions: {self._defined_regions}"
            )

        if isinstance(regions, str):
            regions = [regions]

        if self.not_defined:
            self.regions = self._macro_regions["fv"]
        else:
            self.regions = set()

        self.add_regions(regions)

    def add_regions(self, regions: list) -> None:
        """
        Add regions to the selection.
        """
        for region in regions:
            region = region.lower()
            if region in self._defined_regions:
                if self.not_defined:
                    self.regions = self.regions - set([region])
                else:
                    self.regions.add(region)
            elif region in self._macro_regions:
                if self.not_defined:
                    self.regions = self.regions - self._macro_regions[region]
                else:
                    self.regions = self.regions | self._macro_regions[region]
            else:
                raise AssertionError(
                    f"Got unexpected region: {region}. Allowed: {self._defined_regions} "
                )

    def add_positions(self, positions: list[tuple[int, str]], chain: str) -> None:
        for position in positions:
            self.positions[chain].add(position)

    def exclude_positions(self, positions: list[tuple[int, str]], chain: str) -> None:
        for position in positions:
            self.exclude[chain].add(position)

    def accept(self, position: tuple[int, str], chain: str) -> Optional[int]:
        if position in self.exclude[chain]:
            return None
        if (
            get_region(position, chain, self.numbering_scheme, self.definition)
            in self.regions
            or position in self.positions[chain]
        ):
            return 1
        return None


def get_region(  # pylint: disable=R0911
    position: tuple[int, str],
    chain: str,
    numbering_scheme: str = "imgt",
    definition: str = "imgt",
) -> str:
    """
    Get the region in which the position belongs given the chain, numbering scheme and definition.
    **Note** this function does not know about insertions on the sequence. Therefore, it will get the region annotation
    wrong when using non-equivalent scheme-definitions.
    To get around this please use the annotate_regions function
    which implements heuristics to get the definition correct
    in the scheme.
    """

    if numbering_scheme == "wolfguy" or definition == "wolfguy":
        raise NotImplementedError(
            "Wolguy cdr/framework identification is not implemented"
        )

    index, insertion = position
    chain = chain.upper()

    # Horrible exception cases revolving around the kabat scheme/definition and cdr h1
    # Kabat numbering scheme will be deprecated in ANARCI v1.0.0
    if definition == "kabat":
        if (
            numbering_scheme == "kabat" and chain == "H" and 31 <= index < 36
        ):  # Kabat scheme kabat definition.
            if index == 35:
                if insertion in " AB":  # Position 31 to 35B
                    return "cdrh1"

                return "fwh2"  # 31C would be framework.

            return "cdrh1"
    if numbering_scheme == "kabat":  # Kabat numbering, chothia or imgt definitions.
        if definition == "chothia" and chain == "H" and 33 <= index < 36:
            return "fwh2"
        if definition == "imgt" and chain == "H" and 34 <= index < 36:
            return "fwh2"

    try:
        return (
            _reg_one2three[
                _regions[definition][chain][
                    _index_to_imgt_state[(numbering_scheme, chain)][index]
                ]
            ]
            % chain.lower()
        )
    except KeyError:
        return "?"
