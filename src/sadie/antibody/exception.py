from typing import List


class Error(Exception):
    """Base class for exceptions in this module."""


class LongHCDR3Error(Error):
    """Exception raised for HCDR3 being too long for chosen numbering scheme.

    Attributes:
    """

    def __init__(
        self, sequence_name: str, hcdr3: str, chosen_scheme: str, acceptable_scheme: List[str] = ["imgt", "aho"]
    ):
        super().__init__()
        self.sequence_name = sequence_name
        self.hcdr3 = hcdr3
        self.chosen_scheme = chosen_scheme
        self.acceptable_scheme = acceptable_scheme

    def __str__(self) -> str:
        return f"{self.sequence_name} with an HCDR3 {self.hcdr3} is too long for the chosen numbering scheme {self.chosen_scheme}. Consider using {self.acceptable_scheme}"


class NumberingDecreasing(Error):

    """Exception raised for HCDR3 being too long for chosen numbering scheme.

    Attributes:
    """

    def __init__(self, sequence_name: str, msg: str):
        super().__init__()
        self.sequence_name = sequence_name
        self.msg = msg

    def __str__(self) -> str:
        return f"{self.sequence_name}:{self.msg}"


class BadAASequenceError(Error):
    """Exception raised for bad amino acid sequences.

    Attributes:
        where -- What subclass this occurs in
        position -- what position index is the bad amino acid sequence at
        aa -- The bad amino acid string
        accepted_aa -- a string of accepted aa sequences
    """

    def __init__(self, where: str, position: int, amino_acid: str, accepted_aa: str):
        super().__init__()
        self.where = where
        self.position = position
        self.amino_acid = amino_acid
        self.accepted_aa = accepted_aa

    def __str__(self) -> str:
        return "{}:Position {} is {}, only {} are accepted.".format(
            self.where, self.position, self.amino_acid, self.accepted_aa
        )


class BadAASequenceWarning(UserWarning):
    pass


class BadNTSequenceError(Error):
    """Exception raised for bad nucleotide sequences.

    Attributes:
        where -- What subclass this occurs in
        position -- what position index is the bad amino acid sequence at
        nt -- The bad nt string
        accepted_nt -- a string of accepted nt sequences
    """

    def __init__(self, where: str, position: int, nt: str, accepted_nt: str):
        super().__init__()
        self.where = where
        self.position = position
        self.nt = nt
        self.accepted_nt = accepted_nt

    def __str__(self) -> str:
        return "{}:Position {} is {}, only {} are accepted.".format(
            self.where, self.position, self.nt, self.accepted_nt
        )


class BadNTSequenceWarning(UserWarning):
    pass


class HeavyChainException(Error):
    def __init__(self, msg: str):
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


class LightChainException(Error):
    def __init__(self, msg: str):
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


class BadGene(Error):
    """Exceptions raised for bad Gene input.
    Probably due to the user asking for a v gene or j gene that does not exist for that species

    Attributes:
    """

    def __init__(self, species: str, entered_gene: str, accepted_genes: List[str]):
        super().__init__()
        self.species = species
        self.entered_gene = entered_gene
        self.accepted_genes = accepted_genes

    def __str__(self) -> str:
        difference = list(set(self.entered_gene).difference(self.accepted_genes))
        return f"{difference} doesn't exist for {self.species}. Available genes are f{self.accepted_genes}"


class AmbiguousGene(Error):
    """Exceptions raised for bad an ambiguous Gene.
    Most likely due to the user asking for a gene with multiple alleles.

    Attributes:
    """

    def __init__(self, species: str, entered_gene: str, narrow_gene: str):
        super().__init__()
        self.species = species
        self.entered_gene = entered_gene
        self.narrow_gene = narrow_gene

    def __str__(self) -> str:
        return """Entered Gene {} for {} has ambiguous listing. Try specifying the allele like {}""".format(
            self.entered_gene, self.species, self.narrow_gene
        )


class NoSpecies(Error):
    """Exceptions raised for not having a species"""

    def __init__(self, species: str, available_species: List[str]):
        super().__init__()
        self.species = species
        self.availale = available_species

    def __str__(self) -> str:
        return "Species not found {} in {}".format(self.species, self.availale)


class BadRequstedFileType(Error):
    """Exception raised for not finiding the igblast module

    Attributes:
    """

    def __init__(self, requested_type: str, accepted_types: str):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self) -> str:
        return "{} file passed, only accepts {}".format(self.requested_type, self.accepted_types)


class NumberingDuplicateIdError(Error):
    """Exception raised for having duplicated IDS"""

    def __init__(self, ids: List[str], found: int):
        super().__init__()
        self.ids = ids
        self.found = found

    def __str__(self) -> str:
        return f"Duplicate Ids are found {self.ids} {self.found} times"
