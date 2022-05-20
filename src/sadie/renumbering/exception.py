from typing import List


class Error(Exception):
    """Base class for exceptions in this module."""


class LongHCDR3Error(Error):
    """Exception raised for HCDR3 being too long for chosen numbering scheme.

    Attributes:
    """

    def __init__(self, sequence_name, hcdr3, chosen_scheme, acceptable_scheme=("imgt", "aho")):
        super().__init__()
        self.sequence_name = sequence_name
        self.hcdr3 = hcdr3
        self.chosen_scheme = chosen_scheme
        self.acceptable_scheme = acceptable_scheme

    def __str__(self):
        return f"{self.sequence_name} with an HCDR3 {self.hcdr3} is too long for the chosen numbering scheme {self.chosen_scheme}. Consider using {self.acceptable_scheme}"


class BadRequstedFileType(Error):
    """Exception raised for not finiding the igblast module

    Attributes:
    """

    def __init__(self, requested_type: str, accepted_types: List[str]):
        super().__init__()
        self.requested_type = requested_type
        self.accepted_types = accepted_types

    def __str__(self):
        return "{} file passed, only accepts {}".format(self.requested_type, self.accepted_types)


class BadNumberingArgument(Error):
    """Exception raised for passing incorrect params to an Ancari arguments"""

    def __init__(self, passed_arguments, accepted_argumetns):
        super().__init__()
        self.passed_arguments = passed_arguments
        self.accepted_arguments = accepted_argumetns

    def __str__(self):
        return f"Passed argument {self.passed_arguments}. Only accepts {self.accepted_arguments}"


class NumberingDuplicateIdError(Error):
    """Exception raised for having duplicated IDS"""

    def __init__(self, ids, found):
        super().__init__()
        self.ids = ids
        self.found = found

    def __str__(self):
        return f"Duplicate Ids are found {self.ids} {self.found} times"


class NumberingExecutionError(Error):
    """Exception raised for hmmscan not executin"""

    def __init__(self, path, msg):
        super().__init__()
        self.path = path
        self.msg = msg

    def __str__(self):
        return f"Problem with {self.path} \n{self.msg} "
