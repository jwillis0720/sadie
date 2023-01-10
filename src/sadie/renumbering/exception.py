class Error(Exception):
    """Base class for exceptions in this module."""


class BadNumberingArgument(Error):
    """Exception raised for passing incorrect params to an Ancari arguments"""

    def __init__(self, passed_arguments, accepted_argumetns):
        super().__init__()
        self.passed_arguments = passed_arguments
        self.accepted_arguments = accepted_argumetns

    # def __str__(self):
    #     return f"Passed argument {self.passed_arguments}. Only accepts {self.accepted_arguments}"


class NumberingDuplicateIdError(Error):
    """Exception raised for having duplicated IDS"""

    def __init__(self, ids, found):
        super().__init__()
        self.ids = ids
        self.found = found

    # def __str__(self):
    #     return f"Duplicate Ids are found {self.ids} {self.found} times"
