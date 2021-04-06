class Error(Exception):
    """Base class for exceptions in this module."""


class SadieIOError(Error):
    """Exception raised for SadieIO"""

    def __init__(self, msg):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return f"{self.msg}"


class IOInferError(SadieIOError):
    """Exception raised for inferring a file name that does not exist"""

    def __init__(self, msg):
        super().__init__(msg)
