class Error(Exception):
    """Base class for exceptions in this module."""


class DirectoryExistsError(FileExistsError):
    pass


class NoExtensionNameWarning(UserWarning):
    pass


class NotAValidSequenceFile(NotImplementedError):
    pass


class NotAValidCompression(NotImplementedError):
    pass
