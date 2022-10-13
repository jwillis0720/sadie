from collections.abc import Generator
from typing import Any

from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

from .Interfaces import SequenceIterator as SequenceIterator

class AbiIterator(SequenceIterator):
    trim: Any
    def __init__(self, source, trim: bool = ...) -> None: ...
    def parse(self, handle): ...
    def iterate(self, handle) -> Generator[Any, None, None]: ...
