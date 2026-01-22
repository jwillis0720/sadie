"""
Builders - BLAST Database, Auxiliary File, Gapping, and HMM
============================================================

Utilities for building IgBLAST-compatible databases and auxiliary files.

- BlastDBBuilder: Creates BLAST databases from FASTA files
- AuxFileBuilder: Generates IgBLAST auxiliary files from gapped sequences
- GapperService: Gaps ungapped sequences to IMGT numbering using BioPython
- HMMBuilder: Generates Stockholm alignment files for HMMER
"""

from .blast import BlastDBBuilder
from .aux import AuxFileBuilder
from .gapper import GapperService, gap_sequences_batch
from .hmm import HMMBuilder, get_gapped_sequences

__all__ = [
    "BlastDBBuilder",
    "AuxFileBuilder",
    "GapperService",
    "gap_sequences_batch",
    "HMMBuilder",
    "get_gapped_sequences",
]
