"""
Germline Providers - Data Source Abstraction
=============================================

Providers abstract different germline database sources:
- IMGT: International ImMunoGeneTics information system
- OGRDB: Open Germline Receptor Database
- VDJbase: Population-specific germline alleles
- Custom: User-supplied sequences

All providers implement a common interface defined in base.py
"""

from .base import GermlineProvider
from .imgt import IMGTProvider
from .ogrdb import OGRDBProvider
from .vdjbase import VDJbaseProvider
from .custom import CustomProvider

__all__ = [
    "GermlineProvider",
    "IMGTProvider",
    "OGRDBProvider",
    "VDJbaseProvider",
    "CustomProvider",
]
