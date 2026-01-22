"""
Germline data download scripts.

This package provides scripts for downloading germline data from various sources:
- download_imgt.py: Download IMGT germline data
- download_ogrdb.py: Download OGRDB germline data from Zenodo archive

Usage:
    python -m sadie.germlines.scripts.download_ogrdb --species human
    python -m sadie.germlines.scripts.download_imgt --species human
"""

from .download_ogrdb import OGRDBDownloader

__all__ = ["OGRDBDownloader"]
