#!/usr/bin/env python3
"""
DEPRECATED: This script is deprecated. Please use regenerate_catnap_references.py instead.

The new script provides:
- Backup functionality
- Progress reporting
- Comparison of old vs new files
- Command line options
- Error handling

Usage: python scripts/regenerate_catnap_references.py
"""

import warnings

warnings.warn(
    "remake_integration.py is deprecated. Use regenerate_catnap_references.py instead.",
    DeprecationWarning,
    stacklevel=2
)

# Use this script if you ever have changes with major integration changes
from sadie.airr import Airr

input_heavy = "tests/data/fixtures/fasta_inputs/catnap_nt_heavy.fasta"
input_light = "tests/data/fixtures/fasta_inputs/catnap_nt_light.fasta"
output_heavy = "tests/data/fixtures/airr_tables/catnap_heavy_airrtable.feather"
output_light = "tests/data/fixtures/airr_tables/catnap_light_airrtable.feather"
airr_api = Airr("human", adaptable=True)
catnap_heavy = airr_api.run_fasta(input_heavy)
catnap_light = airr_api.run_fasta(input_light)
catnap_heavy.to_feather(output_heavy)
catnap_light.to_feather(output_light)
