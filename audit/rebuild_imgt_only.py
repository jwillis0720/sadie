#!/usr/bin/env python
"""Rebuild germline databases using only IMGT sources for parity testing."""

import sys
sys.path.insert(0, "src")

from pathlib import Path
from sadie.germlines.pipeline import GermlinePipeline
from sadie.germlines.manager import GermlineManager

# Monkey-patch GermlineManager to use only IMGT
original_init = GermlineManager.__init__

def imgt_only_init(self, providers=None, data_dir=None):
    original_init(self, providers=["imgt"], data_dir=data_dir)

GermlineManager.__init__ = imgt_only_init

# Rebuild
print("Rebuilding germline databases with IMGT only...")
pipeline = GermlinePipeline(Path("src/sadie/germlines"))
pipeline.force_rebuild("human")
print("Done!")
