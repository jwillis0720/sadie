# Germlines Module Integration Guide

This guide shows how to update existing Sadie code to use the new germlines module instead of G3.

**IMPORTANT**: This guide is for reference. The actual integration should be done in a separate PR after the germlines module is populated with data.

## Overview

The germlines module replaces G3 API calls with local database access. Three main files need updating:

1. **IgBLAST Integration** (`src/sadie/airr/igblast/germline.py`)
2. **Reference System** (`src/sadie/reference/reference.py`)
3. **HMM/Renumbering** (`src/sadie/renumbering/aligners/hmmer.py`)

## 1. Update IgBLAST Integration

**File**: `src/sadie/airr/igblast/germline.py`

### Current Code

```python
class GermlineData:
    def __init__(self, name, receptor="Ig", database_dir=None, scheme="imgt"):
        # Currently points to old structure
        if database_dir is None:
            database_dir = Path(__file__).parent / "data" / "germlines"

        self.base_dir = database_dir
        self.v_gene_dir = database_dir / "Ig" / "blastdb" / name / f"{name}_V"
        # ...
```

### New Code

```python
class GermlineData:
    def __init__(self, name, receptor="Ig", database_dir=None, scheme="imgt"):
        # Point to germlines module
        if database_dir is None:
            germlines_base = Path(__file__).parent.parent.parent / "germlines"
            database_dir = germlines_base / "igblast"

        self.base_dir = database_dir
        # New path structure
        self.v_gene_dir = database_dir / "database" / name / f"{name}_V"
        self.d_gene_dir = database_dir / "database" / name / f"{name}_D"
        self.j_gene_dir = database_dir / "database" / name / f"{name}_J"
        self.c_gene_dir = database_dir / "database" / name / f"{name}_C"
        self.aux_path = database_dir / "aux_db" / f"{name}_gl.aux"
        self.igdata = database_dir / "internal_data"
```

### Testing

```python
# Test the update
from sadie.airr.igblast.germline import GermlineData

germline = GermlineData("human")
print(f"V database: {germline.v_gene_dir}")
print(f"Aux file: {germline.aux_path}")
print(f"V exists: {germline.v_gene_dir.with_suffix('.nhr').exists()}")
print(f"Aux exists: {germline.aux_path.exists()}")
```

## 2. Update Reference System

**File**: `src/sadie/reference/reference.py`

### Current Code

```python
class Reference:
    _endpoint = "https://g3.jordanrwillis.com/api/v1/genes"

    def __init__(self, endpoint: str = _endpoint):
        self.data = []
        self.endpoint = endpoint

    def _get_gene(self, gene: GeneEntry) -> Dict[str, str]:
        # Queries G3 API
        gene_url = url_quote(gene.gene)
        query = f"{self.endpoint}?source={gene.source}&common={gene.species}&gene={gene_url}"
        status_code, response_json = self._g3_get(query)
        return response_json[0]
```

### New Code

```python
class Reference:
    def __init__(self):
        self.data = []
        # No endpoint needed - uses local germlines

    def _get_gene(self, gene: GeneEntry) -> Dict[str, str]:
        """Get gene from germlines module instead of G3 API."""
        from sadie.germlines import get_gene_by_name

        # Fetch from germlines
        germline_gene = get_gene_by_name(gene.gene, species=gene.species)

        if not germline_gene:
            raise ValueError(f"Gene {gene.gene} not found in germlines databases")

        # Convert to old G3 format for compatibility
        response_data = {
            'gene': germline_gene.name,
            'sequence': germline_gene.sequence,
            'species': germline_gene.species,
            'segment': germline_gene.segment,
            'chain': germline_gene.chain,
            'source': germline_gene.source,
            'imgt': {
                'sequence_gapped': germline_gene.sequence_gapped,
                'sequence_gapped_aa': germline_gene.sequence_aa_gapped,
                'imgt_functional': germline_gene.functionality,
                'cdr3_aa': germline_gene.regions.get('CDR3', '') if germline_gene.regions else '',
                'fwr4_aa': germline_gene.regions.get('FWR4', '') if germline_gene.regions else '',
            }
        }

        return response_data
```

### Additional Changes

Remove G3-specific code:

```python
# Remove these methods (no longer needed):
# - endpoint property setter (validates G3 API)
# - _g3_get() method
# - G3Error exception (or rename to GermlineError)

# Keep these methods (still needed):
# - add_gene()
# - add_genes()
# - _get_genes()
```

### Testing

```python
from sadie.reference.reference import Reference
from sadie.reference.models import GeneEntry

ref = Reference()

# Test single gene
gene_entry = GeneEntry(species="human", gene="IGHV1-69*01", source="imgt")
gene_data = ref._get_gene(gene_entry)
print(f"Got gene: {gene_data['gene']}")
print(f"Source: {gene_data['source']}")
print(f"Has gapped: {'sequence_gapped' in gene_data['imgt']}")
```

## 3. Update HMM/Renumbering System

**File**: `src/sadie/renumbering/aligners/hmmer.py`

### Current Code

```python
from sadie.renumbering.clients.g3 import G3

class HMMER:
    g3 = G3()

    def get_hmm_models(self, species, chains):
        hmms = {}
        for chain in chains:
            try:
                hmm = self.g3.get_hmm(species=species, chain=chain)
                hmms[f"{species}_{chain}_imgt"] = hmm
            except Exception as e:
                # Fall back to numbering
                pass
        return hmms
```

### New Code

```python
from sadie.germlines.builders.hmm import HMMBuilder

class HMMER:
    def __init__(self):
        self.hmm_builder = HMMBuilder()
        self.numbering = NumberingTranslator()  # Keep as fallback

    def get_hmm_models(self, species, chains):
        """Get HMM models from germlines module."""
        hmms = {}

        for chain in chains:
            try:
                # Try germlines first
                hmm = self.hmm_builder.get_or_build_hmm(species, chain)
                hmms[f"{species}_{chain}_imgt"] = hmm
            except Exception as e:
                # Fall back to legacy numbering
                try:
                    hmm = self.numbering.get_hmm(species, chain)
                    if hmm:
                        hmms[f"{species}_{chain}_numbering"] = hmm
                except Exception:
                    pass

        return hmms
```

### Create HMM Builder (New File)

**File**: `src/sadie/germlines/builders/hmm.py`

```python
"""
HMM Builder for Germline Sequences
===================================

Builds HMM models from germline V and J genes.

Extracted from old G3 class but now provider-agnostic.
"""

import logging
from pathlib import Path
import pyhmmer
from typing import Optional

from ..manager import GermlineManager

logger = logging.getLogger(__name__)


class HMMBuilder:
    """Build HMM models from germline genes."""

    def __init__(self, manager: Optional[GermlineManager] = None):
        if manager is None:
            from .. import get_manager
            manager = get_manager()

        self.manager = manager
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.builder = pyhmmer.plan7.Builder(
            self.alphabet,
            architecture="hand"
        )
        self.background = pyhmmer.plan7.Background(self.alphabet)

        # Cache directory
        base = Path(__file__).parent.parent
        self.hmm_dir = base / "processed" / "hmms"
        self.hmm_dir.mkdir(parents=True, exist_ok=True)

    def get_or_build_hmm(
        self,
        species: str,
        chain: str
    ) -> pyhmmer.plan7.HMM:
        """Get cached HMM or build from germline data."""

        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"

        # Return cached if exists
        if hmm_path.exists():
            with pyhmmer.plan7.HMMFile(hmm_path) as f:
                return next(f)

        # Build new
        return self.build_hmm(species, chain)

    def build_hmm(self, species: str, chain: str) -> pyhmmer.plan7.HMM:
        """Build HMM from germline V and J genes."""

        # Get genes
        v_genes = self.manager.get_genes(species, "V", chain)
        j_genes = self.manager.get_genes(species, "J", chain)

        if not v_genes or not j_genes:
            raise ValueError(f"No genes for {species} {chain}")

        # Build Stockholm alignment
        stockholm_path = self._build_stockholm(v_genes, j_genes, species, chain)

        # Build HMM
        with pyhmmer.easel.MSAFile(
            stockholm_path,
            digital=True,
            alphabet=self.alphabet,
            format="stockholm"
        ) as msa_file:
            msa = next(msa_file)
            hmm, _, _ = self.builder.build_msa(msa, self.background)

        # Cache it
        hmm_path = self.hmm_dir / f"{species}_{chain}.hmm"
        with open(hmm_path, "wb") as f:
            hmm.write(f)

        return hmm

    def _build_stockholm(self, v_genes, j_genes, species, chain):
        """Build Stockholm alignment from V+J genes."""
        # TODO: Implement Stockholm building logic
        # (Copy from old G3 class)
        raise NotImplementedError("Stockholm building not yet implemented")
```

### Testing

```python
from sadie.renumbering.aligners.hmmer import HMMER

hmmer = HMMER()
hmms = hmmer.get_hmm_models("human", ["H", "K", "L"])

print(f"Built {len(hmms)} HMMs:")
for name in hmms:
    print(f"  {name}")
```

## 4. Remove G3 Dependencies

### Files to Remove (Eventually)

```bash
# G3 client (after ensuring everything uses germlines)
src/sadie/renumbering/clients/g3.py

# G3 tests
tests/unit/renumbering/test_g3.py
```

### Files to Update

```bash
# Update imports
src/sadie/app.py                    # Remove G3 imports
src/sadie/typing/species.py          # Update species list
tests/                                # Update tests
```

### Search for G3 References

```bash
# Find all G3 references
grep -r "from.*g3 import" src/
grep -r "G3()" src/
grep -r "g3.jordanrwillis.com" src/
```

## 5. Update Configuration Files

### reference.yml

Update to include source selection:

```yaml
# Old format
clk:
  imgt:
    human:
      - IGHV1-2*02
      - IGHJ2*01

# New format (optional - defaults work)
clk:
  sources: ["custom", "imgt", "ogrdb"]  # Priority order
  human:
    - IGHV1-2*02
    - IGHJ2*01
```

## 6. Testing Strategy

### Unit Tests

```python
# Test germlines module
pytest src/sadie/germlines/tests/

# Test integration
pytest tests/ -k "reference or igblast or hmmer"
```

### Integration Tests

```python
# Full pipeline test
from sadie.germlines import update_databases
from sadie.airr.igblast import IgBLASTN

# Update databases
update_databases("human")

# Run IgBLAST
igblast = IgBLASTN(organism="human")
results = igblast.run_file("test_sequences.fasta")

print(f"Results: {len(results)} sequences annotated")
```

### Regression Tests

Compare results with G3-based version:

```python
# Run same sequences with old (G3) and new (germlines)
# Compare gene calls, should be identical for IMGT

old_results = run_with_g3("test.fasta")
new_results = run_with_germlines("test.fasta")

assert old_results == new_results, "Results differ!"
```

## 7. Migration Checklist

- [x] Populate germlines/sources/imgt/human/ with IMGT data
- [x] Test germlines module independently
- [x] Update IgBLAST integration (germline.py)
- [x] Test IgBLAST with new paths
- [x] Update Reference system (reference.py)
- [x] Test reference building
- [x] Update HMM system (hmmer.py)
- [x] Create HMM builder (renumbering_integration.py)
- [x] Test HMM generation
- [x] Run full Sadie test suite
- [x] Compare results with G3 version
- [x] Update documentation
- [x] Mark G3 as deprecated
- [x] Create migration guide for users

## Integration Status (as of 2026-01-21)

**All integration tasks complete.** The germlines module is fully integrated with:

1. **IgBLAST Integration** - `GermlineData` class now checks `SADIE_USE_GERMLINES_MODULE` feature flag
2. **Reference System** - `Reference` class supports `use_germlines=True` parameter
3. **HMM/Renumbering** - `LocalHMMBuilder` in `renumbering_integration.py` generates HMMs from germlines

### Test Suite

Run the integration tests:
```bash
pytest tests/unit/germlines/ -v
```

All 25 tests should pass, including:
- AIRR annotation with germlines backend
- Renumbering with LocalHMMBuilder
- Reference system with G3-compatible output
- Offline operation (network disabled)

## 8. Rollback Plan

If integration has issues:

1. **Keep G3 as fallback**:
```python
try:
    from sadie.germlines import get_gene_by_name
    gene = get_gene_by_name(name, species)
except Exception:
    # Fall back to G3
    gene = fetch_from_g3(name, species)
```

2. **Feature flag**:
```python
USE_GERMLINES = os.environ.get("SADIE_USE_GERMLINES", "true").lower() == "true"

if USE_GERMLINES:
    from sadie.germlines import get_gene_by_name
else:
    from sadie.renumbering.clients.g3 import G3
```

## 9. Performance Considerations

### First Run
- Downloads/processes data: ~1-2 minutes for human
- Builds BLAST databases: ~30 seconds
- Builds HMMs: ~10 seconds

### Subsequent Runs
- Uses cached data: instant
- No network requests
- Pure local file access

### Optimization
- Pre-build databases for common species
- Include in package distribution
- Lazy loading for less common species

## Questions?

- **Module design**: See `germlines/README.md`
- **Custom sequences**: See `sources/custom/README.md`
- **Data sources**: See `sources/*/README.md`
- **Issues**: https://github.com/jwillis0720/sadie/issues
