# Phase 29 Research: Germline Source Tracking

## Summary

This phase adds provenance columns (`v_call_source`, `d_call_source`, `j_call_source`, `c_call_source`) to AIRR output, leveraging the existing `GermlineGene.source` field in the germlines module. The implementation follows established SADIE patterns for post-IgBLAST column addition.

**Primary Recommendation:** Add source lookup capability to `GermlineData` class and integrate into `Airr.run_fasta()` using the same pattern as `reference_name` and penalty columns.

---

## Standard Stack

| Component | Existing in Project | Notes |
|-----------|---------------------|-------|
| pandas | Yes (core dependency) | DataFrame operations for column addition |
| GermlineManager | Yes (`src/sadie/germlines/manager.py`) | `get_gene_by_name()` for source lookup |
| GermlineGene | Yes (`src/sadie/germlines/models.py`) | `.source` field already tracks provenance |
| AirrTable | Yes (`src/sadie/airr/airrtable/airrtable.py`) | Target for new columns |

No new dependencies required.

---

## Architecture Patterns

### Pattern 1: Post-IgBLAST Column Addition (Use This)

**Location:** `Airr.run_fasta()` in `src/sadie/airr/airr.py`

**Existing Pattern:**
```python
# After IgBLAST returns DataFrame
result.insert(2, "reference_name", pd.Series([self.name] * len(result)))
result = AirrTable(result)
result["v_penalty"] = self._v_gene_penalty
result["d_penalty"] = self._d_gene_penalty
result["j_penalty"] = self._j_gene_penalty
```

**New Pattern to Follow:**
```python
result.insert(2, "reference_name", pd.Series([self.name] * len(result)))
result = AirrTable(result)
result = self._add_source_columns(result)  # NEW
result["v_penalty"] = self._v_gene_penalty
# ... rest unchanged
```

**Confidence:** HIGH - Follows existing codebase patterns exactly.

### Pattern 2: Computed Columns in _verify() (Not Recommended for Source)

**Location:** `AirrTable._verify()` in `src/sadie/airr/airrtable/airrtable.py`

This is where columns like `v_call_top`, `v_mutation`, `liable`, `vdj_nt` are computed. These are transformations of existing columns.

**Why Not for Source:** Source lookup requires external data (GermlineManager), which `AirrTable` should not depend on. Keep AirrTable focused on DataFrame operations.

### Pattern 3: Source Lookup via GermlineManager

**Location:** `src/sadie/germlines/manager.py`

```python
# Existing API
gene = manager.get_gene_by_name("IGHV1-69*01", "human")
source = gene.source  # Returns: "imgt", "vdjbase", "ogrdb", or "custom"
```

**Challenge:** Gene names from IgBLAST may include multiple alleles (comma-separated). Must parse first allele.

---

## Don't Hand-Roll

1. **Gene→Source Mapping**: Don't build custom parsing. Use `GermlineManager.get_gene_by_name()` - it already handles priority-based lookup across providers.

2. **Provider Name Formatting**: Don't transform provider names. `GermlineGene.source` already contains lowercase strings: `"imgt"`, `"vdjbase"`, `"ogrdb"`, `"custom"`.

3. **LinkedAirrTable Suffixes**: Don't manually add `_heavy`/`_light` suffixes. The merge in `_run_scfv()` automatically applies suffixes, and `LinkedAirrTable._verify()` handles suffixed column operations.

4. **NaN Handling**: Don't special-case empty calls. Use pandas' natural NaN propagation:
   ```python
   df[call_col].map(lambda x: lookup.get(x) if pd.notna(x) else np.nan)
   ```

---

## Common Pitfalls

### Pitfall 1: Gene Name Format Mismatch
**Problem:** IgBLAST returns comma-separated alleles: `"IGHV1-69*01,IGHV1-69*02"`
**Solution:** Parse first allele before lookup:
```python
gene_name = call_value.split(",")[0] if isinstance(call_value, str) else call_value
```

### Pitfall 2: Missing Gene in Lookup
**Problem:** Gene name from IgBLAST not found in current GermlineManager configuration
**Solution:** Return `"unknown"` per CONTEXT.md decisions, or NaN if call is NaN
```python
source = lookup.get(gene_name, "unknown") if pd.notna(call_value) else np.nan
```

### Pitfall 3: Performance - Repeated Lookups
**Problem:** Calling `get_gene_by_name()` for every row is expensive
**Solution:** Build lookup dict once per species, cache on GermlineData:
```python
@functools.lru_cache(maxsize=1)
def get_source_lookup(self) -> Dict[str, str]:
    """Build gene name → source mapping."""
```

### Pitfall 4: LinkedAirrTable Column Duplication
**Problem:** Source columns might be added before merge, causing duplication
**Solution:** Add source columns in the same location as other per-table columns (in `Airr.run_fasta()` before `AirrTable()` constructor, not in `_verify()`)

### Pitfall 5: prebuilt Database Path
**Problem:** When `database=` parameter is used (prebuilt database from `sadie reference build`), GermlineManager may not have the same genes
**Solution:** For prebuilt databases, load source info from `.references_dataframe.csv.gz` if available, or fall back to "unknown"

### Pitfall 6: Custom Database Without Source Info
**Problem:** References built via Reference class may not preserve source information
**Solution:** The References dataframe includes `source` column - can be used for lookup

---

## Code Examples

### Example 1: Build Source Lookup (GermlineData)

```python
# In src/sadie/airr/igblast/germline.py

from functools import lru_cache
from typing import Dict, Optional
from sadie.germlines import GermlineManager

class GermlineData:
    # ... existing code ...
    
    @lru_cache(maxsize=1)
    def get_source_lookup(self) -> Dict[str, str]:
        """
        Build gene name → source lookup table.
        
        Returns
        -------
        Dict[str, str]
            Mapping from gene name to source provider
            (imgt, vdjbase, ogrdb, custom)
        """
        lookup: Dict[str, str] = {}
        
        # Skip for prebuilt databases - they may have different genes
        if hasattr(self, '_prebuilt') and self._prebuilt:
            return self._load_prebuilt_sources()
        
        manager = GermlineManager()
        for segment in ['V', 'D', 'J', 'C']:
            for chain in ['H', 'K', 'L']:
                try:
                    genes = manager.get_genes(
                        self.name, segment, chain, 
                        functional_only=False, strict=False
                    )
                    for gene in genes:
                        lookup[gene.name] = gene.source
                except Exception:
                    pass  # Some segment/chain combos may not exist
        
        return lookup
```

### Example 2: Add Source Columns (Airr)

```python
# In src/sadie/airr/airr.py

import numpy as np

def _add_source_columns(self, df: pd.DataFrame) -> pd.DataFrame:
    """
    Add v_call_source, d_call_source, j_call_source, c_call_source columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        AIRR DataFrame with v_call, d_call, j_call, c_call columns
        
    Returns
    -------
    pd.DataFrame
        DataFrame with source columns added
    """
    source_lookup = self.germline_data.get_source_lookup()
    
    for segment in ['v', 'd', 'j', 'c']:
        call_col = f"{segment}_call"
        source_col = f"{segment}_call_source"
        
        if call_col in df.columns:
            df[source_col] = df[call_col].apply(
                lambda x: self._lookup_source(x, source_lookup)
            )
    
    return df

def _lookup_source(
    self, call_value: Optional[str], lookup: Dict[str, str]
) -> Optional[str]:
    """
    Look up source for a gene call.
    
    Parameters
    ----------
    call_value : str or None
        Gene call, possibly comma-separated
    lookup : Dict[str, str]
        Gene name to source mapping
        
    Returns
    -------
    str or None
        Source provider or None if call is NaN
    """
    if pd.isna(call_value) or not call_value:
        return np.nan
    
    # Get first allele from comma-separated list
    first_allele = str(call_value).split(",")[0].strip()
    
    return lookup.get(first_allele, "unknown")
```

### Example 3: Integration Point in run_fasta()

```python
# In Airr.run_fasta(), after result = self.igblast.run_file(Path(file))

result.insert(2, "reference_name", pd.Series([self.name] * len(result)))

# Add source columns BEFORE converting to AirrTable
result = self._add_source_columns(result)

result = AirrTable(result)
```

### Example 4: Test for Source Columns

```python
# In tests/unit/airr/test_airr.py

def test_source_columns_in_output():
    """Verify source tracking columns are present and populated."""
    airr = Airr("human")
    # Use a known sequence that will match germline genes
    result = airr.run_single("test", PG9_SEQUENCE)
    
    # Check columns exist
    assert "v_call_source" in result.columns
    assert "d_call_source" in result.columns
    assert "j_call_source" in result.columns
    assert "c_call_source" in result.columns
    
    # Check valid values
    valid_sources = {"imgt", "vdjbase", "ogrdb", "custom", "unknown"}
    for col in ["v_call_source", "d_call_source", "j_call_source"]:
        values = result[col].dropna().unique()
        for v in values:
            assert v in valid_sources, f"Invalid source: {v}"

def test_source_nan_for_nan_calls():
    """Source should be NaN when call is NaN."""
    airr = Airr("human")
    result = airr.run_single("test", LIGHT_CHAIN_SEQUENCE)  # No D gene
    
    # D call should be NaN, so should D source
    if result["d_call"].isna().all():
        assert result["d_call_source"].isna().all()
```

---

## Sources

| Source | Confidence | Notes |
|--------|------------|-------|
| AIRR Standards 2.0 Documentation | HIGH | Custom columns allowed per spec |
| `src/sadie/germlines/models.py` | HIGH | GermlineGene.source field |
| `src/sadie/germlines/manager.py` | HIGH | get_gene_by_name() API |
| `src/sadie/airr/airr.py` | HIGH | Column addition patterns |
| `src/sadie/airr/airrtable/airrtable.py` | HIGH | _verify() and LinkedAirrTable patterns |
| CONTEXT.md decisions | HIGH | Column naming, edge cases |

---

## Verification Checklist

- [x] AIRR standard permits custom columns (verified via docs.airr-community.org)
- [x] GermlineGene.source exists and contains valid provider names (verified in models.py)
- [x] GermlineManager.get_gene_by_name() returns genes with source populated (verified in manager.py)
- [x] Column addition pattern established in Airr.run_fasta() (verified in airr.py)
- [x] LinkedAirrTable automatically handles suffixes during merge (verified in airrtable.py)
- [x] NaN handling follows pandas conventions (verified in existing code)

---

*Research completed: 2026-01-26*
