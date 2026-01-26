# Phase 29: Add Germline Source Tracking to AIRR Output

## Goal
Add columns to AIRR output showing which germline database each gene call came from (imgt, vdjbase, ogrdb, custom).

## Background

Currently:
- `GermlineGene.source` tracks provenance at the GermlineManager level
- AIRR output only shows `v_call`, `d_call`, `j_call`, `c_call` (gene names)
- Users cannot determine which database each matched gene came from

## Tasks

### Task 1: Build Gene→Source Lookup Table
**File:** `src/sadie/airr/igblast/germline.py`

Create a cached lookup table mapping gene names to their sources:
```python
def build_source_lookup(self) -> Dict[str, str]:
    """Build gene name → source lookup table."""
    lookup = {}
    manager = GermlineManager()
    for segment in ['V', 'D', 'J', 'C']:
        for chain in ['H', 'K', 'L']:
            genes = manager.get_genes(self.species, segment, chain)
            for gene in genes:
                lookup[gene.name] = gene.source
    return lookup
```

### Task 2: Add Source Columns During Parsing
**File:** `src/sadie/airr/igblast/igblast.py`

After parsing IgBLAST output, add source lookup:
```python
def _add_source_columns(self, df: pd.DataFrame) -> pd.DataFrame:
    """Add v_call_source, d_call_source, j_call_source, c_call_source columns."""
    source_lookup = self.germline_data.get_source_lookup()
    
    for call_col in ['v_call', 'd_call', 'j_call', 'c_call']:
        source_col = f"{call_col}_source"
        df[source_col] = df[call_col].map(
            lambda x: source_lookup.get(x, np.nan) if pd.notna(x) else np.nan
        )
    return df
```

### Task 3: Handle LinkedAirrTable
**File:** `src/sadie/airr/airrtable/airrtable.py`

Ensure source columns get proper suffixes (_heavy, _light) when linking tables.

### Task 4: Add Tests
**File:** `tests/unit/airr/test_airr.py`

```python
def test_source_columns_present():
    """Verify source columns exist in AIRR output."""
    airr = Airr('human')
    result = airr.run_single(test_sequence)
    
    assert 'v_call_source' in result.columns
    assert 'd_call_source' in result.columns
    assert 'j_call_source' in result.columns
    assert 'c_call_source' in result.columns

def test_source_values_valid():
    """Verify source values are valid provider names."""
    valid_sources = {'imgt', 'vdjbase', 'ogrdb', 'custom'}
    # ... test logic
```

## Performance Considerations

- Cache source lookup table per species (avoid rebuilding)
- Consider storing lookup in metadata file during database build
- Lazy initialization to avoid overhead when sources not needed

## Success Criteria

1. ✅ AIRR output contains `v_call_source`, `d_call_source`, `j_call_source`, `c_call_source` columns
2. ✅ Source values match germline provider names
3. ✅ Source columns populated for all matched genes
4. ✅ Unmatched genes (NaN calls) have NaN sources
5. ✅ LinkedAirrTable has `_heavy`/`_light` suffixed source columns
