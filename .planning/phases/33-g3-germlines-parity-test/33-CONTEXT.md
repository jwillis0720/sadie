# Phase 33 Context: G3-Germlines Parity Test

**Phase goal:** Create automated parity test validating G3 and Germlines backends produce identical AIRR output when built with the same alleles from `reference.g3.yml`.

---

## Decisions

### Database Build Strategy

| Decision | Rationale |
|----------|-----------|
| **Session-scoped fixture** | Build databases once per pytest session, reuse across all test functions |
| **No disk caching** | Build fresh each session; build step is fast enough |
| **Fail hard on G3 API error** | G3 API is still active; this is a critical parity test, not optional |
| **No --rebuild flag** | Building from scratch is the requirement; no conditional rebuild needed |

### Sequence Selection

| Decision | Rationale |
|----------|-----------|
| **Use all available sequences** | PG9_H, PG9_L, OAS_subsample_1000, catnap_nt_heavy, catnap_nt_light |
| **Human species only** | Focus validation on real sequences; other species in YAML not tested |
| **No filtering** | Include all sequences (productive and unproductive) |
| **No quick mode** | Always run full suite; speed is not a priority for this validation |

### Column Exclusions

| Decision | Rationale |
|----------|-----------|
| **Exclude only source columns** | `v_call_source`, `d_call_source`, `j_call_source`, `c_call_source` |
| **Check column presence first** | Verify both outputs have identical column sets before comparing values |
| **NaN == NaN** | Treat NaN values as equal (both backends returning NaN is valid parity) |
| **`_id` not a concern** | Not present in AIRR output; only in internal Reference structures |

### Mismatch Reporting

| Decision | Rationale |
|----------|-----------|
| **Report single cell** | Row index, column name, sequence_id, G3 value, Germlines value |
| **Include sequence_id** | Enables quick lookup; should be equal anyway, reduces code complexity |
| **Stop at first difference** | Fail hard immediately; don't accumulate mismatches |
| **Stdout only** | No diff files; print only what went wrong, avoid noise |

---

## Test Structure

```
tests/migration/
├── __init__.py
├── conftest.py          # Session-scoped fixtures for database builds
└── valid_parity.py      # Parity test functions
```

### Fixtures (conftest.py)

```python
@pytest.fixture(scope="session")
def g3_database(tmp_path_factory):
    """Build database from reference.g3.yml using G3 backend."""
    # use_germlines=False
    ...

@pytest.fixture(scope="session")
def germlines_database(tmp_path_factory):
    """Build database from reference.g3.yml using Germlines backend."""
    # use_germlines=True
    ...
```

### Test Flow

1. Build G3 database from `reference.g3.yml` (use_germlines=False)
2. Build Germlines database from `reference.g3.yml` (use_germlines=True)
3. For each FASTA file in fixtures:
   - Run `Airr(database=g3_db).run_fasta(fasta)`
   - Run `Airr(database=germlines_db).run_fasta(fasta)`
   - Compare column presence
   - Compare values (excluding source columns, NaN == NaN)
   - Fail fast on first mismatch

---

## Excluded Columns

```python
EXCLUDED_COLUMNS = [
    "v_call_source",
    "d_call_source", 
    "j_call_source",
    "c_call_source",
]
```

---

## Error Output Format

On mismatch:
```
PARITY MISMATCH DETECTED

Sequence: OAS_subsample_1000.fasta
Row: 42
Column: cdr3_aa
Sequence ID: seq_00042

G3 value:      CARDGYSSGWYFDY
Germlines value: CARDGYSSGWYFDX
```

---

## Deferred Ideas

None identified during discussion.

---

*Created: 2026-01-28*
