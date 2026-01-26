---
status: passed
verified_at: "2026-01-25"
gaps: []
human_verification_needed: false
---

# Phase 20 Verification: Integration Foundation

## Verification Summary

- `References.from_yaml(..., use_germlines=True)` routes to germlines backend
- G3 adapter generates deterministic `_id` values
- Explicit provider routing uses `GermlineManager(providers=[source])`

## Evidence

- `src/sadie/reference/reference.py` includes `use_germlines` parameter in `from_yaml`
- `src/sadie/reference/reference.py` uses `GermlineManager(providers=[source])` in `_get_gene`/`_get_genes`
- `src/sadie/germlines/g3_adapter.py` includes `_generate_id()` and `_id` in output
- Tests in `tests/unit/reference/test_reference.py`:
  - `test_reference_use_germlines`
  - `test_references_from_yaml_use_germlines`
