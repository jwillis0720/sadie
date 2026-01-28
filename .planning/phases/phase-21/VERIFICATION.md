---
status: passed
verified_at: "2026-01-25"
gaps: []
human_verification_needed: false
---

# Phase 21 Verification: Build CLI

## Verification Summary

- `sadie reference build` command exists with required options
- Build flow runs end-to-end with germlines backend
- Progress output and error handling present

## Evidence

- `src/sadie/app.py` defines `@reference.command("build")` with `--output` and `--use-germlines`
- Integration test run:
  ```
  PYTHONPATH=src pytest tests/integration/reference/test_reference_integration.py -k "reference and build"
  ```
- Output structure verified in `test_reference_build_with_germlines`
