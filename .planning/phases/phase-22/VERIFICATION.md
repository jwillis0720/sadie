---
status: passed
verified_at: "2026-01-25"
gaps: []
human_verification_needed: false
---

# Phase 22 Verification: Runtime Usage

## Verification Summary

- `Airr(database=...)` supports prebuilt database paths
- Database structure validation runs before annotation
- Germlines/G3 lookup is skipped when prebuilt path is provided

## Evidence

- `src/sadie/airr/igblast/germline.py` includes `validate_prebuilt_database()` and `prebuilt` logic
- `src/sadie/airr/airr.py` accepts `database` parameter and passes `prebuilt=True`
