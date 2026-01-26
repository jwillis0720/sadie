# Phase 29 Context: Germline Source Tracking

## Decisions

### 1. Column Naming Convention
**Decision:** Follow AIRR naming patterns exactly

New columns:
- `v_call_source`
- `d_call_source`
- `j_call_source`
- `c_call_source`

### 2. Source Value Format
**Decision:** Lowercase provider names

Valid values:
- `imgt`
- `vdjbase`
- `ogrdb`
- `custom`

### 3. Edge Case Handling

| Scenario | Behavior |
|----------|----------|
| Unmatched gene (NaN call) | Use AirrTable's existing convention (NaN) |
| Custom provider used | Literal string `custom` |
| Allele exists in multiple DBs | Use priority list fallback (vdjbase > ogrdb > imgt > custom) |

### 4. LinkedAirrTable Integration
**Decision:** Follow existing AIRR suffix pattern

Columns:
- `v_call_source_heavy`, `v_call_source_light`
- `d_call_source_heavy`, `d_call_source_light`
- `j_call_source_heavy`, `j_call_source_light`
- `c_call_source_heavy`, `c_call_source_light`

---

## Implementation Notes

- Source lookup happens at GermlineData initialization (build gene→source table)
- Cache lookup table for performance
- Priority order determines source when allele exists in multiple providers
- No version information in source values (keep simple)

---
*Created: 2026-01-26*
