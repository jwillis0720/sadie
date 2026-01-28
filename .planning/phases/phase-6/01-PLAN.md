---
phase: 6
plan: 01
title: Test Gapped AA Fallback Translation
gap_closure: true
wave: 1
autonomous: true
---

# Plan: Test Gapped AA Fallback Translation

## Objective

Implement test T035a to verify the gapped AA fallback translation mechanism works correctly when genes have only gapped nucleotide sequences (no pre-computed gapped amino acid).

## Context

- LocalHMMBuilder._translate_gapped_nt_to_aa() provides fallback translation
- FR-013 requires fail-fast when both gapped AA and gapped NT are missing
- The fallback removes gaps, translates NT->AA, then reinserts gaps at mapped positions

## Tasks

### T035a: Implement test_gapped_aa_fallback_translation

**File**: `tests/unit/germlines/test_renumbering_integration.py`

**Tests to add**:
1. `test_translate_gapped_nt_to_aa_basic` - Test translation method directly
2. `test_gapped_aa_fallback_in_hmm_builder` - Test fallback is used in HMM building
3. `test_fallback_gap_mapping` - Verify gaps are correctly mapped from NT to AA

**Acceptance**:
- All tests pass
- Coverage of fallback translation path verified
