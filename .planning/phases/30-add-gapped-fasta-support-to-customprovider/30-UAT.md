---
status: complete
phase: 30-add-gapped-fasta-support-to-customprovider
source: 30-01-SUMMARY.md
started: 2026-01-27T17:01:28Z
updated: 2026-01-27T17:28:34Z
---

## Current Test

[testing complete]

## Tests

### 1. Use pre-gapped custom FASTA when present
expected: CustomProvider uses `_gapped.fasta` sequences for matching genes.
result: pass

### 2. Auto-gap fallback when gapped file missing
expected: With no `_gapped.fasta` file, CustomProvider still returns genes with gapped sequences via auto-gapping.
result: pass

### 3. Partial gapped coverage falls back per-gene
expected: Genes in `_gapped.fasta` use pre-gapped sequences while other genes auto-gap successfully.
result: pass

## Summary

total: 3
passed: 3
issues: 0
pending: 0
skipped: 0

## Gaps
