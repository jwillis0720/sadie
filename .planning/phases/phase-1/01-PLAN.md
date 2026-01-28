---
phase: 1
plan: 01
title: Verify Gapped Sequences for V/J Genes
gap_closure: true
wave: 1
autonomous: true
---

# Plan: Verify Gapped AA/NT Sequences

## Objective

Verify that all V and J genes in the germlines module have either gapped AA sequences (`sequence_aa_gapped`) or gapped NT sequences (`sequence_gapped`) available for HMM building.

## Context

- HMM builder requires gapped amino acid sequences for Stockholm alignment
- Two-tier approach: use gapped AA if available, fallback to translating gapped NT
- FR-013 requires fail-fast when both are missing for a gene

## Tasks

### T004a: Verify gapped sequences coverage

**Implementation:**
1. Query GermlineManager for all V and J genes across all species
2. Check each gene for `sequence_aa_gapped` or `sequence_gapped`
3. Report genes missing both
4. Validate the fallback translation works for genes with only gapped NT

**Acceptance:**
- All V/J genes have at least gapped NT available
- Genes with only gapped NT can be translated successfully
- Report showing coverage statistics

## Success Criteria

- [ ] Coverage report generated for all species
- [ ] Any gaps identified and documented
- [ ] Fallback translation verified working
