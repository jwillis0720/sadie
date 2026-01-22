# Research: Germline Database Integration

**Feature**: 002-germline-integration
**Date**: 2026-01-20
**Status**: Complete

## Overview

Research phase to understand existing integration points, G3 API format, and integration patterns before implementing germlines module connection to AIRR/renumbering systems.

## Research Questions Resolved

### R1: G3 API Response Format

**Question**: What is the exact structure of G3 API JSON responses that existing code expects?

**Decision**: Use G3 API format with nested IMGT structure

**Rationale**:
- Examined live G3 data: `/reference/data/imgt-g3.json.gz`
- Found comprehensive structure with regions (CDR1-3, FWR1-4), positions, gapped sequences
- Existing code in `reference/reference.py` expects this format

**Format Structure**:
```json
{
  "source": "imgt",
  "common": "human",
  "gene": "IGHV1-69*01",
  "sequence": "CAGGTG...",
  "imgt": {
    "sequence_gapped": "...with dots...",
    "sequence_gapped_aa": "QVQLVQ...",
    "fwr1": "...", "fwr1_aa": "...", "fwr1_start": 0, "fwr1_end": 74,
    "cdr1": "...", "cdr1_aa": "...",
    "fwr2": "...", "fwr2_aa": "...",
    "cdr2": "...", "cdr2_aa": "...",
    "fwr3": "...", "fwr3_aa": "...",
    "cdr3": "...", "cdr3_aa": "...",
    "imgt_functional": "F"
  }
}
```

**Alternatives Considered**:
- Create new format: Rejected - would require updating all Reference system consumers
- Simplify format: Rejected - loses region information needed by downstream analysis

### R2: IgBLAST Database Structure and Paths

**Question**: Where does IgBLAST expect database files and what structure does it require?

**Decision**: Use separated database structure with makeblastdb-generated files

**Rationale**:
- Current structure in `airr/data/germlines/Ig/blastdb/{species}/`
- Germlines module mirrors this in `germlines/igblast/database/{species}/`
- IgBLAST requires prefixes (e.g., `human_V`) not full file paths
- Validation via `ensure_prefix_to()` checks .nhr, .nin, .nsq files exist

**Path Structure**:
```
germlines/igblast/
├── database/{species}/
│   ├── {species}_V  (prefix for .nhr, .nin, .nsq, .ndb, etc.)
│   ├── {species}_D
│   └── {species}_J
├── aux_db/
│   └── {species}_gl.aux  (J gene CDR3 start positions)
└── internal_data/
    └── {species}.ndm.imgt  (domain classification)
```

**Alternatives Considered**:
- Single combined database: Rejected - IgBLAST requires separate V/D/J
- Different naming: Rejected - follows IgBLAST conventions

### R3: HMM Generation Workflow with pyhmmer

**Question**: How does G3 currently generate HMM models and how can we replicate locally?

**Decision**: Use pyhmmer with Stockholm alignment format, matching G3 workflow exactly

**Rationale**:
- Examined `renumbering/clients/g3.py` complete HMM building process
- G3 workflow: V/J genes → gapped AA sequences → Stockholm file → pyhmmer MSA → HMM binary
- Already have gapped sequences in germlines module from 001-germline-completion
- pyhmmer is existing dependency

**Workflow Steps**:
1. Query germlines for V and J genes (gapped AA sequences)
2. Create Stockholm alignment format:
   ```
   # STOCKHOLM 1.0
   #=GF ID species_chain
   IGHV1-69*01  QVQLVQ...
   IGHJ4*01     WGQGT...
   #=GC RF      xxx...
   //
   ```
3. Use pyhmmer to build HMM:
   ```python
   with pyhmmer.easel.MSAFile(sto_path, digital=True, alphabet=amino) as f:
       msa = next(f)
       hmm, _, _ = builder.build_msa(msa, background)
   ```
4. Cache HMM binary in `germlines/hmms/{species}_{chain}.hmm`

**Alternatives Considered**:
- Use existing G3 HMMs directly: Rejected - defeats purpose of local-first operation
- Different alignment format (FASTA): Rejected - Stockholm is pyhmmer standard
- Skip caching: Rejected - HMM building is expensive (5-10s per species/chain)

**Gapped Sequence Handling** (Clarification 2026-01-21):
The HMM builder requires gapped amino acid sequences for Stockholm alignment. The system handles this via a two-tier approach:
1. **Primary**: Use pre-computed `sequence_aa_gapped` from germlines module
2. **Fallback**: Translate `sequence_gapped` (gapped nucleotide) to gapped AA at runtime

This fallback is implemented in `LocalHMMBuilder._translate_gapped_nt_to_aa()` which:
- Removes gaps from NT sequence for translation
- Maps gap positions back to AA (3 NT gaps = 1 AA gap)
- Returns None on translation failure (logged as warning)

**Researcher Note**: Ingesting gapped AA sequences directly from source (e.g., IMGT) is preferred for accuracy. Translation fallback may introduce minor discrepancies at codon boundaries containing gaps. See FR-013 in spec.md for formal requirement.

### R4: Feature Flag Integration Pattern

**Question**: How to enable/disable germlines integration without breaking existing code?

**Decision**: Environment variable feature flag with graceful fallback

**Rationale**:
- Existing pattern in germlines module: `SADIE_USE_GERMLINES_MODULE`
- Default "true" (use germlines) for forward progress
- Set to "false" for G3 fallback during validation period
- No code changes required for users - just environment variable

**Implementation**:
```python
def use_germlines_module() -> bool:
    flag = os.getenv("SADIE_USE_GERMLINES_MODULE", "true").lower()
    return flag in ("true", "1", "yes", "on")
```

**Integration Points**:
1. IgBLAST: Check flag in `GermlineData.__init__()`, switch base_dir
2. HMM: Check flag in `HMMER.get_hmm_models()`, use LocalHMMBuilder or G3
3. Reference: Add `use_germlines` parameter (opt-in per requirement)

**Alternatives Considered**:
- Function parameter everywhere: Rejected - too invasive, breaks existing code
- Config file: Rejected - adds complexity, environment variable simpler
- Automatic detection: Rejected - implicit behavior causes confusion

### R5: Testing Strategy for Integration Validation

**Question**: How to validate germlines integration produces equivalent results to G3?

**Decision**: Mirror critical path tests in separate directory with feature flag control

**Rationale**:
- Spec requires mirrored tests (FR-007, FR-008)
- Existing test suites in `tests/unit/airr/` and `tests/unit/renumbering/`
- Create `tests/unit/germlines/` to isolate integration tests
- Use same test data, compare outputs between G3 and germlines backends

**Test Structure**:
```python
# tests/unit/germlines/test_airr_integration.py
def test_airr_annotation_with_germlines(monkeypatch):
    # Enable germlines
    monkeypatch.setenv("SADIE_USE_GERMLINES_MODULE", "true")

    # Run AIRR annotation
    result_germlines = run_airr_annotation(test_sequences)

    # Compare with G3 baseline (if available)
    # Or validate structure/completeness
    assert result_germlines.genes_called == expected_genes
```

**Coverage Requirements** (from spec):
- AIRR annotation critical paths
- Renumbering HMM alignment critical paths
- Provider selection logic
- Offline operation
- Error handling (missing data, invalid sequences)

**Alternatives Considered**:
- Modify existing tests: Rejected - could break CI, hard to toggle
- Property-based testing: Rejected - overkill for integration validation
- Manual testing only: Rejected - spec explicitly requires automated tests

## Technology Stack Summary

**Confirmed Technologies**:
- Python 3.10+ (existing SADIE requirement)
- pyhmmer 0.10+ (existing dependency for HMM building)
- Biopython (existing dependency for sequence parsing)
- pydantic (existing dependency in germlines module)
- pytest (existing test framework)

**New Dependencies**: None (all leverage existing)

**File Formats**:
- FASTA (germline sequences)
- Stockholm (MSA format for HMMs)
- BLAST databases (.nhr, .nin, .nsq, etc.)
- HMM binaries (pyhmmer format)
- JSON/Dict (G3 API compatibility)

## Integration Pattern Summary

**Adapter Pattern**: GermlineGene → G3 format transformation
- Maintains compatibility with existing Reference system
- Isolated in `g3_adapter.py` module
- Handles region extraction and formatting

**Builder Pattern**: LocalHMMBuilder for renumbering
- Queries germlines for V/J genes
- Builds Stockholm alignments
- Generates HMM binaries via pyhmmer
- Caches results for performance

**Strategy Pattern**: Feature flag switching
- Runtime selection of germlines vs G3 backend
- Transparent to calling code
- Graceful fallback on errors

## Risks and Mitigations

**Risk 1**: Germlines data quality differs from G3
- **Mitigation**: Validation tests compare results, gapping algorithm matches IMGT standard

**Risk 2**: Performance regression
- **Mitigation**: HMM caching, BLAST database indexing, lazy initialization

**Risk 3**: Breaking existing workflows
- **Mitigation**: Feature flags, G3 remains default, comprehensive backwards compat testing

**Risk 4**: Missing species/chains in germlines
- **Mitigation**: Clear error messages, fallback to G3 if available, documentation

## Next Steps

Phase 0 complete. Proceed to Phase 1:
1. Create data-model.md (entity definitions)
2. Create contracts/ (API interfaces)
3. Create quickstart.md (integration guide)
4. Update agent context

