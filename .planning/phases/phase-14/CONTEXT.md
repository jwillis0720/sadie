# Phase 14: C Region Data Integration — Context

**Created:** 2026-01-22
**Status:** Decisions locked

---

## Decision 1: C Gene Data Sources

**Decision:** Pull C genes from primary sources, not G3.

**Sources to integrate (ALL THREE in Phase 14):**

1. **IMGT/GENE-DB**
   - URL: `https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP`
   - Contains: All C genes (IGHC, IGKC, IGLC) for all species
   - Filter: Parse FASTA headers for species="Homo sapiens" and gene type C
   - Update frequency: Regular IMGT releases
   - **File to modify:** `scripts/download_imgt.py`

2. **OGRDB**
   - URL: `https://ogrdb.airr-community.org/`
   - Contains: 105 full-length IGHC sequences (Jana et al. 2025)
   - Note: Heavy chain only (IGHC), but newer long-read sequencing data
   - **File to modify:** `scripts/download_ogrdb.py`

3. **VDJbase**
   - URL: `https://vdjbase.org/api`
   - Contains: Population-level genomically-derived C data
   - **File to modify:** `providers/vdjbase.py` (API already supports pagination)

**Rationale:** G3 is deprecated (removal after 2026-06-01). Building on primary sources ensures long-term sustainability. All three sources provide complementary C gene data.

---

## Decision 2: C Gene Scope

**Decision:** Include ALL C gene types for complete coverage.

**Heavy chain (IGH):**
- IGHA1, IGHA2 (IgA)
- IGHD (IgD - delta, not D segment)
- IGHE (IgE)
- IGHG1, IGHG2, IGHG3, IGHG4 (IgG subclasses)
- IGHM (IgM)
- IGHGP (pseudogene, if functional)

**Kappa light chain (IGK):**
- IGKC

**Lambda light chain (IGL):**
- IGLC1, IGLC2, IGLC3, IGLC4, IGLC5, IGLC6, IGLC7

**Rationale:** Complete isotype coverage needed for full antibody annotation.

---

## Decision 3: Provider Architecture

**Decision:** Extend existing IMGT provider to also fetch from GENE-DB for C genes.

**Implementation approach:**
- Update `providers/imgt.py` to handle C gene fetching
- Update `scripts/download_imgt.py` (IMGTDownloader) to fetch GENE-DB bulk file
- C genes are fetched **during the normal pipeline build** (not a separate script)
- Single source of truth: `IMGTProvider.download()` handles both V-QUEST (V/D/J) and GENE-DB (C)

**Why extend existing provider (not create new):**
- Keeps all IMGT data under one provider
- Pipeline already calls `IMGTProvider.download()`
- Simpler architecture - one provider per data source

**GENE-DB URL to add:**
```
https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP
```

**Fetching logic:**
1. Download bulk GENE-DB FASTA (all species, all genes)
2. Filter for target species (e.g., "Homo sapiens")
3. Filter for C genes (gene name matches `IG[HKL][ACDEFGM]*` pattern)
4. Write to `sources/imgt/{species}/IG{H|K|L}C.fasta`

---

## Decision 4: Parity Strategy

**Decision:** Prioritize correctness over strict G3 parity.

**Approach:**
1. First: Get C genes working with IMGT data
2. Then: Compare with G3 to understand differences
3. Accept: Results may differ from G3 if IMGT data is more current/correct

**Rationale:** G3 is being deprecated. The goal is a working germlines module, not perfect backward compatibility with a deprecated system.

---

## Decision 5: Species Scope

**Decision:** Human first, then extend.

**Phase 14:** Human C genes only
**Future:** Mouse, macaque, etc. (same infrastructure, different species filter)

---

## Deferred Ideas

*Captured during discussion, out of scope for Phase 14:*

- [ ] Automated GENE-DB update checking — separate phase
- [ ] T cell receptor (TR) constant regions — separate phase

---

## Technical Notes

**IMGT GENE-DB FASTA header format:**
```
>ACCESSION|GENE*ALLELE|SPECIES|FUNCTIONALITY|REGION|POSITIONS|LENGTH|CODON|+5'|+3'|CORRECTIONS|AA_COUNT|TOTAL_LENGTH|PARTIAL|REV_COMP
```

**Example C gene entry:**
```
>M87789|IGHG1*01|Homo sapiens|F|CH1+H+CH2+CH3+CHS|...
```

**Filter criteria for human C genes:**
- Field 3 (species) = "Homo sapiens"
- Field 2 (gene) matches pattern `IG[HKL]C*` or `IG[HKL][ADEGM]*`

---

*Locked: 2026-01-22*
