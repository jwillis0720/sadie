# Phase 15: J Gene Matching & CDR3 Annotation Fix — Context

**Created:** 2026-01-22  
**Status:** Decisions locked

---

## Decision 1: Investigation Approach

**Decision:** Compare with G3 first, then isolated IgBLAST tests.

**Sandbox location:** `audit/` directory

**Approach:**
1. **Compare G3 vs Germlines** — Diff directory structures, BLAST DB formats, aux files, command params
2. **Isolated IgBLAST tests** — Run IgBLAST manually in `audit/` to verify fixes
3. **Iterative debugging** — Fix one issue at a time, verify with isolated test, repeat

**Sandbox files to create:**
- `audit/igblast_debug.py` — Manual IgBLAST execution script
- `audit/compare_configs.py` — G3 vs germlines configuration comparison
- `audit/test_sequences.fasta` — Small test set for quick iteration

---

## Decision 2: Root Cause Investigation Priority

**Decision:** Investigate ALL potential causes systematically.

**Investigation checklist:**

1. **J gene BLAST database**
   - [ ] Compare `germlines/igblast/database/human/human_J.*` vs G3 `airr/data/germlines/Ig/blastdb/human/human_J.*`
   - [ ] Verify BLAST DB format (makeblastdb version, sequence count)
   - [ ] Check sequence content matches

2. **Aux file**
   - [ ] Compare `germlines/igblast/aux_db/human_gl.aux` vs G3 aux file
   - [ ] Verify format (columns, separators, gene names)
   - [ ] Check all J genes are listed

3. **internal_data structure**
   - [ ] Compare `germlines/igblast/Ig/internal_data/human/` vs G3 structure
   - [ ] Verify required files present (human_V.*, human_D.*, human_J.*)
   - [ ] Check `.ndm.imgt` file format

4. **IgBLAST parameters**
   - [ ] Compare command line parameters between G3 and germlines execution
   - [ ] Check `-germline_db_J`, `-auxiliary_data`, `-ig_seqtype` params
   - [ ] Verify environment variables (IGDATA, etc.)

5. **Sequence orientation**
   - [ ] Verify input sequences are sense strand (5' to 3')
   - [ ] Check if reverse complement handling differs

---

## Decision 3: Reference Baseline

**Decision:** Use G3's working IgBLAST configuration as the gold standard.

**G3 paths (reference):**
```
src/sadie/airr/data/germlines/
├── Ig/
│   ├── blastdb/human/
│   │   ├── human_V.*
│   │   ├── human_D.*
│   │   ├── human_J.*
│   │   └── human_C.*
│   └── internal_data/human/
│       ├── human_V.*
│       ├── human_D.*
│       ├── human_J.*
│       └── human.ndm.imgt
└── aux_db/
    └── human_gl.aux
```

**Germlines paths (to fix):**
```
src/sadie/germlines/igblast/
├── database/human/
│   ├── human_V.*
│   ├── human_D.*
│   ├── human_J.*
│   └── human_C.*
├── Ig/internal_data/human/
│   ├── human_V.*
│   ├── human_D.*
│   ├── human_J.*
│   └── human.ndm.imgt
└── aux_db/
    └── human_gl.aux
```

**Comparison strategy:**
1. Diff file counts and sizes
2. Diff BLAST DB metadata (`blastdbcmd -info`)
3. Diff aux file content
4. Diff internal_data file content

---

## Decision 4: Fix Scope

**Decision:** Expect fixes in BOTH germlines module AND airr/igblast integration.

**Likely fix areas:**

**germlines module:**
- Database generation (`builders/blast.py`)
- Aux file generation (may need new builder)
- internal_data file generation
- File naming conventions

**airr/igblast integration:**
- Path resolution (`airr/igblast/germline.py`)
- IgBLAST command parameters (`airr/igblast/igblast.py`)
- Environment variable setup

**Fix approach:**
1. Identify differences via comparison
2. Fix germlines module to match G3 output format
3. Verify airr/igblast integration uses correct paths
4. Run isolated test to verify
5. Run full audit to confirm parity

---

## Deferred Ideas

*Captured during discussion, out of scope for Phase 15:*

- [ ] Automated regression tests for IgBLAST output
- [ ] Performance optimization of annotation pipeline
- [ ] Support for additional species beyond human

---

## Debug Checklist

Quick reference for investigation:

```bash
# Compare BLAST DB info
blastdbcmd -db src/sadie/airr/data/germlines/Ig/blastdb/human/human_J -info
blastdbcmd -db src/sadie/germlines/igblast/database/human/human_J -info

# Compare aux files
diff src/sadie/airr/data/germlines/aux_db/human_gl.aux \
     src/sadie/germlines/igblast/aux_db/human_gl.aux

# Compare internal_data
ls -la src/sadie/airr/data/germlines/Ig/internal_data/human/
ls -la src/sadie/germlines/igblast/Ig/internal_data/human/
```

---

*Locked: 2026-01-22*
