# Feature Specification: Germlines Module Completion

**Feature Branch**: `001-germline-completion`
**Created**: 2026-01-08
**Status**: Draft
**Input**: User description: "Complete germlines module: add VDJbase provider, populate data, integrate with Sadie, remove G3 dependency"

## Clarifications

### Session 2026-01-08

- Q: VDJbase Data Source & Format - Should VDJbase provider use API automation or manual FASTA approach? → A: API automation implemented with `VDJbaseProvider.download()` method; also supports manual FASTA files with documented download instructions
- Q: Auto-Gapping Implementation Method - What tool/approach should be used for IMGT gapping of ungapped sequences? → A: Use Biopython alignment against IMGT-gapped templates (per-gene when available; fallback to per-segment consensus)
- Q: G3 Migration Strategy & Timeline - Should G3 be removed immediately or phased out gradually? → A: Feature flag only (no automatic fallback) during validation, then deprecation after validation period
- Q: Performance Monitoring & Observability - What level of logging/metrics should be implemented? → A: Basic logging (INFO/WARNING/ERROR) with timing metrics for rebuild operations
- Q: Test Data Management - How should test data be managed for CI/CD and unit tests? → A: Small curated test dataset checked into repo (5-10 genes per segment, mock data)

### Session 2026-01-14

- Q: Validation fallback behavior - Should validation use feature-flag only or automatic fallback? → A: Feature-flag only; no automatic fallback
- Q: Duplicate gene names - How to handle same gene name with different sequences? → A: Keep highest-priority sequence; drop lower-priority duplicates with warning
- Q: Auto-gapping template strategy - Which IMGT-gapped templates should alignment use? → A: Per-gene template when available; fallback to per-segment consensus

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Add Custom Germline Sequences (Priority: P1)

A researcher discovers a novel IGHV allele through sequencing and wants to add it to their local germline database for annotation of new sequences without waiting for IMGT/OGRDB approval.

**Why this priority**: Core value proposition - enables immediate use of novel germlines. This is the primary reason for decoupling from G3. Can be tested independently of other data sources.

**Independent Test**: User can add a FASTA file to `src/sadie/germlines/sources/custom/human/IGHV.fasta`, run Sadie annotation, and see the custom allele used in results. Delivers immediate value for researchers with novel alleles.

**Acceptance Scenarios**:

1. **Given** user has a novel IGHV sequence, **When** they add it to `src/sadie/germlines/sources/custom/human/IGHV.fasta` and run Sadie, **Then** the pipeline auto-detects the change, rebuilds databases, and annotates sequences using the custom allele
2. **Given** custom allele has same name as IMGT allele, **When** Sadie runs annotation, **Then** custom version takes priority per constitution principle II
3. **Given** custom FASTA has invalid nucleotides, **When** pipeline loads it, **Then** system logs warning with specific error but continues with valid sequences

---

### User Story 2 - Use Local Germline Databases Offline (Priority: P1)

A bioinformatician working in a secure/air-gapped environment needs to annotate antibody sequences using IMGT/OGRDB data without internet access.

**Why this priority**: Core requirement for local-first operation (Constitution Principle III). Must work for feature to meet success criteria. No dependency on external APIs.

**Independent Test**: Populate IMGT/OGRDB data once, disconnect from network, run full Sadie pipeline. Should complete successfully without any network calls.

**Acceptance Scenarios**:

1. **Given** IMGT and OGRDB data populated in src/sadie/germlines/sources/, **When** user runs Sadie annotation without internet, **Then** annotation completes successfully using local data
2. **Given** src/sadie/germlines/sources/ directories are empty, **When** user attempts annotation, **Then** system provides clear error with instructions to populate data (references download scripts)
3. **Given** data was populated 6 months ago, **When** user runs annotation, **Then** pipeline uses cached data without checking for updates (offline-first)

---

### User Story 3 - Priority-Based Database Selection (Priority: P2)

A researcher wants to use OGRDB alleles preferentially over IMGT when available, with fallback to IMGT for genes not in OGRDB.

**Why this priority**: Demonstrates priority system working correctly. Valuable for researchers who trust OGRDB curation but need IMGT completeness. Tests deduplication logic.

**Independent Test**: Configure manager with priority `["ogrdb", "imgt"]`, run annotation, verify OGRDB alleles used when present and IMGT used for genes absent from OGRDB.

**Note**: This tests user-configured priority override. Default system priority is `custom > ogrdb > vdjbase > imgt`.

**Acceptance Scenarios**:

1. **Given** both OGRDB and IMGT have IGHV1-69*01, **When** priority is `["ogrdb", "imgt"]`, **Then** OGRDB version is used for annotation
2. **Given** IMGT has gene X but OGRDB does not, **When** priority is `["ogrdb", "imgt"]`, **Then** IMGT version of gene X is included in merged database
3. **Given** two providers have exact same sequence but different names, **When** merging, **Then** both names are kept (novel genes principle)

---

### User Story 4 - Integrate with Existing Sadie Workflows (Priority: P1)

An existing Sadie user runs their standard IgBLAST annotation workflow and wants it to transparently use the new germlines module without changing their scripts.

**Why this priority**: Backward compatibility is mandatory per Constitution Principle V. If this fails, existing users' workflows break. Must maintain API compatibility.

**Independent Test**: Run existing Sadie IgBLAST annotation test suite against new germlines module. All tests should pass without modification.

**Acceptance Scenarios**:

1. **Given** existing Sadie code calls `GermlineData("human")`, **When** germlines module is integrated, **Then** IgBLAST paths resolve correctly to new database locations
2. **Given** Reference system requests gene "IGHV1-69*01", **When** using new germlines module, **Then** gene data is returned in same format as G3 API response
3. **Given** HMM builder needs V gene sequences, **When** calling germlines module, **Then** gapped sequences are available for Stockholm alignment building
4. **Given** feature flag `SADIE_USE_GERMLINES_MODULE=false`, **When** user runs annotation, **Then** system falls back to G3 API without errors (validation period only)

---

### User Story 5 - Add VDJbase Provider (Priority: P2)

A researcher wants to use VDJbase genotype data for population-specific germline analysis.

**Why this priority**: Completes the ranked database list requirement. Adds value but not critical for MVP. Can be added after core providers work.

**Independent Test**: Implement VDJbase provider, populate with sample data, configure priority as `["vdjbase", "imgt"]`, verify VDJbase sequences used in annotation.

**Acceptance Scenarios**:

1. **Given** VDJbase FASTA files in src/sadie/germlines/sources/vdjbase/human/, **When** manager loads genes, **Then** VDJbase sequences are parsed and available
2. **Given** VDJbase has population-specific allele, **When** priority includes vdjbase first, **Then** that allele is used instead of IMGT reference
3. **Given** VDJbase API changes format, **When** provider attempts to parse, **Then** clear error message explains format issue with link to documentation

---

### User Story 6 - Populate Reference Data Sources (Priority: P1)

A new Sadie user wants to set up germlines module with standard IMGT and OGRDB data.

**Why this priority**: Without data, module cannot function. First-time setup experience is critical. Must be straightforward and well-documented.

**Independent Test**: Fresh Sadie install, user follows README instructions to download IMGT/OGRDB data, runs validation script, receives confirmation data is ready.

**Acceptance Scenarios**:

1. **Given** empty src/sadie/germlines/sources/ directories, **When** user runs `src/sadie/germlines/scripts/download_imgt.py human`, **Then** IMGT FASTA files are downloaded and validated with INFO-level logging showing progress
2. **Given** IMGT data is downloaded, **When** user runs first annotation, **Then** pipeline builds BLAST databases automatically (~1-2 minutes) with timing metrics logged
3. **Given** download script fails partway, **When** user re-runs it, **Then** script resumes from checkpoint and doesn't re-download completed files, logging resume event

---

### Edge Cases

- **What happens when src/sadie/germlines/sources/ files change mid-annotation?** Pipeline checks hashes at start; if files change during run, current run completes with old data, next run triggers rebuild
- **How does system handle corrupted BLAST databases?** Validation detects corruption via makeblastdb exit codes; automatically triggers rebuild with clear error message
- **What if two providers have sequences that are 99% similar but different names?** Both included (not considered duplicates); deduplication only applies to identical gene names across providers
- **What if two providers have the same gene name but different sequences?** Keep highest-priority provider's sequence; drop lower-priority duplicates and log a warning
- **How to handle mixed gapped/ungapped sources?** Auto-detect gaps (contains "." or "-"); convert ungapped to gapped using Biopython alignment against IMGT-gapped templates (per-gene fallback to per-segment consensus) for V/J segments
- **What if user has 10GB of custom sequences?** Pipeline warns about database size; BLAST building may take longer; consider optimization in future (chunking, indexing)
- **VDJbase API deprecation?** Provider gracefully degrades; if API unavailable, uses cached data with warning; user can manually download FASTAs

## Requirements *(mandatory)*

### Functional Requirements

#### VDJbase Provider
- **FR-001**: System MUST implement VDJbaseProvider following base provider interface (GermlineProvider)
- **FR-002**: VDJbaseProvider MUST parse VDJbase FASTA format and create GermlineGene objects
- **FR-002a**: VDJbase FASTA format MUST follow standard FASTA specification: header line starting with ">", followed by gene identifier, optional metadata (population, genotype) separated by "|", and sequence on subsequent lines
- **FR-002b**: VDJbase FASTA header format MUST be: `>{gene_name}|{species}|{segment}|{chain}[|population={pop}][|genotype={gt}]` where fields in brackets are optional
- **FR-002c**: VDJbaseProvider MUST extract population and genotype metadata from FASTA headers when present and store in GermlineGene.metadata dictionary
- **FR-003**: VDJbaseProvider MUST handle manually-added FASTA files in src/sadie/germlines/sources/vdjbase/ directory
- **FR-004**: System MUST support VDJbase in priority ordering (e.g., `["ogrdb", "vdjbase", "imgt"]`)
- **FR-004b**: When multiple providers supply the same gene name with different sequences, system MUST keep the highest-priority provider's sequence, drop lower-priority duplicates, and log a WARNING with gene name and provider sources
- **FR-004a**: When VDJbase data is unavailable or parsing fails, VDJbaseProvider MUST log WARNING with message "VDJbase data not found at {path}. Skipping VDJbase provider." and continue without error

#### Data Population
- **FR-005**: System MUST provide download scripts for IMGT data (`src/sadie/germlines/scripts/download_imgt.py`)
- **FR-006**: System MUST provide download scripts for OGRDB data (`src/sadie/germlines/scripts/download_ogrdb.py`)
- **FR-006a**: OGRDB download script MUST fetch data from Zenodo archive (https://zenodo.org/records/18145568/files/ogrdb_archive.tgz) and extract sequences from SQL dump's `gene_description` table, outputting both ungapped (`IG{chain}{segment}.fasta`) and gapped (`IG{chain}{segment}_gapped.fasta`) FASTA files
- **FR-007**: VDJbase automated download is implemented via `VDJbaseProvider.download()` method; manual download instructions also available per FR-010
- **FR-008**: Download scripts MUST support species parameter (e.g., `--species human mouse`)
- **FR-009**: Download scripts MUST validate downloaded data (valid FASTA, expected file structure)
- **FR-009a**: Valid FASTA validation MUST check: (1) file contains at least one header line starting with ">", (2) sequences contain only valid nucleotides (ACGT) or IUPAC ambiguity codes (RYSWKMBDHVN), (3) no empty sequences, (4) headers contain parseable gene identifiers
- **FR-009b**: Expected file structure validation MUST verify: (1) species subdirectory exists (e.g., src/sadie/germlines/sources/imgt/human/), (2) segment files exist (IGHV.fasta, IGHD.fasta, IGHJ.fasta minimum for heavy chain), (3) file permissions are readable
- **FR-009c**: Download scripts MUST implement resume capability using checkpoint file (.download_progress.json) storing: last_downloaded_file, total_files, completed_files, timestamp
- **FR-009d**: Download scripts MUST log progress at INFO level every 10 files: "Downloaded {completed}/{total} files ({percentage}%)"
- **FR-010**: System MUST document manual download instructions in src/sadie/germlines/sources/*/README.md for each provider

#### Sadie Integration
- **FR-011**: System MUST update `src/sadie/airr/igblast/germline.py` to use new germlines module paths
- **FR-011a**: IgBLAST integration MUST resolve database paths via GermlineManager.get_blast_db_path(species, segment) returning absolute path to .ndb/.nsq/.nin files
- **FR-012**: System MUST update `src/sadie/reference/reference.py` to query germlines module instead of G3 API
- **FR-012a**: G3 API response format MUST be replicated as: `{"gene": str, "sequence": str, "sequence_gapped": str, "species": str, "segment": str, "chain": str, "source": str, "functional": bool, "regions": {"fwr1": str, "cdr1": str, "fwr2": str, "cdr2": str, "fwr3": str}}`
- **FR-012b**: Reference system adapter MUST call GermlineManager.get_gene_by_name(name, species) and transform GermlineGene object to G3 response format structure
- **FR-013**: System MUST implement HMM builder in germlines module for renumbering system
- **FR-013a**: HMM builder MUST accept gapped V/J sequences and generate Stockholm format alignment files for HMMER hmmbuild
- **FR-013b**: HMM builder MUST call GermlineManager.get_gapped_sequences(species, segment) returning List[Tuple[gene_name, gapped_sequence]]
- **FR-013c**: Stockholm alignment MUST include minimum 3 sequences per segment/chain combination. Sequences selected by: (1) all functional alleles from highest-priority provider, (2) if <3, add from next provider in priority order. Maximum 100 sequences per alignment to limit file size.
- **FR-014**: System MUST maintain existing API signatures for GermlineData class (backward compatibility)
- **FR-014a**: GermlineData class MUST support legacy methods: get_v_genes(species), get_d_genes(species), get_j_genes(species) returning List[str] of gene names
- **FR-015**: System MUST return gene data in G3-compatible format from Reference system for existing code

#### G3 Dependency Removal
- **FR-016**: System MUST implement feature flag `SADIE_USE_GERMLINES_MODULE` as environment variable (default: "true")
- **FR-016a**: When SADIE_USE_GERMLINES_MODULE="true", system MUST use germlines module exclusively (no G3 API calls)
- **FR-016b**: When SADIE_USE_GERMLINES_MODULE="false", system MUST use G3 API exclusively (for rollback capability during validation period only)
- **FR-017**: System MUST NOT implement automatic fallback; validation uses explicit feature flag toggle only
- **FR-017a**: Validation period duration MUST be: minimum 30 days of production use OR 100% test pass rate on 3 consecutive releases, whichever comes first
- **FR-017b**: Validation period success criteria: (1) SC-004 maintained (100% test pass), (2) zero critical bugs filed against germlines module, (3) performance metrics within 10% of G3 baseline
- **FR-018**: System MUST maintain G3 API response format compatibility for Reference system (adapter pattern)
- **FR-019**: After validation period and deprecation notice, system MUST remove G3 dependencies including `src/sadie/renumbering/clients/g3.py` and all g3.jordanrwillis.com API calls
- **FR-019a**: Deprecation notice MUST be issued minimum 60 days before G3 code removal, communicated via: (1) CHANGELOG entry, (2) runtime WARNING log when SADIE_USE_GERMLINES_MODULE="false", (3) GitHub issue/discussion
- **FR-019b**: G3 removal scope MUST include: src/sadie/renumbering/clients/g3.py, all imports of G3 class, SADIE_USE_GERMLINES_MODULE feature flag code, G3-related tests

#### Auto-Gapping
- **FR-020**: System MUST automatically detect ungapped sequences (no "." or "-" characters)
- **FR-021**: System MUST gap V and J segment sequences using IMGT numbering scheme via Biopython pairwise alignment against IMGT-gapped templates
- **FR-021a**: Template selection MUST use per-gene IMGT-gapped reference from IMGT provider FASTA when available; otherwise derive per-segment consensus from IMGT-gapped sequences for the same segment/chain
- **FR-021b**: Gapper MUST translate nucleotide sequence to amino acid before alignment, then map gaps back to nucleotide sequence using codon positions
- **FR-021c**: If alignment fails, system MUST log WARNING "Failed to gap sequence {gene_name}: {error}" and store ungapped version only
- **FR-021d**: Gap characters MUST use "." (period) for consistency with IMGT gapped format, inserted at codon boundaries in nucleotide sequences
- **FR-021e**: If no IMGT-gapped template exists for a gene (e.g., novel custom allele), system MUST use the per-segment consensus template from IMGT provider and log WARNING "No gene-specific template for {gene_name}, using {segment} consensus template"
- **FR-022**: System MUST store both gapped and ungapped versions in normalized/
- **FR-022a**: Normalized directory structure MUST be: normalized/{species}/gapped/IG{chain}{segment}.fasta and normalized/{species}/ungapped/IG{chain}{segment}.fasta
- **FR-023**: D segments MUST be stored under normalized/{species}/ungapped/ only (no gapping required)
- **FR-024**: Gapping implementation MUST reuse existing germlines pipeline utilities and Biopython to avoid duplication

#### Testing & Validation
- **FR-025**: System MUST provide regression tests comparing germlines module output vs G3 output
- **FR-025a**: Regression test pass criteria: sequence identity ≥99.9% for matched genes, same functional status, CDR/FWR boundaries within ±2 positions
- **FR-025b**: Standard test sequences MUST include: IGHV1-69*01, IGHV3-23*01, IGHD3-3*01, IGHJ4*01 (minimum representative set)
- **FR-026**: System MUST provide unit tests for each provider in isolation using curated test dataset
- **FR-027**: System MUST provide integration tests for full pipeline with mock data
- **FR-027a**: Mock data structure MUST include: src/sadie/germlines/tests/data/{provider}/{species}/{segment}.fasta with 2-3 genes per segment for fast CI/CD execution
- **FR-028**: System MUST validate IgBLAST compatibility with standard test sequences
- **FR-028a**: IgBLAST compatibility test MUST verify: (1) makeblastdb runs successfully, (2) blastn search returns hits, (3) alignment scores match expected values
- **FR-029**: System MUST include curated test dataset (5-10 genes per segment) checked into repository under `src/sadie/germlines/tests/data/`
- **FR-030**: Test data MUST include examples from each provider (custom, IMGT, OGRDB, VDJbase) with both gapped and ungapped variants
- **FR-031**: Test data MUST include edge cases: invalid sequences, duplicate names with different sequences, identical sequences with different names

#### IgBLAST Database Requirements
- **FR-037**: System MUST generate IgBLAST auxiliary files for V and J segments containing CDR/FWR region annotations
- **FR-037a**: Auxiliary file format MUST follow IgBLAST specification: tab-separated with columns: gene_name, FWR1_start, FWR1_end, CDR1_start, CDR1_end, FWR2_start, FWR2_end, CDR2_start, CDR2_end, FWR3_start, FWR3_end
- **FR-037b**: CDR/FWR boundaries MUST be determined from IMGT-gapped sequences aligned to IMGT templates, using IMGT positions: FWR1(1-26), CDR1(27-38), FWR2(39-55), CDR2(56-65), FWR3(66-104)
- **FR-037c**: If CDR/FWR boundaries cannot be determined for a gene, auxiliary file MUST omit that gene with WARNING log
- **FR-038**: System MUST generate makeblastdb command with parameters: `-dbtype nucl -parse_seqids -hash_index` for nucleotide databases
- **FR-038a**: BLAST database files MUST be stored in igblast/database/{species}/ with naming convention: {species}_{segment}.nhr, .nin, .nsq (input FASTA: {species}_{segment}.fasta)
- **FR-039**: System MUST generate IgBLAST internal_data files: organism.yaml containing species metadata and available segments
- **FR-039a**: organism.yaml format MUST be:
  ```yaml
  {species}:
    display_name: "{Species Name}"  # e.g., "Homo sapiens"
    chains:
      - name: "IGH"
        segments: [V, D, J]
        path: "igblast/database/{species}/"
      - name: "IGK"
        segments: [V, J]
        path: "igblast/database/{species}/"
      - name: "IGL"
        segments: [V, J]
        path: "igblast/database/{species}/"
    version: "{data_version}"  # ISO date of last rebuild (YYYY-MM-DD)
  ```

#### Observability & Logging
- **FR-032**: System MUST implement structured logging at INFO, WARNING, and ERROR levels using Python logging module
- **FR-032a**: Log format MUST be: `{timestamp} - {logger_name} - {level} - {message}` using ISO 8601 timestamp format (YYYY-MM-DDTHH:MM:SS.mmmZ)
- **FR-032b**: Structured logging MUST use key-value pairs for machine-readable fields: `operation={op_name} duration_ms={duration} status={success|failure} provider={name}`
- **FR-033**: System MUST log timing metrics for key operations: data download duration, database rebuild duration, validation time
- **FR-033a**: Timing metrics MUST be logged at INFO level with millisecond precision: "Operation {operation_name} completed in {duration_ms}ms"
- **FR-034**: System MUST log provider load events with gene counts (e.g., "Loaded 450 IGHV genes from IMGT provider")
- **FR-035**: System MUST log change detection events indicating which sources triggered rebuild
- **FR-035a**: Change detection log MUST include: file path, hash value (first 8 chars), change type (new|modified|deleted)
- **FR-036**: System MUST provide clear error messages with remediation guidance for common failure modes
- **FR-036a**: Error message format MUST be: "{Error description}. {Cause}. {Remediation: specific steps to resolve}" with max 200 characters
- **FR-036b**: Common failure mode error messages MUST include: (1) Missing data: "IMGT data not found at {path}. Cause: src/sadie/germlines/sources/ directory empty. Remediation: Run src/sadie/germlines/scripts/download_imgt.py human", (2) Invalid FASTA: "Invalid sequence in {file}:{line}. Cause: Contains non-ACGT characters. Remediation: Check sequence contains only ACGT", (3) BLAST build fail: "makeblastdb failed with exit code {code}. Cause: {stderr}. Remediation: Check disk space and file permissions"

### Key Entities

- **VDJbaseProvider**: Provider implementation for VDJbase database. Inherits from GermlineProvider. Handles VDJbase-specific FASTA parsing and metadata extraction. Reads manually-added FASTA files from src/sadie/germlines/sources/vdjbase/ directory. API automation is deferred as optional future enhancement.

- **GermlineGene**: Unified data model (already exists). Must be compatible with VDJbase data fields including population-specific metadata if available.

- **Download Scripts**: Python scripts for automated data fetching. IMGT and OGRDB have automated download scripts. VDJbase uses documented manual download instructions. Each script supports species filtering, resume capability, validation.

- **Integration Adapters**: Code in existing Sadie modules that bridges old G3 API calls to new germlines module. Includes path mapping, response format conversion, backward compatibility shims.

- **Test Dataset**: Curated collection of germline sequences (5-10 genes per segment) stored in `src/sadie/germlines/tests/data/`. Includes representative examples from all providers with both gapped/ungapped variants. Covers edge cases for comprehensive unit testing without requiring full reference data download.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: User can add custom germline FASTA and see it used in annotation within 5 minutes (including auto-rebuild time)
- **SC-002**: Full Sadie pipeline completes successfully without network access after initial data population
- **SC-003**: Download scripts successfully populate IMGT and OGRDB data for human species in under 10 minutes (excluding network speed); VDJbase setup via documented manual instructions
- **SC-004**: Existing Sadie IgBLAST test suite passes 100% with germlines module (zero regression)
- **SC-005**: Germlines module uses <500MB disk space for human reference data (IMGT + OGRDB + VDJbase)
- **SC-006**: Priority ordering produces different gene selections: `["ogrdb", "imgt"]` vs `["imgt", "ogrdb"]` yields measurably different results when providers overlap
- **SC-007**: Feature flag successfully controls germlines vs G3 usage; germlines module is default and fully functional; G3 code isolated behind feature flag for validation period
- **SC-008**: Documentation provides step-by-step setup instructions that new user can complete in under 30 minutes (measured from git clone to first successful annotation run, including download script execution and validation)
- **SC-009**: Change detection correctly triggers rebuild only when src/sadie/germlines/sources/ files change (verified via timestamp checks)
- **SC-010**: System provides clear, actionable error messages for 100% of common failure modes (missing data, corrupt files, invalid sequences) with structured logging at appropriate levels
- **SC-011**: Timing metrics are logged for all major operations (download, rebuild, validation) allowing performance monitoring
- **SC-012**: Unit and integration tests run successfully in CI/CD using curated test dataset without requiring full reference data download (tests complete in <5 minutes)

### Qualitative Measures

- **SC-013**: Users report setup process is straightforward (validation via user testing with 3+ scientists)
- **SC-014**: Germlines module can be extracted as standalone package with <4 hours effort (measured by actual extraction attempt)
- **SC-015**: Code reviewers confirm alignment with Constitution principles (all 5 core principles verified in PR review)

## Assumptions

1. **IMGT data license**: Assumed IMGT reference data can be downloaded and used locally per their terms of service. User responsibility to comply with licensing.

2. **VDJbase data availability**: Assumed VDJbase provides downloadable FASTA files via their website or data portal. Users follow documented manual download instructions. API automation deferred to reduce external dependency risk.

3. **Sequence format consistency**: Assumed IMGT uses IMGT-gapped format, OGRDB may be ungapped or gapped, VDJbase format follows standard FASTA. Auto-detection handles mixed formats.

4. **Species scope**: Initial implementation focuses on human. Mouse and other species follow same pattern but may require additional validation.

5. **IgBLAST version**: Assumed IgBLAST version bundled with Sadie is compatible with generated BLAST databases. No version-specific issues expected.

6. **Migration strategy**: Phased migration approach with feature flag for controlled rollout. Germlines module is primary (default); G3 remains available only via explicit feature flag toggle during validation period (no automatic fallback). Validation period duration: minimum 30 days of production use OR 100% test pass rate on 3 consecutive releases, whichever comes first. Success criteria: (1) SC-004 maintained (100% test pass), (2) zero critical bugs, (3) performance within 10% of G3 baseline. After successful validation and 60-day deprecation notice (CHANGELOG + runtime WARNING + GitHub discussion), G3 dependencies removed completely.

7. **Performance requirements**: Assumed ~1-2 minute rebuild time is acceptable for interactive use. Batch processing users can pre-build databases.

8. **Test data scope**: Curated test dataset (5-10 genes per segment) sufficient for unit and integration testing. CI/CD does not require full reference data download. Full-scale testing with complete reference data performed manually or in dedicated test environment.

## Out of Scope

- **Multi-species simultaneous support**: Initial implementation handles one species at a time. Multi-species in single run is future enhancement.

- **Automatic update checking**: System does not check for new IMGT/OGRDB releases. User manually re-runs download scripts to update.

- **GUI for database management**: Command-line interface only. GUI is future enhancement if user demand exists.

- **Cloud storage integration**: All data stored locally. S3/cloud backup is user responsibility via standard tools.

- **Germline alignment visualization**: System provides sequences; visualization tools are separate.

- **Allele inference**: System uses provided germlines; inferring novel alleles from repertoire data is separate tool (IgDiscover).

## Constitution Alignment

This specification aligns with the Sadie Germlines Module Constitution v1.0.0:

- **Principle I (Provider-Based)**: VDJbase implemented as provider following base interface ✓
- **Principle II (Priority-Based Merging)**: Priority ordering system tested and validated ✓
- **Principle III (Local-First)**: Offline operation verified as success criterion ✓
- **Principle IV (Staged Pipeline)**: Uses existing sources → normalized → igblast flow ✓
- **Principle V (Integration Compatibility)**: Backward compatibility and regression testing required ✓

All functional requirements map to constitutional principles. Testing requirements enforce compliance.
