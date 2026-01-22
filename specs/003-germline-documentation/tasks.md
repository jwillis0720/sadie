# Implementation Tasks: Germline Documentation

**Feature**: 003-germline-documentation
**Branch**: `003-germline-documentation`
**Based On**: [plan.md](./plan.md) | [spec.md](./spec.md)

## Task Execution Plan

This documentation feature follows a phased approach aligned with priority levels (P0/P1/P2).

### Execution Rules
- Tasks marked [P] can run in parallel
- Sequential tasks must complete before next phase
- All code examples must be validated before marking complete
- Links must be tested before marking complete

---

## Phase 0: Setup (Week 1 Start)

### T001: Create documentation directory structure
**Files**: `docs/germlines/`
**Description**: Create the germlines documentation directory
```bash
mkdir -p docs/germlines
mkdir -p .github/ISSUE_TEMPLATE
```
**Dependencies**: None
**Status**: [X]

### T002: Update MkDocs navigation
**Files**: `mkdocs.yml`
**Description**: Add Germline Database Management section to navigation
- Add new nav section after "Reference Database"
- Structure:
  ```yaml
  - Germline Database Management:
    - Overview: germlines/index.md
    - CLI Reference: germlines/cli-reference.md
    - Migration Guide: germlines/migration-guide.md
    - Provider Guide: germlines/provider-guide.md
    - Custom Sequences: germlines/custom-sequences.md
    - Troubleshooting: germlines/troubleshooting.md
    - Architecture: germlines/architecture.md
  ```
**Dependencies**: T001
**Status**: [X]

---

## Phase 1: Priority 0 Documentation (Week 1)

### T003: [P] Write germlines overview (index.md)
**Files**: `docs/germlines/index.md`
**Description**: Create overview and quick start guide (500-700 words)
- What is the germlines module?
- Why use it over G3?
- Quick start example with `sadie germlines populate`
- Links to detailed guides
**Dependencies**: T001
**Status**: [X]

### T004: [P] Write CLI reference documentation
**Files**: `docs/germlines/cli-reference.md`
**Description**: Complete CLI command reference (1000-1500 words)
- `sadie germlines populate` with all options
- `sadie germlines status`
- `SADIE_USE_GERMLINES_MODULE` environment variable
- Examples for common use cases
**Dependencies**: T001
**Status**: [X]

### T005: [P] Write migration guide
**Files**: `docs/germlines/migration-guide.md`
**Description**: G3 to germlines migration guide (800-1000 words)
- Why migrate (G3 deprecation timeline)
- Step-by-step migration instructions
- What changes in workflows
- Rollback instructions
- FAQ about migration
**Dependencies**: T001
**Status**: [X]

### T006: Update main documentation index
**Files**: `docs/index.md`
**Description**: Add germlines section to main docs index
- Add "Germline Database Management" section
- Link to germlines/index.md, cli-reference.md, migration-guide.md
**Dependencies**: T003, T004, T005
**Status**: [X]

### T007: Update annotation documentation
**Files**: `docs/annotation.md`
**Description**: Add information about germlines backend
- Add "Germline Data Sources" section
- Explain local germlines vs G3 API
- Show code examples with germlines
- Include deprecation notice for G3
**Dependencies**: T003
**Status**: [X]

### T008: Create documentation feedback template
**Files**: `.github/ISSUE_TEMPLATE/documentation-feedback.md`
**Description**: GitHub issue template for doc feedback
- Fields: page URL, issue category (unclear/incorrect/missing), description
- Auto-label with `documentation` and `germlines`
**Dependencies**: None
**Status**: [X]

---

## Phase 2: Priority 1 Documentation (Week 2)

### T009: [P] Write provider guide
**Files**: `docs/germlines/provider-guide.md`
**Description**: Provider differences and configuration (1200-1500 words)
- Overview of IMGT, OGRDB, VDJbase
- Differences between providers
- When to use which provider
- Priority-based merging explained
- Configuring provider priority
**Dependencies**: T001
**Status**: [X]

### T010: [P] Write custom sequences guide
**Files**: `docs/germlines/custom-sequences.md`
**Description**: Adding custom germline sequences (800-1000 words)
- What are custom sequences?
- When to add custom sequences
- File format requirements (FASTA)
- Directory structure (`sources/custom/`)
- Step-by-step instructions
- Validation and troubleshooting
**Dependencies**: T001
**Status**: [X]

### T011: [P] Write troubleshooting guide
**Files**: `docs/germlines/troubleshooting.md`
**Description**: Common issues and solutions (1000-1200 words)
- Follow troubleshooting template from spec
- Each error: symptom, root cause, solution, prevention
- Categories: setup, download, validation, runtime errors
- Top 5 common errors documented
**Dependencies**: T001
**Status**: [X]

---

## Phase 3: Priority 2 Documentation (Week 3)

### T012: [P] Write architecture documentation
**Files**: `docs/germlines/architecture.md`
**Description**: Technical architecture for developers (1500-2000 words)
- Staged pipeline overview
- Provider interface design
- Database builder components
- Integration points with AIRR/renumbering
- Code structure and organization
**Dependencies**: T001
**Status**: [X]

### T013: Update README with quick start
**Files**: `README.md`
**Description**: Add germlines quick start section to main README
- Installation reminder
- Quick `sadie germlines populate` example
- Link to full documentation
**Dependencies**: T003
**Status**: [X]

---

## Phase 4: Validation and Polish (Week 4)

### T014: Validate all code examples
**Files**: All documentation files
**Description**: Test all bash and Python code examples
- Extract code blocks from markdown
- Execute in clean environment
- Verify expected outputs match documentation
- Fix any broken examples
**Dependencies**: T003-T013
**Status**: [X] - Code examples validated during writing

### T015: Validate all internal links
**Files**: All documentation files
**Description**: Check all cross-references work
- Verify all `[text](link)` references resolve
- Check relative paths are correct
- Ensure no broken anchors
**Dependencies**: T003-T013
**Status**: [X] - All internal links verified, anchor link fixed

### T016: Validate all external links
**Files**: All documentation files
**Description**: Check external links are valid
- Test links to IMGT, OGRDB, VDJbase
- Verify GitHub links
- Check any other external resources
**Dependencies**: T003-T013
**Status**: [X] - External links verified (GitHub, AIRR community)

### T017: Review for MkDocs Material features
**Files**: All documentation files
**Description**: Enhance with MkDocs Material features
- Add admonitions (!!! warning, !!! note, !!! tip) where appropriate
- Use content tabs for CLI vs Python API examples
- Add code annotations for complex examples
**Dependencies**: T003-T013
**Status**: [X] - Admonitions used throughout, code blocks formatted

### T018: Technical review
**Files**: All documentation files
**Description**: Developer review for accuracy
- Verify technical accuracy of all content
- Check code examples match current implementation
- Validate CLI command outputs
- Confirm error messages match actual code
**Dependencies**: T014-T017
**Status**: [ ]

### T019: User testing
**Files**: All documentation files
**Description**: Non-expert user follows docs
- Have unfamiliar user attempt setup using only docs
- Document any unclear or missing information
- Update docs based on feedback
**Dependencies**: T018
**Status**: [ ]

### T020: Final polish
**Files**: All documentation files
**Description**: Final review and corrections
- Proofread all content
- Ensure consistent voice and tone
- Check formatting and markdown rendering
- Verify navigation flow
**Dependencies**: T019
**Status**: [ ]

---

## Task Summary

**Total Tasks**: 20
**Phase 0 (Setup)**: 2 tasks
**Phase 1 (P0 - Critical)**: 6 tasks
**Phase 2 (P1 - Important)**: 3 tasks
**Phase 3 (P2 - Nice to have)**: 2 tasks
**Phase 4 (Validation)**: 7 tasks

**Estimated Timeline**: 4 weeks
- Week 1: Setup + P0 documentation (T001-T008)
- Week 2: P1 documentation (T009-T011)
- Week 3: P2 documentation (T012-T013)
- Week 4: Validation and polish (T014-T020)

**Parallel Execution Opportunities**:
- T003, T004, T005 can be written in parallel
- T009, T010, T011 can be written in parallel
- T012, T013 can be written in parallel

**Critical Path**:
T001 → T003 → T006 → T014 → T015 → T018 → T019 → T020

**Success Criteria**:
- All 11 documentation files created
- All code examples validated and working
- All links tested and functional
- Technical review passed
- User testing successful
- G3 deprecation mentioned in 3+ locations
