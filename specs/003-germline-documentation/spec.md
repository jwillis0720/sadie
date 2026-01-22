# Feature Specification: Germline Database Documentation

**Feature Branch**: `sadie_docs`
**Created**: 2026-01-22
**Status**: Planning
**Based On**: [002-germline-integration](../002-germline-integration/spec.md)

## Overview

Create comprehensive user-facing documentation for the germline database integration feature implemented in [002-germline-integration](../002-germline-integration/). This documentation will enable users to:
- Understand and use the new `sadie germlines` CLI commands
- Migrate from G3 API to the local germlines module
- Add custom germline sequences
- Troubleshoot common issues
- Understand the architecture for advanced use cases

## User Scenarios *(mandatory)*

### User Story 1 - New User Setup (Priority: P0)

As a new SADIE user, I want clear instructions on how to set up the germline databases, so that I can start using SADIE for AIRR annotation without confusion.

**Why this priority**: First-time setup is the critical path for all new users. Poor documentation here blocks adoption.

**Acceptance Scenarios**:

1. **Given** I've installed SADIE, **When** I read the germlines documentation, **Then** I understand I need to run `sadie germlines populate` before my first analysis.

2. **Given** I run `sadie germlines populate`, **When** I check the documentation for what's happening, **Then** I understand the download progress, build steps, and expected completion time.

3. **Given** the download completes, **When** I run `sadie germlines status`, **Then** I can verify my setup is correct by comparing with documented expected output.

---

### User Story 2 - G3 Migration (Priority: P0)

As an existing SADIE user using G3, I want a clear migration guide, so that I can transition to the germlines module before G3 is deprecated.

**Why this priority**: Existing users need to migrate before 2026-06-01 G3 deprecation deadline. This is time-sensitive.

**Acceptance Scenarios**:

1. **Given** I'm using G3 API, **When** I see the deprecation warning, **Then** I can find clear migration documentation explaining why and how to migrate.

2. **Given** I follow the migration guide, **When** I complete the steps, **Then** my existing workflows continue to work unchanged with the germlines module.

3. **Given** I encounter issues during migration, **When** I check the troubleshooting guide, **Then** I find solutions for common migration problems.

---

### User Story 3 - CLI Command Reference (Priority: P0)

As a user running germline commands, I want complete CLI reference documentation, so that I understand all available options and their effects.

**Why this priority**: Without CLI reference, users cannot use the new commands effectively or understand available options.

**Acceptance Scenarios**:

1. **Given** I need to download specific species, **When** I check the CLI reference, **Then** I find examples showing how to use `--species` option.

2. **Given** I want to understand what `--dry-run` does, **When** I read the CLI reference, **Then** I see a clear explanation with example output.

3. **Given** I need to force re-download, **When** I check the documentation, **Then** I find the `--force` flag documented with use cases.

---

### User Story 4 - Custom Sequences (Priority: P1)

As a researcher with novel germline sequences, I want documentation on adding custom sequences, so that I can include my sequences in the analysis pipeline.

**Why this priority**: Custom sequences are a key differentiator from G3, but require clear documentation to use correctly.

**Acceptance Scenarios**:

1. **Given** I have novel sequences, **When** I read the custom sequences guide, **Then** I understand the file format, naming conventions, and directory structure required.

2. **Given** I add custom sequences, **When** I consult the documentation, **Then** I understand how priority-based merging works and that my custom sequences override IMGT.

3. **Given** my custom sequences fail validation, **When** I check the troubleshooting guide, **Then** I find common formatting issues and solutions.

---

### User Story 5 - Troubleshooting (Priority: P1)

As a user encountering errors, I want a troubleshooting guide, so that I can resolve issues without external help.

**Why this priority**: Self-service troubleshooting reduces support burden and improves user experience.

**Acceptance Scenarios**:

1. **Given** I see "Species not found" error, **When** I search the troubleshooting guide, **Then** I find this specific error with step-by-step resolution.

2. **Given** downloads are failing, **When** I check troubleshooting, **Then** I find common network-related issues and solutions.

3. **Given** I want to debug issues, **When** I read the troubleshooting guide, **Then** I learn how to enable verbose logging for detailed diagnostics.

---

### Edge Cases

- User has no internet connection during initial setup: Documentation should clarify that initial `populate` requires internet
- User's disk is full: Troubleshooting should mention disk space requirements (500MB-3GB)
- User tries to use germlines before populating: Error message should point to documentation
- Multiple SADIE versions installed: Documentation should clarify version compatibility
- User wants to use specific provider only: CLI reference should show provider filtering examples

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Documentation MUST include complete reference for `sadie germlines populate` command with all options, examples, and expected output.

- **FR-002**: Documentation MUST include complete reference for `sadie germlines status` command with output interpretation.

- **FR-003**: Documentation MUST include step-by-step migration guide from G3 API to germlines module with verification steps.

- **FR-004**: Documentation MUST explain the `SADIE_USE_GERMLINES_MODULE` environment variable and its usage.

- **FR-005**: Documentation MUST include guide for adding custom germline sequences with file format specifications and validation rules.

- **FR-006**: Documentation MUST include troubleshooting guide covering common errors: "Species not found", "Download failed", "Validation error", "Disk space". Each error entry MUST include: exact error message text, root cause explanation, step-by-step fix instructions, and prevention tips.

- **FR-007**: Documentation MUST explain provider differences (IMGT vs OGRDB vs VDJbase) and priority-based merging.

- **FR-008**: Documentation MUST include examples showing both command-line and Python API usage where applicable.

- **FR-009**: Documentation MUST be structured with clear navigation: overview → quick start → detailed guides → troubleshooting.

- **FR-010**: Documentation MUST include code examples that are tested and working.

- **FR-011**: Documentation MUST mention G3 deprecation timeline (removal after 2026-06-01) in appropriate locations.

- **FR-012**: Documentation MUST be integrated into existing docs structure with cross-links to related pages.

- **FR-013**: Documentation MUST include a GitHub issue template for user feedback (`.github/ISSUE_TEMPLATE/documentation-feedback.md`) with fields for page URL, issue category, and description.

### Non-Functional Requirements

- **NFR-001**: Documentation MUST be written in clear, concise language accessible to users with basic command-line experience.

- **NFR-002**: All code examples MUST be complete, runnable, and include expected output where helpful.

- **NFR-003**: Documentation MUST follow existing SADIE documentation style and formatting conventions.

- **NFR-004**: All internal documentation links MUST be valid and tested.

- **NFR-005**: Documentation pages MUST load quickly and be readable on both desktop and mobile devices, leveraging MkDocs Material's responsive design.

- **NFR-006**: Documentation MUST use MkDocs Material features (admonitions, code tabs, content tabs) where appropriate to enhance readability.

- **NFR-007**: Documentation updates MUST be included in PRs that change CLI commands, options, error messages, or user-facing behavior. A designated maintainer MUST review documentation monthly for accuracy and completeness.

- **NFR-008**: All code examples in documentation MUST be automatically validated via pytest-based doctest runner in CI/CD pipeline. Examples that fail validation MUST block PR merges.

### Key Entities

- **CLI Command Documentation**: Reference material for `sadie germlines` commands
- **Migration Guide**: Step-by-step instructions for G3 → germlines transition
- **Provider Guide**: Explanation of IMGT/OGRDB/VDJbase differences
- **Custom Sequences Guide**: Instructions for adding user sequences
- **Troubleshooting Guide**: Common issues and solutions
- **Architecture Documentation**: Technical overview for developers

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: New users can complete initial setup (populate + verify) following documentation without external help.

- **SC-002**: Existing G3 users can migrate to germlines module following the migration guide with no workflow changes.

- **SC-003**: All CLI commands (`sadie germlines populate` and `sadie germlines status`) are documented with at least 3 usage examples each.

- **SC-004**: All code examples in documentation execute without errors when tested.

- **SC-005**: Troubleshooting guide addresses the top 5 most common error messages with solutions.

- **SC-006**: Documentation passes technical review by developer who implemented the feature.

- **SC-007**: Documentation passes user testing by someone unfamiliar with the feature who can successfully complete setup.

- **SC-008**: All internal documentation links are valid (no 404 errors).

- **SC-009**: G3 deprecation notice is mentioned in at least 3 relevant locations (migration guide, CLI reference, main README).

- **SC-010**: Documentation feedback mechanism (GitHub issue template) is established and documented, with quarterly analytics reviews showing user engagement metrics (page views, time on page, bounce rate).

## Documentation Deliverables

### Priority 0 (Critical - Week 1)

1. **docs/germlines/index.md** (500-700 words)
   - Overview, benefits, quick start

2. **docs/germlines/cli-reference.md** (1000-1500 words)
   - Complete command reference with examples

3. **docs/germlines/migration-guide.md** (800-1000 words)
   - Step-by-step G3 migration instructions

4. **Update docs/index.md** (add germlines section)

5. **Update docs/annotation.md** (add germlines backend info)

6. **Create `.github/ISSUE_TEMPLATE/documentation-feedback.md`**
   - GitHub issue template for doc feedback

### Priority 1 (Important - Week 2)

7. **docs/germlines/provider-guide.md** (1200-1500 words)
   - Provider differences, priority merging

8. **docs/germlines/custom-sequences.md** (800-1000 words)
   - Adding custom sequences step-by-step

9. **docs/germlines/troubleshooting.md** (1000-1200 words)
   - Common issues and solutions
   - Each error entry includes: exact error message, root cause, step-by-step fix, prevention tips
   - Organized by category: setup errors, download errors, validation errors, runtime errors

### Priority 2 (Nice to have - Week 3)

10. **docs/germlines/architecture.md** (1500-2000 words)
   - Technical architecture for developers

11. **Update README.md** (add germlines quick start)

## Style and Formatting Guidelines

### Voice and Tone
- Clear, concise, professional
- Action-oriented (imperatives: "Download...", "Check...")
- Friendly but not casual

### Code Examples
- Always complete and runnable
- Include both CLI and Python API where applicable
- Show expected output for verification

### Structure
- Brief overview at page top
- Prerequisites section when needed
- Step-by-step instructions
- Examples and use cases
- Next steps or related links

### Formatting
- Code blocks with syntax highlighting
- Command prompts (`$` for bash)
- Tables for comparing options
- MkDocs Material admonitions for warnings/notes (e.g., `!!! warning`, `!!! note`, `!!! tip`)
- Content tabs for alternative approaches (e.g., CLI vs Python API)
- Code annotations for complex examples

## Testing and Validation

### Documentation Testing

1. **Technical Accuracy**
   - Developer reviews for correctness
   - All code examples are tested via automated pytest-based doctest runner
   - CI/CD pipeline extracts code blocks from markdown and validates execution
   - Both bash commands and Python code examples are validated

2. **User Testing**
   - Non-expert follows instructions
   - Identifies unclear or missing information

3. **User Feedback Collection**
   - GitHub issue template for documentation feedback (`.github/ISSUE_TEMPLATE/documentation-feedback.md`)
   - Template includes: page URL, issue category (unclear, incorrect, missing), description
   - Issues tagged with `documentation` and `germlines` labels for tracking

4. **Link Validation**
   - All internal links work
   - External links are valid

5. **Completeness Check**
   - All CLI options documented
   - All common errors covered
   - Migration path is clear

### Automated Code Example Validation

**Setup:**
- pytest plugin (pytest-markdown-docs or similar) configured to extract code blocks
- CI/CD job runs validation on every PR and merge to main
- Code blocks tagged with language identifier (bash, python) are validated
- Examples marked with special comments can be excluded from validation if needed

**Validation Process:**
1. Extract all bash and python code blocks from documentation markdown files
2. Set up isolated test environment with SADIE installed
3. Execute bash commands in sequence (capture stdout/stderr)
4. Execute Python code blocks (validate expected output if specified)
5. Report failures with file, line number, and error message

**Failure Handling:**
- PR checks fail if any code example fails validation
- Monthly maintenance includes fixing broken examples
- Breaking changes require doc updates in same PR

## Maintenance Strategy

### Continuous Updates
- Developers MUST update documentation in the same PR when changing:
  - CLI command options or behavior
  - Error messages or codes
  - Environment variables
  - User-facing APIs

### Monthly Reviews
- Designated documentation maintainer reviews all germlines documentation monthly
- Checklist includes:
  - Code examples still execute correctly
  - CLI command output matches current version
  - Links are valid (internal and external)
  - Version-specific content is accurate
  - Deprecation notices are current

### Review Process
- Maintainer creates GitHub issue for any found discrepancies
- Issues tagged with `documentation` and `germlines` labels
- Priority based on impact: P0 (blocks users), P1 (misleading), P2 (cosmetic)

### Quarterly Analytics Review
- Review documentation analytics (page views, time on page, bounce rate, search queries)
- Identify high-traffic pages that need enhancement
- Identify low-traffic pages that may need better discovery/links
- Review GitHub issues tagged `documentation` + `germlines` for common themes
- Create improvement backlog based on quantitative and qualitative feedback
- Report findings to team with recommended improvements

## Timeline

- **Week 1**: P0 documentation (critical setup and migration)
- **Week 2**: P1 documentation (advanced guides)
- **Week 3**: P2 documentation (architecture)
- **Week 4**: Review, testing, polish

**Target Completion**: Before 2026-02-28 (well ahead of 2026-06-01 G3 deprecation)

## Dependencies

- Requires completed [002-germline-integration](../002-germline-integration/) implementation
- Access to test environment for validating examples
- Review from original feature developer
- MkDocs Material documentation platform (already in project dependencies: mkdocs-material in pyproject.toml)
- pytest with markdown doctest plugin (e.g., pytest-markdown-docs or pytest-codeblocks) for automated code example validation

## Clarifications

### Session 2026-01-22

- Q: What documentation hosting platform will be used for these docs? → A: MkDocs Material
- Q: Who is responsible for keeping documentation up-to-date when code changes? → A: PR requirement + dedicated maintainer reviews docs monthly for accuracy
- Q: What method should be used to validate code examples are working? → A: Automated pytest-based doctest runner extracts and validates examples
- Q: How should documentation measure user success and gather feedback? → A: GitHub issue template + quarterly analytics review (page views, time on page)
- Q: What level of detail should troubleshooting guide include for errors? → A: Error message text + root cause + step-by-step fix + prevention tips

## Open Questions

- Should we create video tutorials for complex workflows? (Deferred to future work)
- Should we include performance benchmarks in documentation? (Yes, in architecture.md)
- Do we need translations for non-English speakers? (No, English only for now)

## Related Specifications

- [002-germline-integration/spec.md](../002-germline-integration/spec.md) - Feature specification
- [002-germline-integration/plan.md](../002-germline-integration/plan.md) - Implementation plan
- [src/sadie/germlines/INTEGRATION_GUIDE.md](../../src/sadie/germlines/INTEGRATION_GUIDE.md) - Developer integration guide
- [src/sadie/germlines/README.md](../../src/sadie/germlines/README.md) - Module README

## Troubleshooting Entry Template

Each troubleshooting entry should follow this structure:

```markdown
### Error: [Exact error message text]

**Symptom**:
[Describe what the user sees, including full error output]

**Root Cause**:
[Explain why this error occurs - the underlying issue]

**Solution**:
1. [Step-by-step fix instructions]
2. [Each step should be actionable and specific]
3. [Include verification step to confirm fix worked]

**Prevention**:
[Tips to avoid this error in the future]

**Related Issues**:
- [Links to related errors or edge cases]
```

**Example**:
```markdown
### Error: "Species 'human' not found in germlines module"

**Symptom**:
```
ValueError: Species 'human' not found in germlines module at /path/to/internal_data/human.
Build germlines databases with: update_databases('human').
```

**Root Cause**:
Germlines data has not been downloaded and built for the requested species.

**Solution**:
1. Run `sadie germlines populate -p imgt -s human`
2. Wait for download to complete (2-5 minutes)
3. Verify with `sadie germlines status` - should show "Up-to-date"
4. Retry your original command

**Prevention**:
Run `sadie germlines populate` after installing SADIE for the first time. This downloads data for all species.

**Related Issues**:
- Download failures → see "Download failed for species"
- Network issues → see "Connection timeout during populate"
```

## Notes

- Documentation is user-facing and should be accessible to researchers without deep technical background
- Examples should use realistic scenarios (e.g., "human", "mouse") not abstract placeholders
- Error messages in troubleshooting should match actual error text from code
- Consider adding diagrams for architecture and data flow (can be added in revision)
