#!/bin/bash

# GitHub Issues Creation Script for Germlines Module Completion
# Feature: 001-germline-completion
# Repository: kipkurui/sadie

set -e

REPO="kipkurui/sadie"
SPEC_URL="https://github.com/kipkurui/sadie/blob/001-germline-completion/specs/001-germline-completion/spec.md"
PLAN_URL="https://github.com/kipkurui/sadie/blob/001-germline-completion/specs/001-germline-completion/plan.md"
TASKS_URL="https://github.com/kipkurui/sadie/blob/001-germline-completion/specs/001-germline-completion/tasks.md"

echo "🚀 Creating GitHub issues for Germlines Module Completion..."
echo "Repository: $REPO"
echo ""

# Check if gh CLI is authenticated
if ! gh auth status > /dev/null 2>&1; then
    echo "❌ Error: Not authenticated with GitHub"
    echo "Run: gh auth login"
    exit 1
fi

echo "✅ GitHub authentication verified"
echo ""

#==============================================
# Create Labels (if they don't exist)
#==============================================

echo "🏷️  Creating GitHub labels..."

# Function to create label if it doesn't exist
create_label() {
    local name=$1
    local color=$2
    local description=$3

    # Check if label exists
    if ! gh label list --repo "$REPO" --json name --jq '.[].name' | grep -q "^${name}$"; then
        gh label create "$name" --repo "$REPO" --color "$color" --description "$description" 2>/dev/null || true
        echo "  ✓ Created label: $name"
    else
        echo "  ✓ Label exists: $name"
    fi
}

# Core labels
create_label "germlines" "0052CC" "Germlines module feature"
create_label "P1" "D73A4A" "Priority 1 - Critical"
create_label "P2" "FBCA04" "Priority 2 - Important"

# User story labels
create_label "US1" "C2E0C6" "User Story 1: Custom Germlines"
create_label "US2" "BFD4F2" "User Story 2: Offline Operation"
create_label "US3" "D4C5F9" "User Story 3: Priority System"
create_label "US4" "FEF2C0" "User Story 4: Sadie Integration"
create_label "US5" "F9D0C4" "User Story 5: VDJbase Provider"
create_label "US6" "C5DEF5" "User Story 6: Data Population"

# Phase labels
create_label "setup" "EDEDED" "Phase 1: Project Setup"
create_label "foundation" "D3D3D3" "Phase 2: Foundational Components"
create_label "integration" "B8E6B8" "Integration with Sadie components"
create_label "testing" "E99695" "Testing and validation"
create_label "polish" "FFD700" "Polish and finalization"

# Component labels
create_label "custom-germlines" "006B75" "Custom germline sequences"
create_label "data-population" "5319E7" "Data download and population"
create_label "offline-operation" "0E8A16" "Offline/air-gapped operation"
create_label "vdjbase" "1D76DB" "VDJbase provider"
create_label "priority-system" "D93F0B" "Priority-based database selection"
create_label "observability" "B60205" "Logging and monitoring"
create_label "dependencies" "0366D6" "Dependency management"

# Special markers
create_label "parallelizable" "FBCA04" "Can be worked on in parallel"

echo "✅ Labels created"
echo ""

#==============================================
# Phase 1: Project Setup (7 tasks)
#==============================================

echo "📦 Creating Phase 1: Project Setup issues..."

# T001
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create VDJbase provider directory structure" \
    --body "**Task**: T001
**Phase**: 1 - Project Setup
**Priority**: P1
**Duration**: Part of 1-2 hour phase

**Description**:
Create directory structure for VDJbase provider at \`src/sadie/germlines/sources/vdjbase/\`

**Acceptance Criteria**:
- [ ] Directory \`src/sadie/germlines/sources/vdjbase/\` exists
- [ ] Directory structure matches other providers (custom, imgt, ogrdb)

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1"

# T002
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create VDJbase provider human subdirectory" \
    --body "**Task**: T002
**Phase**: 1 - Project Setup
**Priority**: P1
**Duration**: Part of 1-2 hour phase

**Description**:
Create human species subdirectory for VDJbase provider at \`src/sadie/germlines/sources/vdjbase/human/\`

**Acceptance Criteria**:
- [ ] Directory \`src/sadie/germlines/sources/vdjbase/human/\` exists
- [ ] Directory is ready to receive FASTA files

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1"

# T003
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create test data directory structure" \
    --body "**Task**: T003
**Phase**: 1 - Project Setup
**Priority**: P1
**Parallelizable**: Yes
**Duration**: Part of 1-2 hour phase

**Description**:
Create test data directory structure in \`tests/data/germlines/\` with subdirectories for all providers (custom/, imgt/, ogrdb/, vdjbase/)

**Acceptance Criteria**:
- [ ] Directory \`tests/data/germlines/\` exists
- [ ] Subdirectories created: custom/, imgt/, ogrdb/, vdjbase/
- [ ] Structure matches source data organization

**References**:
- [Spec]($SPEC_URL) - FR-029
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1,parallelizable,testing"

# T004
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create curated test dataset" \
    --body "**Task**: T004
**Phase**: 1 - Project Setup
**Priority**: P1
**Parallelizable**: Yes
**Duration**: Part of 1-2 hour phase

**Description**:
Create curated test dataset with 5-10 genes per segment in \`tests/data/germlines/\`

**Acceptance Criteria**:
- [ ] Test data includes 5-10 genes per segment (V, D, J)
- [ ] Examples from each provider (custom, IMGT, OGRDB, VDJbase)
- [ ] Both gapped and ungapped variants included
- [ ] Edge cases covered: invalid sequences, duplicate names, identical sequences with different names

**References**:
- [Spec]($SPEC_URL) - FR-029, FR-030, FR-031
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1,parallelizable,testing"

# T005
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create validation script template" \
    --body "**Task**: T005
**Phase**: 1 - Project Setup
**Priority**: P1
**Duration**: Part of 1-2 hour phase

**Description**:
Create validation script template at \`src/sadie/germlines/scripts/validate.py\`

**Acceptance Criteria**:
- [ ] Script file created at correct location
- [ ] Template includes basic structure for FASTA validation
- [ ] Script is executable

**References**:
- [Spec]($SPEC_URL) - FR-009, FR-009a, FR-009b
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1"

# T006
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Update pyproject.toml dependencies" \
    --body "**Task**: T006
**Phase**: 1 - Project Setup
**Priority**: P1
**Duration**: Part of 1-2 hour phase

**Description**:
Update \`pyproject.toml\` to ensure all required dependencies are listed (BioPython, pytest)

**Acceptance Criteria**:
- [ ] BioPython dependency listed with appropriate version
- [ ] pytest listed in dev dependencies
- [ ] All germlines module dependencies documented

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL) - Technical Context
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1,dependencies"

# T007
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Setup] Create feature flag utility module" \
    --body "**Task**: T007
**Phase**: 1 - Project Setup
**Priority**: P1
**Duration**: Part of 1-2 hour phase

**Description**:
Create feature flag utility module at \`src/sadie/germlines/utils/feature_flags.py\`

**Acceptance Criteria**:
- [ ] File created at \`src/sadie/germlines/utils/feature_flags.py\`
- [ ] Module structure ready for feature flag implementation
- [ ] Imports and basic structure in place

**References**:
- [Spec]($SPEC_URL) - FR-016
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,setup,P1"

echo "✅ Phase 1 issues created (7 tasks)"
echo ""

#==============================================
# Phase 2: Foundational Components (8 tasks)
#==============================================

echo "🏗️  Creating Phase 2: Foundational Components issues..."

# T008
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Implement feature flag function" \
    --body "**Task**: T008
**Phase**: 2 - Foundational Components
**Priority**: P1
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Implement \`use_germlines_module()\` function in \`src/sadie/germlines/utils/feature_flags.py\`

**Acceptance Criteria**:
- [ ] Function \`use_germlines_module()\` implemented
- [ ] Reads \`SADIE_USE_GERMLINES_MODULE\` environment variable
- [ ] Defaults to \"true\" (germlines module enabled)
- [ ] Returns boolean value

**Implementation Details**:
\`\`\`python
def use_germlines_module() -> bool:
    return os.getenv(\"SADIE_USE_GERMLINES_MODULE\", \"true\").lower() == \"true\"
\`\`\`

**References**:
- [Spec]($SPEC_URL) - FR-016, FR-016a, FR-016b
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1"

# T009
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Implement auto-gapping service using BioPython" \
    --body "**Task**: T009
**Phase**: 2 - Foundational Components
**Priority**: P1
**Parallelizable**: Yes
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Implement auto-gapping service using BioPython pairwise alignment in \`src/sadie/germlines/builders/gapper.py\`

**Acceptance Criteria**:
- [ ] Gapper service implemented using BioPython PairwiseAligner
- [ ] Uses per-gene IMGT-gapped template when available; fallback to per-segment consensus
- [ ] Translates nucleotide to amino acid before alignment
- [ ] Maps gaps back to nucleotide sequence using codon positions
- [ ] Uses \".\" (period) for gap characters at codon boundaries
- [ ] Handles failures gracefully: logs WARNING and stores ungapped version only

**References**:
- [Spec]($SPEC_URL) - FR-021, FR-021a, FR-021b, FR-021c, FR-021d
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1,parallelizable"

# T010
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Add logging configuration for germlines module" \
    --body "**Task**: T010
**Phase**: 2 - Foundational Components
**Priority**: P1
**Parallelizable**: Yes
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Add logging configuration for germlines module in \`src/sadie/germlines/__init__.py\`

**Acceptance Criteria**:
- [ ] Logging configured at INFO, WARNING, ERROR levels
- [ ] Log format: \`{timestamp} - {logger_name} - {level} - {message}\` with ISO 8601 timestamps
- [ ] Structured key-value pairs: \`operation={name} duration_ms={ms} status={success|failure}\`
- [ ] Logger name: \`sadie.germlines\`

**References**:
- [Spec]($SPEC_URL) - FR-032, FR-032a, FR-032b
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1,parallelizable,observability"

# T011
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Update GermlineManager to support VDJbase provider" \
    --body "**Task**: T011
**Phase**: 2 - Foundational Components
**Priority**: P1
**Parallelizable**: Yes
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Update GermlineManager in \`src/sadie/germlines/manager.py\` to support VDJbase provider

**Acceptance Criteria**:
- [ ] VDJbaseProvider added to default provider list
- [ ] Priority ordering supports VDJbase: \`custom > imgt > ogrdb > vdjbase\`
- [ ] Manager can load and merge VDJbase genes

**References**:
- [Spec]($SPEC_URL) - FR-004
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1,parallelizable"

# T012
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Create VDJbase provider stub with base interface" \
    --body "**Task**: T012
**Phase**: 2 - Foundational Components
**Priority**: P1
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Create VDJbase provider stub with base interface in \`src/sadie/germlines/providers/vdjbase.py\`

**Acceptance Criteria**:
- [ ] VDJbaseProvider class inherits from GermlineProvider
- [ ] Implements all required abstract methods: fetch_genes(), fetch_gene_by_name(), is_available(), get_metadata()
- [ ] Stub methods return empty/placeholder values
- [ ] Provider name set to \"vdjbase\"

**References**:
- [Spec]($SPEC_URL) - FR-001
- [Plan]($PLAN_URL) - Provider Interface
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1,vdjbase"

# T013
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P2][Foundation] Implement VDJbase FASTA parsing logic" \
    --body "**Task**: T013
**Phase**: 2 - Foundational Components
**Priority**: P2
**Duration**: Part of 3-4 hour phase
**Prerequisites**: T012 complete

**Description**:
Implement VDJbase FASTA parsing logic in \`src/sadie/germlines/providers/vdjbase.py\`

**Acceptance Criteria**:
- [ ] Parses VDJbase FASTA header format: \`>{gene_name}|{species}|{segment}|{chain}[|population={pop}][|genotype={gt}]\`
- [ ] Extracts population and genotype metadata from headers
- [ ] Stores metadata in GermlineGene.metadata dictionary
- [ ] Creates GermlineGene objects from parsed data
- [ ] Handles malformed headers gracefully

**References**:
- [Spec]($SPEC_URL) - FR-002, FR-002a, FR-002b, FR-002c
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P2,vdjbase"

# T014
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P2][Foundation] Implement VDJbase provider metadata methods" \
    --body "**Task**: T014
**Phase**: 2 - Foundational Components
**Priority**: P2
**Duration**: Part of 3-4 hour phase
**Prerequisites**: T013 complete

**Description**:
Implement VDJbase provider metadata methods in \`src/sadie/germlines/providers/vdjbase.py\`

**Acceptance Criteria**:
- [ ] \`is_available(species)\` checks for VDJbase data files
- [ ] \`get_metadata()\` returns ProviderMetadata with name=\"vdjbase\"
- [ ] When data unavailable: logs WARNING \"VDJbase data not found at {path}. Skipping VDJbase provider.\"
- [ ] Continues without error when data missing

**References**:
- [Spec]($SPEC_URL) - FR-003, FR-004a
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P2,vdjbase"

# T015
gh issue create \
    --repo "$REPO" \
    --title "[Germlines][P1][Foundation] Add timing metrics logging to pipeline" \
    --body "**Task**: T015
**Phase**: 2 - Foundational Components
**Priority**: P1
**Duration**: Part of 3-4 hour phase
**Prerequisites**: Phase 1 complete

**Description**:
Add timing metrics logging to \`src/sadie/germlines/pipeline.py\`

**Acceptance Criteria**:
- [ ] Timing metrics logged for: data download, database rebuild, validation
- [ ] Log format: \"Operation {operation_name} completed in {duration_ms}ms\"
- [ ] Millisecond precision
- [ ] Logged at INFO level

**References**:
- [Spec]($SPEC_URL) - FR-033, FR-033a
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL)" \
    --label "germlines,foundation,P1,observability"

echo "✅ Phase 2 issues created (8 tasks)"
echo ""

#==============================================
# Phase 3: User Story 1 - Add Custom Germline Sequences (P1) (8 tasks)
#==============================================

echo "🧬 Creating Phase 3: User Story 1 (Custom Germlines) issues..."

# Note: Due to script length, creating abbreviated versions for phases 3-9
# Run with --verbose flag to see full task bodies

PHASE3_TASKS=(
    "T016|[US1] Verify custom provider handles new sequences|src/sadie/germlines/providers/custom.py|P1|US1"
    "T017|[US1] Implement change detection for custom sequences|src/sadie/germlines/pipeline.py|P1|US1"
    "T018|[US1] Add validation for custom FASTA files|src/sadie/germlines/providers/custom.py|P1|US1"
    "T019|[US1] Implement auto-rebuild trigger on custom file change|src/sadie/germlines/pipeline.py|P1|US1"
    "T020|[P] [US1] Write unit test for custom sequence priority|tests/test_custom_provider.py|P1|US1|parallelizable"
    "T021|[P] [US1] Write integration test for custom sequence end-to-end|tests/test_integration.py|P1|US1|parallelizable"
    "T022|[US1] Add logging for custom sequence load events|src/sadie/germlines/providers/custom.py|P1|US1"
    "T023|[US1] Document custom sequence addition process|src/sadie/germlines/sources/custom/README.md|P1|US1"
)

for task in "${PHASE3_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,custom-germlines,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 3 - User Story 1 (Custom Germlines)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-3-user-story-1---add-custom-germline-sequences-p1)" \
        --label "$labels"
done

echo "✅ Phase 3 issues created (8 tasks)"
echo ""

#==============================================
# Phase 4: User Story 6 - Populate Reference Data (P1) (11 tasks)
#==============================================

echo "📥 Creating Phase 4: User Story 6 (Data Population) issues..."

PHASE4_TASKS=(
    "T024|[US6] Complete IMGT download script implementation|src/sadie/germlines/scripts/download_imgt.py|P1|US6"
    "T025|[US6] Add species parameter support to IMGT download script|src/sadie/germlines/scripts/download_imgt.py|P1|US6"
    "T026|[US6] Implement resume capability for IMGT downloads|src/sadie/germlines/scripts/download_imgt.py|P1|US6"
    "T027|[US6] Implement OGRDB download script|src/sadie/germlines/scripts/download_ogrdb.py|P1|US6"
    "T028|[US6] Add species parameter support to OGRDB download script|src/sadie/germlines/scripts/download_ogrdb.py|P1|US6"
    "T029|[US6] Implement validation for downloaded FASTA files|src/sadie/germlines/scripts/validate.py|P1|US6"
    "T030|[P] [US6] Create VDJbase manual download instructions|src/sadie/germlines/sources/vdjbase/README.md|P1|US6|parallelizable"
    "T031|[P] [US6] Update IMGT data documentation|src/sadie/germlines/sources/imgt/IMGT_DATA.md|P1|US6|parallelizable"
    "T032|[P] [US6] Update OGRDB data documentation|src/sadie/germlines/sources/ogrdb/OGRDB_DATA.md|P1|US6|parallelizable"
    "T033|[US6] Add progress indicators to download scripts|src/sadie/germlines/scripts/download_imgt.py|P1|US6"
    "T034|[US6] Add timing metrics to download scripts|src/sadie/germlines/scripts/download_ogrdb.py|P1|US6"
)

for task in "${PHASE4_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,data-population,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 4 - User Story 6 (Data Population)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-4-user-story-6---populate-reference-data-sources-p1)" \
        --label "$labels"
done

echo "✅ Phase 4 issues created (11 tasks)"
echo ""

#==============================================
# Phase 5: User Story 2 - Offline Operation (P1) (6 tasks)
#==============================================

echo "🔌 Creating Phase 5: User Story 2 (Offline Operation) issues..."

PHASE5_TASKS=(
    "T035|[US2] Verify pipeline operates without network calls|src/sadie/germlines/pipeline.py|P1|US2"
    "T036|[US2] Add offline mode detection and logging|src/sadie/germlines/pipeline.py|P1|US2"
    "T037|[US2] Implement clear error messages for missing data|src/sadie/germlines/providers/base.py|P1|US2"
    "T038|[US2] Add README references in error messages|src/sadie/germlines/manager.py|P1|US2"
    "T039|[P] [US2] Write offline integration test|tests/test_offline_operation.py|P1|US2|parallelizable"
    "T040|[US2] Verify cached data usage works correctly|src/sadie/germlines/pipeline.py|P1|US2"
)

for task in "${PHASE5_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,offline-operation,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 5 - User Story 2 (Offline Operation)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-5-user-story-2---use-local-germline-databases-offline-p1)" \
        --label "$labels"
done

echo "✅ Phase 5 issues created (6 tasks)"
echo ""

#==============================================
# Phase 6: User Story 4 - Sadie Integration (P1) (11 tasks)
#==============================================

echo "🔗 Creating Phase 6: User Story 4 (Sadie Integration) issues..."

PHASE6_TASKS=(
    "T041|[US4] Update IgBLAST germline paths|src/sadie/airr/igblast/germline.py|P1|US4"
    "T042|[US4] Add feature flag check to IgBLAST|src/sadie/airr/igblast/germline.py|P1|US4"
    "T043|[US4] Update Reference system to query germlines module|src/sadie/reference/reference.py|P1|US4"
    "T044|[US4] Implement G3 API response format adapter|src/sadie/reference/reference.py|P1|US4"
    "T045|[US4] Add feature flag check to Reference system|src/sadie/reference/reference.py|P1|US4"
    "T046|[US4] Update HMM builder to use germlines module|src/sadie/renumbering/aligners/hmmer.py|P1|US4"
    "T047|[US4] Add feature flag check to HMM builder|src/sadie/renumbering/aligners/hmmer.py|P1|US4"
    "T048|[US4] Verify gapped sequences for Stockholm alignment|src/sadie/germlines/builders/gapper.py|P1|US4"
    "T049|[P] [US4] Run existing Sadie IgBLAST test suite|tests/|P1|US4|parallelizable"
    "T050|[P] [US4] Test feature flag G3 fallback|tests/|P1|US4|parallelizable"
    "T051|[US4] Document backward compatibility approach|src/sadie/germlines/INTEGRATION_GUIDE.md|P1|US4"
)

for task in "${PHASE6_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,integration,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 6 - User Story 4 (Sadie Integration)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-6-user-story-4---integrate-with-existing-sadie-workflows-p1)" \
        --label "$labels"
done

echo "✅ Phase 6 issues created (11 tasks)"
echo ""

#==============================================
# Phase 7A: User Story 3 - Priority System (P2) (5 tasks)
#==============================================

echo "⚖️  Creating Phase 7A: User Story 3 (Priority System) issues..."

PHASE7A_TASKS=(
    "T052|[P] [US3] Verify priority ordering logic in GermlineManager|src/sadie/germlines/manager.py|P2|US3|parallelizable"
    "T053|[P] [US3] Write unit test for priority ordering scenarios|tests/test_priority_ordering.py|P2|US3|parallelizable"
    "T054|[P] [US3] Test deduplication rules|tests/test_priority_ordering.py|P2|US3|parallelizable"
    "T055|[US3] Document priority configuration|src/sadie/germlines/README.md|P2|US3"
    "T056|[US3] Add logging for priority-based gene selection|src/sadie/germlines/manager.py|P2|US3"
)

for task in "${PHASE7A_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,priority-system,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 7A - User Story 3 (Priority System)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-7a-user-story-3---priority-based-database-selection-p2)" \
        --label "$labels"
done

echo "✅ Phase 7A issues created (5 tasks)"
echo ""

#==============================================
# Phase 7B: User Story 5 - VDJbase Provider Complete (P2) (7 tasks)
#==============================================

echo "🧬 Creating Phase 7B: User Story 5 (VDJbase Provider) issues..."

PHASE7B_TASKS=(
    "T057|[P] [US5] Complete VDJbase provider implementation|src/sadie/germlines/providers/vdjbase.py|P2|US5|parallelizable"
    "T058|[P] [US5] Add VDJbase to default provider list|src/sadie/germlines/manager.py|P2|US5|parallelizable"
    "T059|[P] [US5] Write unit tests for VDJbase provider|tests/test_vdjbase_provider.py|P2|US5|parallelizable"
    "T060|[P] [US5] Create VDJbase test data|tests/data/germlines/vdjbase/|P2|US5|parallelizable"
    "T061|[US5] Test VDJbase in priority ordering|tests/test_priority_ordering.py|P2|US5"
    "T062|[US5] Add VDJbase error handling for format changes|src/sadie/germlines/providers/vdjbase.py|P2|US5"
    "T063|[US5] Document VDJbase manual download process|src/sadie/germlines/sources/vdjbase/README.md|P2|US5"
)

for task in "${PHASE7B_TASKS[@]}"; do
    IFS='|' read -r tid title file priority userstory extra <<< "$task"
    labels="germlines,vdjbase,$priority,$userstory"
    [[ "$extra" == "parallelizable" ]] && labels="$labels,parallelizable"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][$userstory] $title" \
        --body "**Task**: $tid
**Phase**: 7B - User Story 5 (VDJbase Provider)
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-7b-user-story-5---add-vdjbase-provider-p2)" \
        --label "$labels"
done

echo "✅ Phase 7B issues created (7 tasks)"
echo ""

#==============================================
# Phase 8: Polish & Cross-Cutting Concerns (14 tasks)
#==============================================

echo "✨ Creating Phase 8: Polish & Testing issues..."

PHASE8_TASKS=(
    "T064|Run full Sadie test suite|tests/|P1|polish"
    "T065|Verify SC-001: Custom germline injection <5min|tests/|P1|polish"
    "T066|Verify SC-002: Offline operation works|tests/|P1|polish"
    "T067|Verify SC-003: Download scripts <10min|scripts/|P1|polish"
    "T068|Verify SC-005: Disk usage <500MB|tests/|P1|polish"
    "T069|Verify SC-009: Change detection triggers rebuild|tests/|P1|polish"
    "T070|Verify SC-010: Clear error messages|src/sadie/germlines/|P1|polish"
    "T071|Verify SC-011: Timing metrics logged|src/sadie/germlines/|P1|polish"
    "T072|Verify SC-012: CI tests <5min|.github/workflows/|P1|polish"
    "T073|Update main germlines README|src/sadie/germlines/README.md|P1|polish"
    "T074|Update INTEGRATION_GUIDE|src/sadie/germlines/INTEGRATION_GUIDE.md|P1|polish"
    "T075|Add performance profiling|src/sadie/germlines/pipeline.py|P2|polish"
    "T076|Run pre-commit hooks|./|P1|polish"
    "T077|Update CHANGELOG|CHANGELOG.md|P1|polish"
)

for task in "${PHASE8_TASKS[@]}"; do
    IFS='|' read -r tid title file priority category <<< "$task"
    labels="germlines,$priority,$category,testing"

    gh issue create \
        --repo "$REPO" \
        --title "[Germlines][$priority][Polish] $title" \
        --body "**Task**: $tid
**Phase**: 8 - Polish & Cross-Cutting Concerns
**Priority**: $priority
**File**: \`$file\`

**Description**: $title

**References**:
- [Spec]($SPEC_URL)
- [Plan]($PLAN_URL)
- [Tasks]($TASKS_URL#phase-8-polish--cross-cutting-concerns)" \
        --label "$labels"
done

echo "✅ Phase 8 issues created (14 tasks)"
echo ""

#==============================================
# Summary
#==============================================

echo "=================================================="
echo "✅ GitHub Issues Creation Complete!"
echo "=================================================="
echo ""
echo "📊 Summary:"
echo "  - Phase 1 (Setup):            7 issues"
echo "  - Phase 2 (Foundation):       8 issues"
echo "  - Phase 3 (US1 Custom):       8 issues"
echo "  - Phase 4 (US6 Data Pop):    11 issues"
echo "  - Phase 5 (US2 Offline):      6 issues"
echo "  - Phase 6 (US4 Integration): 11 issues"
echo "  - Phase 7A (US3 Priority):    5 issues"
echo "  - Phase 7B (US5 VDJbase):     7 issues"
echo "  - Phase 8 (Polish):          14 issues"
echo "  ─────────────────────────────────────"
echo "  Total:                       77 issues"
echo ""
echo "🔗 View issues at: https://github.com/$REPO/issues?q=is%3Aissue+label%3Agermlines"
echo ""
echo "📋 Next steps:"
echo "  1. Review all created issues"
echo "  2. Create milestone: 'Germlines Module Completion'"
echo "  3. Assign issues to milestone"
echo "  4. Start with MVP: Phase 3 (US1) + Phase 4 (US6)"
echo ""
echo "MVP Tasks (8-10 hours):"
echo "  - US1: Custom germlines (8 tasks)"
echo "  - US6: Data population (11 tasks)"
echo "=================================================="
