# Implementation Plan: Germlines Module Completion

**Branch**: `001-germline-completion` | **Date**: 2026-01-08 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/001-germline-completion/spec.md`

## Summary

Complete the germlines module by implementing VDJbase provider, populating reference data sources (IMGT/OGRDB/VDJbase), integrating with existing Sadie components (IgBLAST, Reference, HMM), and removing G3 API dependency. The implementation uses a phased migration approach with feature flag, uses Biopython alignment-based auto-gapping, and employs a provider-based architecture for multi-source germline database management with priority-based merging.

## Technical Context

**Language/Version**: Python 3.11 (existing Sadie requirement)
**Primary Dependencies**:
- BioPython (FASTA parsing, SeqIO, pairwise alignment)
- HMMER (existing Sadie dependency for hmmbuild; not used for gapping)
- subprocess (makeblastdb for BLAST database generation)
- Python logging module (structured logging)
- hashlib (change detection via file hashing)

**Storage**: Local filesystem (FASTA files, BLAST databases, auxiliary files)
**Testing**: pytest (existing Sadie test framework) with curated test dataset
**Target Platform**: Cross-platform (Linux, macOS) - matches existing Sadie support
**Project Type**: Python library module within larger Sadie package
**Performance Goals**:
- Database rebuild <2 minutes for human reference data
- Change detection <1 second
- Test suite <5 minutes in CI/CD

**Constraints**:
- Offline-capable after initial setup (no runtime network calls)
- <500MB disk space for human reference data (IMGT + OGRDB + VDJbase)
- Backward compatible with existing Sadie API
- Must reuse existing Sadie infrastructure (logging patterns, path conventions)

**Scale/Scope**:
- ~450 IGHV genes from IMGT
- OGRDB data: Human (~47 genes), Mouse (~1,771 genes including 1,230 IGKV)
- VDJbase population-specific alleles (variable, API-based download)
- 4 providers (custom, IMGT, OGRDB, VDJbase)
- 3 integration points (IgBLAST, Reference, HMM)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Principle I: Provider-Based Architecture ✓
- **Requirement**: Each data source as independent provider
- **Implementation**: VDJbaseProvider follows GermlineProvider base class
- **Status**: PASS - Follows existing pattern (CustomProvider, IMGTProvider, OGRDBProvider already implemented)

### Principle II: Priority-Based Merging (NON-NEGOTIABLE) ✓
- **Requirement**: Explicit priority ordering with conflict resolution
- **Implementation**: GermlineManager handles default priority `custom > ogrdb > vdjbase > imgt` (user-configurable)
- **Status**: PASS - Deduplication logic already implemented in manager.py:145-285

### Principle III: Local-First Operation ✓
- **Requirement**: No runtime external API dependencies
- **Implementation**: Feature flag explicitly selects G3 or germlines path; no automatic fallback during validation
- **Status**: PASS - VDJbase uses manual FASTA files, IMGT/OGRDB automated download stores locally

### Principle IV: Staged Pipeline Architecture ✓
- **Requirement**: Three-stage processing (sources → normalized → igblast)
- **Implementation**: Existing pipeline.py implements staged architecture
- **Status**: PASS - VDJbase integrates into existing pipeline stages

### Principle V: Integration Compatibility ✓
- **Requirement**: Backward compatibility with existing Sadie components
- **Implementation**: Adapter pattern maintains G3 API response format; feature flag enables gradual migration without automatic fallback
- **Status**: PASS - Phased migration strategy with explicit feature flag minimizes breaking changes

**GATE RESULT: PASS** - All constitutional principles satisfied. No complexity justification required.

## Project Structure

### Documentation (this feature)

```text
specs/001-germline-completion/
├── spec.md              # Feature specification (completed)
├── plan.md              # This file
├── research.md          # Phase 0 output (to be created)
├── data-model.md        # Phase 1 output (to be created)
├── quickstart.md        # Phase 1 output (to be created)
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code (repository root)

```text
src/sadie/germlines/
├── __init__.py                           # Public API (exists)
├── models.py                             # GermlineGene model (exists)
├── manager.py                            # Priority-based manager (exists)
├── pipeline.py                           # Staged processing (exists)
│
├── providers/
│   ├── base.py                          # Provider interface (exists)
│   ├── custom.py                        # Custom provider (exists)
│   ├── imgt.py                          # IMGT provider (exists)
│   ├── ogrdb.py                         # OGRDB provider (exists)
│   └── vdjbase.py                       # VDJbase provider (NEW - to implement)
│
├── builders/
│   ├── blast.py                         # BLAST DB builder (exists)
│   ├── aux.py                           # Auxiliary file builder (exists)
│   └── gapper.py                        # Auto-gapping using BioPython alignment (NEW - to implement)
│
├── sources/
│   ├── custom/                          # Custom sequences (exists)
│   ├── imgt/                            # IMGT data (exists, needs population)
│   ├── ogrdb/                           # OGRDB data (exists, needs population)
│   └── vdjbase/                         # VDJbase data (NEW - directory to create)
│       ├── README.md                    # Manual download instructions (NEW)
│       └── human/                       # Human VDJbase data (NEW)
│
├── normalized/                           # Processed data (Stage 2)
│   └── {species}/
│       ├── gapped/                      # IMGT-gapped sequences
│       └── ungapped/                    # Ungapped sequences
│
├── igblast/                              # IgBLAST format (Stage 3)
│   ├── database/
│   │   └── {species}/
│   ├── aux_db/
│   └── internal_data/
│
├── scripts/
│   ├── download_imgt.py                 # IMGT downloader (exists, needs completion)
│   ├── download_ogrdb.py                # OGRDB downloader (NEW - to implement)
│   └── validate.py                      # Validation script (NEW - to implement)
│
└── tests/
    ├── test_manager.py                  # Manager tests (exists, minimal)
    ├── test_vdjbase_provider.py         # VDJbase tests (NEW - to implement)
    ├── test_gapper.py                   # Gapping tests (NEW - to implement)
    ├── test_integration.py              # Full pipeline tests (NEW - to implement)
    └── data/
        ├── custom/                      # Test dataset (NEW - to create)
        ├── imgt/
        ├── ogrdb/
        └── vdjbase/

src/sadie/airr/igblast/
└── germline.py                          # IgBLAST integration (MODIFY - update paths)

src/sadie/reference/
└── reference.py                         # Reference system (MODIFY - use germlines module)

src/sadie/renumbering/aligners/
└── hmmer.py                             # HMM builder (MODIFY - use germlines module)
```

**Structure Decision**: Single project structure. Germlines module is a sub-package within existing `src/sadie/` hierarchy. This maintains consistency with Sadie's existing architecture and simplifies imports. No new top-level directories needed.

## Complexity Tracking

> **No violations** - Constitution check passed completely.

## Phase 0: Research & Unknowns Resolution

**Prerequisites**: Constitution check PASSED

### Research Questions

1. **VDJbase Data Format & Access**
   - Question: What is VDJbase's current FASTA format and download mechanism?
   - Research needed: VDJbase website/API documentation, example files
   - Impact: VDJbase provider implementation

2. **Biopython Alignment Strategy for Gapping**
   - Decision: Use per-gene IMGT-gapped templates when available; fallback to per-segment consensus template
   - Template source: IMGT provider gapped FASTA in `src/sadie/germlines/sources/imgt/`; consensus derived from those gapped sequences
   - Research needed: Identify available gapped references in IMGT sources; evaluate Bio.Align/PairwiseAligner settings for AA alignment
   - Impact: Gapper module implementation

3. **G3 API Response Format**
   - Question: What exact JSON format does G3 return for gene queries?
   - Research needed: Inspect existing G3 client code, sample responses
   - Impact: Reference system adapter design

4. **IgBLAST Auxiliary File Format**
   - Question: What CDR/FWR annotations are required in auxiliary files?
   - Research needed: IgBLAST documentation, existing Sadie aux files
   - Impact: Auxiliary file builder completion

5. **IMGT/OGRDB Download URLs**
   - Question: Current stable URLs for bulk FASTA downloads?
   - Research needed: IMGT/OGRDB websites, download APIs
   - Impact: Download script implementation

**Output**: `research.md` with findings and decisions

## Phase 1: Design & Contracts

**Prerequisites**: `research.md` complete

### Data Model

**File**: `data-model.md`

**Entities** (from spec):

1. **VDJbaseProvider**
   - Inherits: GermlineProvider
   - Fields: data_dir (Path), name (str = "vdjbase")
   - Methods: fetch_genes(), fetch_gene_by_name(), is_available(), get_metadata()
   - Relationships: Produces GermlineGene objects

2. **GermlineGene** (exists, document)
   - Fields: name, species, segment, chain, sequence, sequence_gapped, is_functional, functionality, source, regions (dict)
   - Validation: Valid nucleotides (ACGT), non-empty name
   - State: Immutable dataclass

3. **GapperService** (new)
   - Methods: gap_sequence(sequence: str, segment: str) → str
   - Dependencies: BioPython alignment utilities (pairwise AA alignment)
   - Error handling: Returns original if gapping fails, logs warning

4. **FeatureFlag** (new)
   - Fields: SADIE_USE_GERMLINES_MODULE (env var, default=True)
   - Methods: use_germlines_module() → bool
   - Purpose: Control migration strategy

### API Contracts

**File**: `contracts/`

Since this is a Python library module (not REST/GraphQL API), contracts are Python interface signatures:

**Contract 1: Provider Interface** (`contracts/provider_interface.py`)
```python
from abc import ABC, abstractmethod
from typing import List, Optional
from pathlib import Path

class GermlineProvider(ABC):
    """Base interface all providers must implement."""

    @abstractmethod
    def fetch_genes(self, species: str, segment: str, chain: str) -> List[GermlineGene]:
        """Fetch genes for species/segment/chain."""
        pass

    @abstractmethod
    def fetch_gene_by_name(self, name: str, species: str) -> Optional[GermlineGene]:
        """Fetch specific gene by name."""
        pass

    @abstractmethod
    def is_available(self, species: str) -> bool:
        """Check if data available for species."""
        pass

    @abstractmethod
    def get_metadata(self) -> ProviderMetadata:
        """Get provider metadata."""
        pass
```

**Contract 2: Gapper Interface** (`contracts/gapper_interface.py`)
```python
from typing import Optional

class GapperInterface:
    """Interface for sequence gapping service."""

    def gap_sequence(self, sequence: str, segment: str, scheme: str = "imgt") -> str:
        """
        Gap sequence using specified numbering scheme.

        Args:
            sequence: Ungapped nucleotide sequence
            segment: V, D, or J
            scheme: Numbering scheme (default: imgt)

        Returns:
            Gapped sequence with dots/dashes at appropriate positions

        Raises:
            ValueError: If sequence invalid or segment unsupported
        """
        pass
```

**Contract 3: Migration Adapter** (`contracts/g3_adapter_interface.py`)
```python
from typing import Dict, Any

class G3AdapterInterface:
    """Interface for G3 API compatibility layer."""

    def get_gene(self, gene_name: str, species: str, source: str) -> Dict[str, Any]:
        """
        Get gene in G3 API response format.

        Returns:
            {
                'gene': str,
                'sequence': str,
                'sequence_gapped': str,
                'species': str,
                'segment': str,
                'chain': str,
                'source': str,
                'functional': bool,
                'regions': {
                    'fwr1': str,
                    'cdr1': str,
                    'fwr2': str,
                    'cdr2': str,
                    'fwr3': str
                }
            }
        """
        pass
```

### Quickstart Guide

**File**: `quickstart.md`

```markdown
# Germlines Module Development Quickstart

## Prerequisites

- Python 3.11+
- Sadie development environment set up
- Poetry for dependency management

## Development Setup

1. **Clone and setup**:
   ```bash
   cd /path/to/sadie
   git checkout 001-germline-completion
   poetry install
   ```

2. **Populate test data**:
   ```bash
  # Test dataset already in repo: src/sadie/germlines/tests/data/
   pytest src/sadie/germlines/tests/ --co  # Verify tests discovered
   ```

3. **Run tests**:
   ```bash
   poetry run pytest src/sadie/germlines/tests/ -v
   ```

## Implementation Order

1. **VDJbase Provider** (2-3 hours)
   - File: `src/sadie/germlines/providers/vdjbase.py`
   - Pattern: Copy ogrdb.py structure, modify for VDJbase format
  - Test: `src/sadie/germlines/tests/test_vdjbase_provider.py`

2. **Auto-Gapping Service** (2-3 hours)
   - File: `src/sadie/germlines/builders/gapper.py`
   - Integration: Use BioPython AA alignment against IMGT-gapped templates (per-gene fallback to per-segment consensus) to derive gap positions
  - Test: `src/sadie/germlines/tests/test_gapper.py`

3. **Download Scripts** (3-4 hours)
   - Files: `src/sadie/germlines/scripts/download_ogrdb.py`, complete `src/sadie/germlines/scripts/download_imgt.py`
   - Features: Species filtering, resume, validation
   - Test: Manual testing + validation script

4. **Integration Updates** (4-5 hours)
   - Files: `src/sadie/airr/igblast/germline.py`, `src/sadie/reference/reference.py`, `src/sadie/renumbering/aligners/hmmer.py`
   - Pattern: Add feature flag check, update paths, maintain backward compatibility
   - Test: Existing Sadie test suite (regression)

5. **Feature Flag & Migration** (2-3 hours)
   - Files: Environment variable handling, adapter pattern
   - Testing: Both modes (germlines and G3 via feature flag)

## Key Patterns

### Provider Implementation

```python
class VDJbaseProvider(GermlineProvider):
    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.name = "vdjbase"

    def fetch_genes(self, species: str, segment: str, chain: str) -> List[GermlineGene]:
        fasta_path = self.get_fasta_path(species, segment, chain)
        if not fasta_path.exists():
            logger.debug(f"No VDJbase file: {fasta_path}")
            return []
        return self._parse_fasta(fasta_path, species, segment, chain)
```

### Feature Flag Usage

```python
import os

def use_germlines_module() -> bool:
    return os.getenv("SADIE_USE_GERMLINES_MODULE", "true").lower() == "true"

# In Reference system:
if use_germlines_module():
    from sadie.germlines import get_gene_by_name
    gene = get_gene_by_name(name, species)
else:
    from sadie.renumbering.clients.g3 import G3
    gene = G3().get_gene(name, species)
```

## Testing

### Unit Tests
```bash
pytest src/sadie/germlines/tests/test_vdjbase_provider.py -v
```

### Integration Tests
```bash
pytest src/sadie/germlines/tests/test_integration.py -v
```

### Regression Tests
```bash
# Run full Sadie test suite
pytest tests/ -k igblast
```

## Debugging

Enable debug logging:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

Check provider loading:
```python
from sadie.germlines import GermlineManager
manager = GermlineManager()
for provider in manager.providers:
    print(f"{provider.name}: {provider.is_available('human')}")
```

## Common Issues

1. **FASTA not found**: Check `src/sadie/germlines/sources/{provider}/{species}/` exists
2. **Gapping fails**: Verify IMGT gapped reference exists and alignment settings are valid
3. **Tests fail**: Ensure test dataset populated in `src/sadie/germlines/tests/data/`
```

### Agent Context Update

**Action**: Run `.specify/scripts/bash/update-agent-context.sh claude`

**Content to add** (new technologies from this plan):
- VDJbase FASTA format parsing
- BioPython alignment-based gapping
- Feature flag pattern for phased migration
- Curated test dataset approach

## Post-Design Constitution Re-Check

*Re-evaluate after Phase 1 design complete*

### Principle I: Provider-Based Architecture ✓
- Design: VDJbaseProvider follows GermlineProvider interface exactly
- Status: PASS

### Principle II: Priority-Based Merging ✓
- Design: No changes to priority system; VDJbase added as 4th provider; default order `custom > ogrdb > vdjbase > imgt`
- Status: PASS

### Principle III: Local-First Operation ✓
- Design: VDJbase uses manual FASTA files; feature flag explicitly selects G3 during validation (no automatic fallback)
- Status: PASS

### Principle IV: Staged Pipeline Architecture ✓
- Design: VDJbase integrates into existing pipeline stages; gapper fits in normalized stage
- Status: PASS

### Principle V: Integration Compatibility ✓
- Design: Feature flag + adapter pattern maintains backward compatibility
- Status: PASS

**GATE RESULT: PASS** - All principles maintained through design phase.

## Implementation Phases Summary

### Phase 0: Research ✓ (plan created)
- Outputs: `research.md`
- Estimated: 2-3 hours

### Phase 1: Design ✓ (plan created)
- Outputs: `data-model.md`, `contracts/`, `quickstart.md`
- Estimated: 1-2 hours

### Phase 2: Task Breakdown (next: `/speckit.tasks`)
- Outputs: `tasks.md`
- Estimated: Implementation begins after task creation

## Risk Assessment

| Risk | Impact | Mitigation |
|------|---------|------------|
| Biopython gapping accuracy | High | Align to IMGT-gapped templates per gene (fallback to per-segment consensus); fallback to ungapped on failure |
| VDJbase format changes | Medium | Manual FASTA approach; well-documented README |
| G3 parity during migration | High | Feature flag + extensive regression testing |
| Test data completeness | Medium | Curated test dataset covers edge cases; document limitations |
| Performance degradation | Low | Profile critical paths; 2-minute rebuild budget |

## Success Metrics Tracking

From spec success criteria (SC-001 to SC-015):

- **SC-001**: Custom germline in <5 minutes → Test with stopwatch during implementation
- **SC-002**: Offline operation → Integration test with network disabled
- **SC-003**: Download scripts <10 minutes → Time during manual testing
- **SC-004**: 100% test pass rate → CI/CD verification
- **SC-005**: <500MB disk usage → du -sh after population
- **SC-006**: Priority ordering → Unit test comparing provider orders
- **SC-007**: Feature flag works → Test both true/false modes
- **SC-008**: Setup in <30 minutes → User testing with fresh environment
- **SC-009**: Change detection → Unit test file modification triggers
- **SC-010**: Clear error messages → Code review + error handling tests
- **SC-011**: Timing metrics logged → Verify logs contain duration fields
- **SC-012**: CI tests <5 minutes → GitHub Actions timing
- **SC-013-015**: Qualitative → Post-implementation user testing + code review

## Next Steps

1. Execute `/speckit.plan` Phase 0: Create `research.md`
2. Execute `/speckit.plan` Phase 1: Create design artifacts
3. Execute `/speckit.tasks` to generate task breakdown
4. Begin implementation following task sequence
5. Run regression tests continuously
6. Update documentation as implementation progresses

**Plan Status**: Ready for Phase 0 execution
**Estimated Total Implementation Time**: 24-28 hours (VDJbase 3h + Gapping 3h + Downloads 4h + Integration 5h + Testing 5h + Documentation 2h + Buffer 6h)
