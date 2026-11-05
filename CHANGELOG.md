# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **OGRDB germline database updated** to Zenodo record `23947887` (published 2026-11-04).
  Version is now tracked by Zenodo record ID (`release-23947887`) and auto-checked on each
  `sadie germlines populate` run.

  **Human** — constant-region isotype genes added for the first time:

  | Locus | Sequences |
  |-------|-----------|
  | IGHV  | 229       |
  | IGKV  | 83        |
  | IGLV  | 88        |
  | IGHD  | 53        |
  | IGHJ  | 7         |
  | IGKJ  | 7         |
  | IGLJ  | 10        |
  | IGHC  | 328 *(IGHA×102, IGHE×36, IGHG×152, IGHM×28, IGHD-delta×10)* |

  **Mouse** — strain-specific alleles now fully populated across 45 germline sets covering
  18 inbred strains (AKR/J, A/J, 129S1/SvImJ, BALB/c/ByJ, C3H/HeJ, C57BL/6J, CAST/EiJ,
  CBA/J, DBA/1J, DBA/2J, LEWES/EiJ, MRL/MpJ, MSM/MsJ, NOD/ShiLtJ, NOR/LtJ, NZB/BlNJ,
  PWD/PhJ, SJL/J):

  | Locus | Sequences |
  |-------|-----------|
  | IGHV  | 1,117     |
  | IGKV  | 2,810     |
  | IGLV  | 66        |
  | IGHD  | 132       |
  | IGHJ  | 56        |
  | IGKJ  | 22        |
  | IGLJ  | 14        |

- **Custom mouse provider** — 29 strain-specific S-alleles (`IGHV1S*`, `IGKV6-B/C/D*01`,
  etc.) extracted from bundled G3 reference data and registered under
  `sources/custom/mouse/`; `reference.yml` mouse section restored to 635 genes.

### Changed

- **`sadie germlines populate`** — OGRDB now auto-detects new Zenodo releases at run time;
  `is_up_to_date()` compares Zenodo record IDs rather than calendar month.

### Added
- **MOTIF_PROVENANCE API**: New provenance metadata system for J-gene FWR4 motif patterns
  - Access pattern origins and validation status via `sadie.reference.MOTIF_PROVENANCE`
  - Distinguish between IMGT-validated (2 species) and legacy patterns (35 species)
  - Comprehensive metadata including source, last_reviewed date, and validation status

- **JSON-based MOTIF_LOOKUP**: Migrated from Python dict to JSON with enhanced capabilities
  - Centralized motif registry in `src/sadie/reference/data/j_gene_motif.json`
  - Schema validation with `specs/001-germline-completion/contracts/motif-registry.schema.json`
  - Backward compatibility maintained - existing code continues to work unchanged

- **Empirical Coverage Testing**: Framework for validating motif patterns against real sequences
  - Test suite in `tests/unit/reference/test_motif_coverage.py`
  - Species-specific thresholds: 95% for IMGT-validated, 80% for legacy patterns
  - Biopython integration for sequence validation
  - JSON schema-compliant coverage reports

- **Migration Validation Tools**: Comprehensive validation and migration tools
  - Migration validation script: `scripts/validate_motif_migration.py`
  - Schema validation for all 37 species
  - Performance benchmarking and regression testing
  - Integration testing with existing j_gene_data module

- **Enhanced Error Handling**: Improved error messages and data integrity validation
  - Helpful error messages for missing species with suggestions
  - Comprehensive JSON validation and structure checks
  - Graceful handling of malformed data
  - Performance-optimized loading with caching

### Changed
- **j_gene_data.py**: Updated to use centralized motif lookup from `sadie.reference`
  - Removed local JSON loading code
  - Improved integration with new provenance system
  - Maintained all existing functionality

### Technical Details
- **Species Coverage**: Support for 37 species with comprehensive motif patterns
- **Performance**: <50ms loading time, cached for subsequent access
- **Compatibility**: Full backward compatibility via `_MotifLookupProxy` class
- **Validation**: 17 test classes with comprehensive coverage validation
- **Documentation**: Updated quickstart guide with migration examples

### Migration Guide

#### For Existing Users (No Changes Required)
```python
# This code continues to work unchanged
from sadie.reference import MOTIF_LOOKUP

pattern = MOTIF_LOOKUP['human']['IGHJ']  # Still works exactly the same
```

#### For New Features
```python
# New provenance API
from sadie.reference import MOTIF_LOOKUP, MOTIF_PROVENANCE

# Check if patterns are validated
human_prov = MOTIF_PROVENANCE['human']
if human_prov['imgt_validated']:
    pattern = MOTIF_LOOKUP['human']['IGHJ']  # High confidence
    print(f"IMGT-validated pattern from: {human_prov['source']}")
```

#### Validation Tools
```bash
# Validate complete migration
python scripts/validate_motif_migration.py --verify-all

# Run empirical coverage tests
pytest tests/unit/reference/test_motif_coverage.py -v
```

### Compatibility
- ✅ **Python 3.10+**: Full compatibility maintained
- ✅ **Existing API**: No breaking changes to `MOTIF_LOOKUP` usage
- ✅ **j_gene_data.py**: Seamless integration with updated module
- ✅ **Performance**: <50ms loading requirement met
- ✅ **Testing**: All existing tests continue to pass

### Dependencies
- **BioPython**: Required for empirical coverage testing (optional feature)
- **jsonschema**: Required for schema validation (development/testing)
- **pytest**: Required for running test suite

---

## Development Notes

This release implements the complete MOTIF_LOOKUP provenance and coverage testing system as specified in the project requirements. The implementation maintains full backward compatibility while adding comprehensive metadata tracking and empirical validation capabilities.

Key architectural decisions:
1. **JSON-first approach**: Central registry with schema validation
2. **Backward compatibility**: Proxy pattern preserves existing API
3. **Performance optimization**: Lazy loading with caching
4. **Comprehensive testing**: 60+ test cases covering all aspects
5. **Production readiness**: Full error handling and validation

The system supports 37 species with over 130 motif patterns, providing a solid foundation for J-gene FWR4 boundary detection across diverse biological systems.