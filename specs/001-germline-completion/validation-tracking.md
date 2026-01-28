# Validation Period Tracking

## Overview

This document tracks the validation period for the germlines module migration from G3 API.

## Timeline

| Milestone | Date | Status |
|-----------|------|--------|
| Validation Start | 2026-01-18 | Active |
| Release 1 | - | Pending |
| Release 2 | - | Pending |
| Release 3 | - | Pending |
| Deprecation Notice | - | Pending |
| G3 Removal | - | Pending |

## Success Criteria (FR-017b)

- [ ] 3 releases without critical bugs
- [ ] Zero data loss incidents
- [ ] Performance parity with G3 API
- [ ] User feedback collection complete

## Release Tracking

### Release 1

- Version: -
- Date: -
- Critical bugs: -
- Performance notes: -

### Release 2

- Version: -
- Date: -
- Critical bugs: -
- Performance notes: -

### Release 3

- Version: -
- Date: -
- Critical bugs: -
- Performance notes: -

## Bug Tracker

| ID | Severity | Description | Status | Resolution |
|----|----------|-------------|--------|------------|
| - | - | - | - | - |

## Performance Baseline

### G3 API Baseline

| Operation | Time (ms) | Memory (MB) |
|-----------|-----------|-------------|
| Gene lookup | - | - |
| Full rebuild | - | - |
| IgBLAST init | - | - |

### Germlines Module

| Operation | Time (ms) | Memory (MB) | Delta |
|-----------|-----------|-------------|-------|
| Gene lookup | - | - | - |
| Full rebuild | - | - | - |
| IgBLAST init | - | - | - |

## Deprecation Schedule

Per FR-019a:

1. CHANGELOG entry announcing deprecation
2. GitHub discussion for feedback
3. Runtime warning when G3 used (SADIE_USE_GERMLINES_MODULE=false)

## Notes

- Validation period: Minimum 3 releases or 3 months
- Feature flag remains available during validation
- G3 removal requires explicit approval
