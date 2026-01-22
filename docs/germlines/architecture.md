# Architecture: Technical Design

This document provides technical details about the germlines module architecture for developers contributing to SADIE or building tools that integrate with the germlines system.

---

## Overview

The germlines module is a staged pipeline that downloads, processes, merges, and builds germline reference databases for AIRR annotation. It replaces the legacy G3 API with a local, multi-source, extensible system.

**Key Design Principles:**

1. **Provider Abstraction**: Generic provider interface supports multiple data sources
2. **Staged Pipeline**: Clear separation of download, merge, build, and validation
3. **Priority-Based Merging**: Deterministic conflict resolution
4. **Checkpoint Recovery**: Automatic resume after failures
5. **Integration Layer**: Clean integration with AIRR and renumbering modules

---

## System Architecture

### High-Level Pipeline

```
┌─────────────────────────────────────────────────────────────┐
│                     User Command                             │
│              sadie germlines populate                        │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                  Stage 1: Download                           │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐   │
│  │   IMGT   │  │  OGRDB   │  │ VDJbase  │  │  Custom  │   │
│  │ Provider │  │ Provider │  │ Provider │  │ Provider │   │
│  └────┬─────┘  └────┬─────┘  └────┬─────┘  └────┬─────┘   │
│       │             │              │             │          │
│       ▼             ▼              ▼             ▼          │
│  ~/.sadie/germlines/sources/                                │
│       imgt/       ogrdb/       vdjbase/      custom/        │
└─────────────────────────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                  Stage 2: Merge                              │
│  Priority-based deduplication:                               │
│  Custom > OGRDB > VDJbase > IMGT                            │
│                                                              │
│  Output: Merged FASTA per species/locus                     │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                  Stage 3: Build                              │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐        │
│  │   BLAST     │  │  Auxiliary  │  │  Internal   │        │
│  │  Database   │  │   Files     │  │    Data     │        │
│  │ (makeblastdb)│  │   (.aux)    │  │ (.ndm.imgt) │        │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘        │
│         │                │                 │                │
│         ▼                ▼                 ▼                │
│  ~/.sadie/germlines/databases/                              │
│       human/         mouse/         ...                     │
│        ├── IGHV.ndb  ├── IGHV.ndb                          │
│        ├── IGHV.aux  ├── IGHV.aux                          │
│        └── ...       └── ...                                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│                  Stage 4: Validate                           │
│  - Database integrity checks                                 │
│  - File completeness verification                            │
│  - Version metadata tracking                                 │
└─────────────────────────────────────────────────────────────┘
```

---

## Component Design

### 1. Provider Interface

**Abstract Base Class: `GermlineProvider`**

```python
from abc import ABC, abstractmethod
from typing import List, Dict
from pathlib import Path

class GermlineProvider(ABC):
    """Abstract base class for germline data providers."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Provider name (e.g., 'imgt', 'ogrdb')."""
        pass

    @abstractmethod
    def get_available_species(self) -> List[str]:
        """Return list of available species."""
        pass

    @abstractmethod
    def download(self, species: str, output_dir: Path) -> bool:
        """Download germline data for species to output_dir."""
        pass

    @abstractmethod
    def get_version(self) -> str:
        """Return current version/release identifier."""
        pass

    @abstractmethod
    def needs_update(self, local_version: str) -> bool:
        """Check if local version needs update."""
        pass
```

**Concrete Implementations:**

- `IMGTProvider` - Downloads from IMGT FTP server
- `OGRDBProvider` - Downloads from OGRDB API/repository
- `VDJbaseProvider` - Downloads from VDJbase API
- `CustomProvider` - Reads from local filesystem

**Provider Registry:**

```python
PROVIDERS = {
    'imgt': IMGTProvider(),
    'ogrdb': OGRDBProvider(),
    'vdjbase': VDJbaseProvider(),
    'custom': CustomProvider(),
}
```

---

### 2. Download Stage

**Component: `GermlineDownloader`**

**Responsibilities:**
- Coordinate downloads across multiple providers
- Manage checkpoint files for resume capability
- Validate downloaded FASTA files
- Track provider versions

**Key Methods:**

```python
class GermlineDownloader:
    def download_provider(
        self,
        provider: GermlineProvider,
        species: List[str],
        force: bool = False
    ) -> Dict[str, bool]:
        """
        Download germlines from provider for specified species.

        Args:
            provider: Provider instance
            species: List of species names
            force: Force re-download even if up-to-date

        Returns:
            Dict mapping species -> success status
        """
        pass

    def resume_from_checkpoint(self) -> Dict[str, Any]:
        """Load checkpoint file and return resume state."""
        pass

    def save_checkpoint(self, state: Dict[str, Any]):
        """Save current state to checkpoint file."""
        pass
```

**Checkpoint File Format:**

```json
{
  "provider": "imgt",
  "total_species": 33,
  "completed_species": ["human", "mouse"],
  "failed_species": [],
  "current_species": "rat",
  "timestamp": "2026-01-22T10:30:00Z"
}
```

---

### 3. Merge Stage

**Component: `GermlineMerger`**

**Responsibilities:**
- Load FASTA files from all providers
- Apply priority-based deduplication
- Normalize sequence headers
- Output merged FASTA per species/locus

**Priority Algorithm:**

```python
def merge_sequences(
    providers: List[str],
    species: str,
    locus: str
) -> Dict[str, str]:
    """
    Merge sequences from multiple providers.

    Priority order: custom > ogrdb > vdjbase > imgt

    Args:
        providers: List of provider names
        species: Species name
        locus: Locus (IGHV, IGHD, IGHJ, etc.)

    Returns:
        Dict mapping sequence_name -> sequence
    """
    merged = {}

    # Process in reverse priority order (lowest first)
    for provider in ['imgt', 'vdjbase', 'ogrdb', 'custom']:
        if provider not in providers:
            continue

        fasta_path = get_provider_fasta(provider, species, locus)
        sequences = load_fasta(fasta_path)

        for name, seq in sequences.items():
            # Later providers override earlier ones
            merged[name] = seq

    return merged
```

**Deduplication Strategy:**

- **Exact name match**: Sequences with identical names are considered duplicates
- **No sequence comparison**: Sequence content is NOT compared for deduplication
- **Last-wins**: Higher priority provider overwrites lower priority

**Example:**

```
IMGT:   IGHV3-30*01 → CAGGTGCAGCT...  (495bp)
OGRDB:  IGHV3-30*01 → CAGGTGCAGCT...  (500bp)
Custom: IGHV3-30*01 → CAGGTGCAGCT...  (500bp)

Result: IGHV3-30*01 → CAGGTGCAGCT...  (500bp from Custom)
```

---

### 4. Build Stage

**Component: `DatabaseBuilder`**

Three sub-builders work in sequence:

#### 4.1 BLAST Database Builder

```python
class BlastDatabaseBuilder:
    """Build NCBI BLAST databases using makeblastdb."""

    def build(self, fasta_path: Path, output_dir: Path) -> bool:
        """
        Build BLAST database from FASTA file.

        Executes:
            makeblastdb -parse_seqids \
                        -dbtype nucl \
                        -in <fasta_path> \
                        -out <output_dir>/<locus>

        Produces:
            - IGHV.ndb (database)
            - IGHV.nhr (headers)
            - IGHV.nin (index)
            - IGHV.nsq (sequences)
        """
        pass
```

#### 4.2 Auxiliary File Builder

```python
class AuxiliaryFileBuilder:
    """Build .aux files for CDR3 reconstruction."""

    def build(self, fasta_path: Path, output_path: Path) -> bool:
        """
        Generate auxiliary file with CDR/FWR annotations.

        Format:
            >IGHV1-2*01
            7       23
            39      57
            66      104
            110     118

        Where each line pair is: start, end positions for:
            FWR1, CDR1, FWR2, CDR2, FWR3
        """
        pass
```

#### 4.3 Internal Data Builder

```python
class InternalDataBuilder:
    """Build IgBLAST internal_data files."""

    def build(self, fasta_path: Path, output_path: Path) -> bool:
        """
        Generate .ndm.imgt file for IgBLAST.

        Format: Binary file containing:
            - Gene names
            - Sequence lengths
            - CDR3 start positions
        """
        pass
```

**Build Pipeline Orchestration:**

```python
def build_databases(species: str, merged_dir: Path, output_dir: Path):
    """Build all database files for a species."""

    blast_builder = BlastDatabaseBuilder()
    aux_builder = AuxiliaryFileBuilder()
    internal_builder = InternalDataBuilder()

    for locus in ['IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ']:
        fasta_path = merged_dir / f"{locus}.fasta"

        if not fasta_path.exists():
            continue

        # Build BLAST database
        blast_builder.build(fasta_path, output_dir)

        # Build auxiliary files (V genes only)
        if locus.endswith('V'):
            aux_path = output_dir / f"{locus}.aux"
            aux_builder.build(fasta_path, aux_path)

            # Build internal data
            internal_path = output_dir / f"{locus}.ndm.imgt"
            internal_builder.build(fasta_path, internal_path)
```

---

### 5. Validation Stage

**Component: `DatabaseValidator`**

```python
class DatabaseValidator:
    """Validate database integrity and completeness."""

    def validate(self, species: str, db_dir: Path) -> List[str]:
        """
        Validate database for species.

        Checks:
            - All expected files exist (.ndb, .aux, .ndm.imgt)
            - BLAST database is queryable
            - Auxiliary files are parseable
            - Sequence counts match between sources and databases

        Returns:
            List of validation errors (empty if all checks pass)
        """
        errors = []

        # Check BLAST database
        if not self._validate_blast_db(db_dir / "IGHV.ndb"):
            errors.append("Invalid BLAST database: IGHV")

        # Check auxiliary files
        if not self._validate_aux_file(db_dir / "IGHV.aux"):
            errors.append("Invalid auxiliary file: IGHV.aux")

        # Check sequence counts
        if not self._validate_counts(species, db_dir):
            errors.append("Sequence count mismatch")

        return errors
```

---

## Integration Points

### Integration with AIRR Module

**Interface: `GermlineLocator`**

```python
class GermlineLocator:
    """Locate germline databases for AIRR annotation."""

    def get_database_path(self, species: str, locus: str) -> Path:
        """
        Get path to BLAST database for species/locus.

        Args:
            species: Species name (e.g., "human", "mouse")
            locus: Locus name (e.g., "IGHV", "IGHD")

        Returns:
            Path to BLAST database (without .ndb extension)

        Raises:
            FileNotFoundError: If database not found
        """
        db_dir = Path.home() / ".sadie" / "germlines" / "databases" / species
        db_path = db_dir / locus

        if not (db_path.with_suffix(".ndb")).exists():
            raise FileNotFoundError(
                f"Germline database not found: {species}/{locus}\n"
                f"Run: sadie germlines populate -s {species}"
            )

        return db_path
```

**Usage in AIRR:**

```python
# In sadie/airr/airr.py

from sadie.germlines import GermlineLocator

class Airr:
    def __init__(self, species: str):
        self.species = species
        self.locator = GermlineLocator()

        # Get database paths
        self.ighv_db = self.locator.get_database_path(species, "IGHV")
        self.ighd_db = self.locator.get_database_path(species, "IGHD")
        self.ighj_db = self.locator.get_database_path(species, "IGHJ")
```

### Integration with Renumbering Module

Renumbering uses the same germline databases via `GermlineLocator`:

```python
# In sadie/renumbering/renumber.py

from sadie.germlines import GermlineLocator

class Renumber:
    def __init__(self, species: str):
        self.locator = GermlineLocator()
        self.germline_db = self.locator.get_database_path(species, "IGHV")
```

### Environment Variable Control

**Feature Flag: `SADIE_USE_GERMLINES_MODULE`**

```python
import os

def use_germlines_module() -> bool:
    """Check if germlines module should be used."""
    env_var = os.getenv("SADIE_USE_GERMLINES_MODULE", "true").lower()
    return env_var in ("true", "1", "yes")

# In AIRR module:
if use_germlines_module():
    locator = GermlineLocator()
    db_path = locator.get_database_path(species, locus)
else:
    # Fall back to G3 API (deprecated)
    db_path = g3_api_get_database(species, locus)
```

---

## Code Structure

### Directory Organization

```
sadie/
├── germlines/
│   ├── __init__.py
│   ├── cli.py                  # CLI commands (populate, status)
│   ├── providers/
│   │   ├── __init__.py
│   │   ├── base.py            # GermlineProvider ABC
│   │   ├── imgt.py            # IMGTProvider
│   │   ├── ogrdb.py           # OGRDBProvider
│   │   ├── vdjbase.py         # VDJbaseProvider
│   │   └── custom.py          # CustomProvider
│   ├── pipeline/
│   │   ├── __init__.py
│   │   ├── downloader.py      # GermlineDownloader
│   │   ├── merger.py          # GermlineMerger
│   │   ├── builder.py         # DatabaseBuilder
│   │   └── validator.py       # DatabaseValidator
│   ├── locator.py             # GermlineLocator
│   └── utils.py               # Utility functions
├── airr/
│   ├── __init__.py
│   └── airr.py                # AIRR annotation (uses germlines)
└── renumbering/
    ├── __init__.py
    └── renumber.py            # Renumbering (uses germlines)
```

### Module Dependencies

```
germlines (self-contained)
    ↓
airr (depends on germlines)
    ↓
renumbering (depends on germlines)
```

---

## Testing Strategy

### Unit Tests

```python
# tests/germlines/test_providers.py
def test_imgt_provider_download():
    """Test IMGT provider downloads correctly."""
    pass

# tests/germlines/test_merger.py
def test_priority_merging():
    """Test priority-based sequence merging."""
    pass

# tests/germlines/test_builder.py
def test_blast_database_build():
    """Test BLAST database building."""
    pass
```

### Integration Tests

```python
# tests/integration/test_germlines_airr.py
def test_airr_with_germlines():
    """Test AIRR annotation uses germlines module."""
    # 1. Populate test germlines
    # 2. Run AIRR annotation
    # 3. Verify results use local germlines
    pass
```

### End-to-End Tests

```bash
# tests/e2e/test_full_pipeline.sh

# Full pipeline test
sadie germlines populate -s human --force
sadie airr test_sequences.fasta -o results.tsv

# Verify results
if [ -f results.tsv ]; then
    echo "✓ Full pipeline test passed"
else
    echo "✗ Full pipeline test failed"
    exit 1
fi
```

---

## Performance Considerations

### Download Optimization

- **Parallel downloads**: Each species downloads independently (potential for parallelization)
- **Checkpoint resume**: Avoid re-downloading successful species
- **Compression**: FASTA files stored compressed when possible

### Build Optimization

- **Incremental builds**: Only rebuild changed species
- **Cached databases**: Reuse existing databases if sources unchanged
- **Parallel builds**: Build multiple loci concurrently

### Runtime Optimization

- **Local access**: No network latency (main benefit)
- **BLAST indexing**: Pre-built BLAST databases for fast queries
- **Memory efficiency**: Stream large FASTA files rather than loading entirely

---

## Future Extensibility

### Adding New Providers

```python
# Example: Add AlphaFold Germline Provider

from sadie.germlines.providers.base import GermlineProvider

class AlphaFoldProvider(GermlineProvider):
    @property
    def name(self) -> str:
        return "alphafold"

    def get_available_species(self) -> List[str]:
        return ["human", "mouse"]

    def download(self, species: str, output_dir: Path) -> bool:
        # Implementation
        pass

    # ... implement other abstract methods

# Register provider
PROVIDERS['alphafold'] = AlphaFoldProvider()
```

### Custom Build Steps

```python
# Add custom builder to pipeline

class CustomBuilder:
    def build(self, merged_dir: Path, output_dir: Path):
        # Custom build logic
        pass

# In pipeline orchestration:
builders = [
    BlastDatabaseBuilder(),
    AuxiliaryFileBuilder(),
    InternalDataBuilder(),
    CustomBuilder(),  # Add custom builder
]
```

---

## See Also

- [CLI Reference](cli-reference.md) - Command-line interface
- [Provider Guide](provider-guide.md) - Provider details
- [Custom Sequences](custom-sequences.md) - Adding custom data
- [Troubleshooting](troubleshooting.md) - Common issues
