# Integration API Contracts

**Feature**: 002-germline-integration
**Date**: 2026-01-20

## Overview

This document defines the programmatic interfaces for integrating the germlines module with existing SADIE components. These are Python API contracts, not REST/GraphQL endpoints.

## Public Integration APIs

### 1. Feature Flag API

**Module**: `sadie.germlines.utils.feature_flags`

```python
def use_germlines_module() -> bool:
    """
    Determine whether to use germlines module or fall back to G3 API.

    Environment Variable:
        SADIE_USE_GERMLINES_MODULE: "true" (default) or "false"

    Returns:
        bool: True to use germlines module, False to use G3 API

    Examples:
        >>> use_germlines_module()
        True

        >>> os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
        >>> use_germlines_module()
        False
    """
```

**Contract**:
- **Input**: None (reads environment)
- **Output**: Boolean
- **Side Effects**: Logs deprecation warning if false
- **Thread Safety**: Yes (read-only environment variable)
- **Performance**: O(1), single environment variable lookup

### 2. G3 Adapter API

**Module**: `sadie.germlines.g3_adapter`

```python
class GermlineToG3Adapter:
    def to_g3_format(self, gene: GermlineGene) -> Dict[str, Any]:
        """
        Convert GermlineGene to G3 API response format.

        Args:
            gene: GermlineGene object from germlines module

        Returns:
            Dictionary matching G3 API JSON structure

        Raises:
            ValueError: If gene has invalid/missing required fields

        Examples:
            >>> adapter = GermlineToG3Adapter()
            >>> gene = get_gene_by_name("human", "IGHV1-69*01")
            >>> g3_dict = adapter.to_g3_format(gene)
            >>> g3_dict["gene"]
            'IGHV1-69*01'
            >>> g3_dict["imgt"]["imgt_functional"]
            'F'
        """

    def to_g3_format_batch(
        self,
        genes: List[GermlineGene]
    ) -> List[Dict[str, Any]]:
        """
        Convert multiple GermlineGene objects to G3 format.

        Args:
            genes: List of GermlineGene objects

        Returns:
            List of G3-formatted dictionaries

        Performance:
            O(n) where n is number of genes
        """
```

**Contract**:
- **Input**: GermlineGene object(s) (validated Pydantic model)
- **Output**: Dictionary matching G3 API structure
- **Guarantees**:
  - All required G3 fields present
  - IMGT regions included if available in gene
  - Latin names mapped correctly for common species
- **Error Handling**: ValueError on invalid input
- **Thread Safety**: Yes (stateless transformation)
- **Performance**: O(1) per gene, no I/O

### 3. Local HMM Builder API

**Module**: `sadie.germlines.renumbering_integration`

```python
class LocalHMMBuilder:
    def get_hmm(
        self,
        species: str,
        chain: str,
        source: str = "imgt"
    ) -> pyhmmer.plan7.HMM:
        """
        Get or build HMM model for species/chain combination.

        Args:
            species: Species name (e.g., "human", "mouse")
            chain: Chain type ("H", "K", or "L")
            source: Data source (default: "imgt")

        Returns:
            Compiled HMM model from pyhmmer

        Raises:
            ValueError: If no sequences found for species/chain
            FileNotFoundError: If germlines data not populated

        Caching:
            HMMs cached in germlines/hmms/{species}_{chain}.hmm
            Cached HMMs loaded on subsequent calls

        Performance:
            - Cached: O(1), ~10ms (file load)
            - Build: O(n) where n=sequences, ~5-10s (Stockholm + pyhmmer)

        Examples:
            >>> builder = LocalHMMBuilder()
            >>> hmm = builder.get_hmm("human", "H")
            >>> hmm.name
            b'human_H'
        """
```

**Contract**:
- **Input**: Species (str), Chain (str), Source (str)
- **Output**: pyhmmer HMM object
- **Side Effects**: Writes cache file on first build
- **Error Handling**: Raises ValueError/FileNotFoundError with clear messages
- **Thread Safety**: Yes (LRU cache thread-safe, file writes atomic)
- **Performance**: Fast on cache hit, slow on cache miss

### 4. IgBLAST Path Integration API

**Module**: `sadie.airr.igblast.germline`

```python
class GermlineData:
    def __init__(
        self,
        name: str,
        receptor: str = "Ig",
        database_dir: Optional[Union[str, Path]] = None,
        scheme: str = "imgt"
    ):
        """
        Initialize germline database paths for IgBLAST.

        Args:
            name: Species name (e.g., "human")
            receptor: Receptor type (default: "Ig")
            database_dir: Override database directory (optional)
            scheme: Numbering scheme (default: "imgt")

        Behavior:
            - If database_dir provided: Use custom path
            - Elif SADIE_USE_GERMLINES_MODULE=true: Use germlines/igblast/
            - Else: Use airr/data/germlines/ (legacy G3 paths)

        Attributes:
            base_dir: Base directory for databases
            v_gene_dir: Prefix for V gene BLAST database
            d_gene_dir: Prefix for D gene BLAST database
            j_gene_dir: Prefix for J gene BLAST database
            c_gene_dir: Prefix for C gene BLAST database
            aux_path: Auxiliary file path for J genes
            igdata: Internal data directory

        Raises:
            FileNotFoundError: If database files not found at determined path

        Examples:
            >>> gd = GermlineData("human")
            >>> gd.v_gene_dir
            PosixPath('.../germlines/igblast/database/human/human_V')
        """
```

**Contract**:
- **Input**: Species name, optional overrides
- **Output**: Initialized object with validated paths
- **Side Effects**: None (read-only path validation)
- **Error Handling**: FileNotFoundError if databases missing
- **Thread Safety**: Yes (immutable after construction)
- **Performance**: O(1), file existence checks only

## Integration Contracts for Existing APIs

### 5. Reference System Integration

**Module**: `sadie.reference.reference`

**Proposed Addition**:
```python
class Reference:
    def __init__(
        self,
        endpoint: str = _endpoint,
        use_germlines: bool = False  # NEW PARAMETER
    ):
        """
        Initialize reference object.

        Args:
            endpoint: G3 API endpoint (ignored if use_germlines=True)
            use_germlines: Use local germlines module instead of G3 API

        Backwards Compatibility:
            - Default use_germlines=False maintains existing behavior
            - Existing code without parameter continues using G3
        """

    def _get_gene(self, gene: GeneEntry) -> Dict[str, str]:
        """
        Get single gene (from germlines or G3).

        Returns:
            Dictionary in G3 API format (same structure regardless of source)
        """
```

**Contract**:
- **Input**: GeneEntry model (species, gene name, source)
- **Output**: G3-format dictionary
- **Guarantees**:
  - Output format identical between germlines and G3 backends
  - All required fields present
  - IMGT regions included for V genes
- **Error Handling**: G3Error if gene not found (consistent with existing)
- **Backwards Compatibility**: ✅ Default behavior unchanged
- **Performance**: Local query ~10x faster than G3 API

### 6. HMMER Integration

**Module**: `sadie.renumbering.aligners.hmmer`

**Existing API** (Enhanced):
```python
class HMMER:
    def get_hmm_models(
        self,
        species: Optional[Union[List[Species], Species]] = None,
        chains: Optional[Union[List[Chain], Chain]] = None,
        source: Source = "imgt",
        use_numbering_hmms: bool = False
    ) -> List[pyhmmer.plan7.HMM]:
        """
        Get HMM models for species/chains.

        Priority Order (NEW):
            1. LocalHMMBuilder (if SADIE_USE_GERMLINES_MODULE=true)
            2. Numbering legacy HMMs (fallback)
            3. G3 API HMMs (legacy fallback)

        Returns:
            List of pyhmmer HMM objects

        Backwards Compatibility:
            - Existing behavior preserved when feature flag disabled
            - Graceful fallback on LocalHMMBuilder failure
        """
```

**Contract**:
- **Input**: Species list, chain list, source, flags
- **Output**: List of HMM objects
- **Guarantees**:
  - Returns valid HMMs for all requested species/chains (or skips unavailable)
  - HMM format identical across backends
  - Performance equivalent or better than G3
- **Error Handling**: Logs warnings, falls back gracefully
- **Backwards Compatibility**: ✅ No breaking changes
- **Performance**: Cached HMMs ~100x faster than G3 API builds

## Contract Validation

### Backwards Compatibility Tests

```python
# Test 1: Default behavior unchanged
def test_backwards_compatibility_igblast():
    # Without feature flag
    os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
    gd = GermlineData("human")
    # Should use legacy paths
    assert "airr/data/germlines" in str(gd.base_dir)

# Test 2: Reference format consistency
def test_reference_format_consistency():
    # Compare G3 vs germlines output
    ref_g3 = Reference(use_germlines=False)
    ref_local = Reference(use_germlines=True)

    gene_g3 = ref_g3._get_gene(GeneEntry(...))
    gene_local = ref_local._get_gene(GeneEntry(...))

    # Same keys, values may differ slightly (versions)
    assert gene_g3.keys() == gene_local.keys()

# Test 3: HMM output equivalence
def test_hmm_output_equivalence():
    # G3 vs LocalHMMBuilder
    hmmer_g3 = HMMER(species="human")  # Uses G3
    hmmer_local = HMMER(species="human")  # Uses LocalHMMBuilder

    # HMMs should have same structure
    assert hmmer_g3.hmms[0].name == hmmer_local.hmms[0].name
```

### Performance Contracts

| Operation | G3 Backend | Germlines Backend | Target |
|-----------|------------|-------------------|--------|
| Gene lookup | 100-500ms | 10-50ms | <200ms |
| HMM build (uncached) | 5-10s | 5-10s | <15s |
| HMM load (cached) | N/A (no cache) | 10-20ms | <100ms |
| IgBLAST execution | Same | Same | <5s per 1000 seqs |

### Error Contract Matrix

| Scenario | Expected Behavior | Error Type |
|----------|-------------------|------------|
| Germlines not populated | Clear error message with setup instructions | FileNotFoundError |
| Species not available | List available species in error | ValueError |
| Network unavailable (G3 mode) | Connection error, suggest germlines mode | G3Error |
| Invalid gene name | "Gene not found" message | G3Error / ValueError |
| Corrupted HMM cache | Rebuild HMM automatically | Logged warning |

## Summary

All integration APIs maintain backwards compatibility while providing new germlines-based functionality. Key contracts:

1. **Feature flags** control backend selection transparently
2. **G3 format adapter** ensures format consistency
3. **Local HMM builder** matches G3 performance with caching benefits
4. **Path integration** seamlessly switches database locations
5. **Error handling** provides clear, actionable messages
6. **Performance** meets or exceeds G3 API benchmarks

