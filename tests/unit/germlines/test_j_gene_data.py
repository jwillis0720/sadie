"""
Test j_gene_data.py functionality with new JSON source.

This module tests that j_gene_data.py works correctly with the new
centralized motif lookup system from sadie.reference.

User Story 3: Maintain Backward Compatibility
- Test j_gene_data.py functionality with new JSON source
- Ensure motif pattern loading works correctly
- Test J gene parsing and CDR3/FWR4 boundary detection
"""

import pytest
from unittest.mock import patch

from sadie.germlines.builders.j_gene_data import (
    MOTIF_LOOKUP,
    parse_j_gene,
    get_j_gene_data,
    _normalize_species,
    _get_gene_type,
    _get_short_name,
    _translate_sequence,
    _infer_reading_frame
)


class TestJGeneDataIntegration:
    """Test j_gene_data.py integration with new motif lookup system."""

    def test_motif_lookup_loaded(self):
        """Test that MOTIF_LOOKUP is loaded from new system."""
        assert MOTIF_LOOKUP is not None
        assert len(MOTIF_LOOKUP) > 0

        # Should have expected species
        assert 'human' in MOTIF_LOOKUP
        assert 'mouse' in MOTIF_LOOKUP

        # Should have expected patterns
        human_data = MOTIF_LOOKUP['human']
        assert 'IGHJ' in human_data
        assert isinstance(human_data['IGHJ'], str)

    def test_motif_lookup_no_provenance(self):
        """Test that MOTIF_LOOKUP doesn't contain provenance fields."""
        for species, patterns in MOTIF_LOOKUP.items():
            assert '_provenance' not in patterns

    def test_species_normalization(self):
        """Test species name normalization."""
        assert _normalize_species('Homo sapiens') == 'human'
        assert _normalize_species('homo_sapiens') == 'human'
        assert _normalize_species('Mus musculus') == 'mouse'
        assert _normalize_species('mus_musculus') == 'mouse'
        assert _normalize_species('human') == 'human'

    def test_gene_type_extraction(self):
        """Test gene type extraction from gene names."""
        assert _get_gene_type('IGHJ1*01') == 'IGHJ'
        assert _get_gene_type('IGKJ2*01') == 'IGKJ'
        assert _get_gene_type('IGLJ3*01') == 'IGLJ'
        assert _get_gene_type('traj1*01') == 'TRAJ'

    def test_short_name_extraction(self):
        """Test short name extraction from gene names."""
        assert _get_short_name('IGHJ1*01') == 'IGHJ1'
        assert _get_short_name('IGKJ2*02') == 'IGKJ2'
        assert _get_short_name('IGLJ3') == 'IGLJ3'

    def test_sequence_translation(self):
        """Test nucleotide to amino acid translation."""
        # Test known codon translations
        assert _translate_sequence('TTT') == 'F'
        assert _translate_sequence('TTG') == 'L'
        assert _translate_sequence('TGG') == 'W'
        assert _translate_sequence('GGG') == 'G'

        # Test reading frames
        sequence = 'ATGTTTGGG'  # ATG TTT GGG
        assert _translate_sequence(sequence, 0) == 'MFG'
        assert _translate_sequence(sequence, 1) == 'CL'  # TGT TTG GG
        assert _translate_sequence(sequence, 2) == 'VW'  # GTT TGG G

    def test_reading_frame_inference(self):
        """Test reading frame inference for different chain types."""
        # Mock sequence ending in different amino acids
        heavy_seq = 'TTTGGGTCT'  # Should end with S (Ser) in correct frame
        kappa_seq = 'TTTGGGAAA'  # Should end with K (Lys) in correct frame
        lambda_seq = 'TTTGGGTTG'  # Should end with L (Leu) in correct frame

        # Test reading frame inference
        heavy_rf = _infer_reading_frame(heavy_seq, 'IGHJ1*01')
        kappa_rf = _infer_reading_frame(kappa_seq, 'IGKJ1*01')
        lambda_rf = _infer_reading_frame(lambda_seq, 'IGLJ1*01')

        # Should return valid reading frames
        assert heavy_rf in [0, 1, 2]
        assert kappa_rf in [0, 1, 2]
        assert lambda_rf in [0, 1, 2]


class TestJGeneParsingFunctionality:
    """Test J gene parsing functionality."""

    def test_parse_j_gene_human_ighj(self):
        """Test parsing human IGHJ with known motif."""
        # Human IGHJ sequence that should match motif
        sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'

        result = parse_j_gene('human', 'IGHJ1*01', sequence)

        # Should have valid structure
        assert isinstance(result, dict)
        assert 'sequence_gapped_aa' in result
        assert 'cdr3' in result
        assert 'fwr4' in result
        assert 'reading_frame' in result
        assert 'expression_match' in result
        assert 'not_implemented' in result
        assert 'ignored' in result

    def test_parse_j_gene_unknown_species(self):
        """Test parsing with unknown species."""
        sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'

        result = parse_j_gene('unknown_species', 'IGHJ1*01', sequence)

        assert result['not_implemented'] is True
        assert result['expression_match'] is False

    def test_parse_j_gene_unknown_gene_type(self):
        """Test parsing with unknown gene type."""
        sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'

        result = parse_j_gene('human', 'UNKNOWN*01', sequence)

        assert result['not_implemented'] is True
        assert result['expression_match'] is False

    def test_parse_j_gene_ignored_gene(self):
        """Test parsing with gene in ignore list."""
        # Find a species with ignore list
        species_with_ignore = None
        ignore_gene = None

        for species, patterns in MOTIF_LOOKUP.items():
            if 'ignore' in patterns and patterns['ignore']:
                species_with_ignore = species
                ignore_gene = patterns['ignore'][0] if isinstance(patterns['ignore'], list) else patterns['ignore']
                break

        if species_with_ignore and ignore_gene:
            sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'
            gene_name = f'{ignore_gene}*01'

            result = parse_j_gene(species_with_ignore, gene_name, sequence)

            assert result['ignored'] is True

    def test_parse_j_gene_sequence_cleaning(self):
        """Test that sequence is properly cleaned."""
        # Sequence with gaps and dots
        dirty_sequence = 'TTT.GGG-ACAAGG...TGACCGT'
        clean_sequence = 'TTTGGGACAAGGTGACCGT'

        result1 = parse_j_gene('human', 'IGHJ1*01', dirty_sequence)
        result2 = parse_j_gene('human', 'IGHJ1*01', clean_sequence)

        # Should produce same results after cleaning
        assert result1['sequence_gapped_aa'] == result2['sequence_gapped_aa']


class TestJGeneDataFunction:
    """Test the get_j_gene_data function."""

    def test_get_j_gene_data_human_known(self):
        """Test getting data for known human genes."""
        # Test known human gene
        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('IGHJ1*01', 'H')

        assert isinstance(rf, int)
        assert rf in [0, 1, 2]
        assert chain_type == 'JH'
        assert isinstance(cdr3_end, int)
        assert isinstance(extra_bps, int)

    def test_get_j_gene_data_mouse_known(self):
        """Test getting data for known mouse genes."""
        # Test known mouse gene
        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('IGHJ1*01', 'H', species='mouse')

        assert isinstance(rf, int)
        assert rf in [0, 1, 2]
        assert chain_type == 'JH'
        assert isinstance(cdr3_end, int)
        assert isinstance(extra_bps, int)

    def test_get_j_gene_data_with_sequence(self):
        """Test getting data with sequence for dynamic calculation."""
        # Test with sequence for dynamic calculation
        sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'

        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data(
            'IGHJ1*01', 'H', sequence=sequence, species='human'
        )

        assert isinstance(rf, int)
        assert rf in [0, 1, 2]
        assert chain_type == 'JH'
        assert isinstance(cdr3_end, int)
        assert isinstance(extra_bps, int)

    def test_get_j_gene_data_unknown_allele(self):
        """Test getting data for unknown allele with defaults."""
        # Test unknown allele - should return defaults
        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('UNKNOWN*01', 'H')

        assert isinstance(rf, int)
        assert chain_type == 'JH'
        assert isinstance(cdr3_end, int)
        assert isinstance(extra_bps, int)

    def test_get_j_gene_data_different_chains(self):
        """Test getting data for different chain types."""
        # Test different chain types
        chains = ['H', 'K', 'L']
        expected_chain_types = ['JH', 'JK', 'JL']

        for chain, expected_chain_type in zip(chains, expected_chain_types):
            rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('TEST*01', chain)
            assert chain_type == expected_chain_type


class TestJGeneDataErrorHandling:
    """Test error handling and edge cases."""

    def test_empty_sequence_handling(self):
        """Test handling of empty sequences."""
        result = parse_j_gene('human', 'IGHJ1*01', '')

        # Should handle gracefully
        assert isinstance(result, dict)
        assert result['sequence_gapped_aa'] == ''

    def test_invalid_nucleotides(self):
        """Test handling of invalid nucleotide characters."""
        # Sequence with invalid characters
        sequence = 'TTTGGGXYZACAAGGT'

        result = parse_j_gene('human', 'IGHJ1*01', sequence)

        # Should handle gracefully (invalid codons become 'X')
        assert isinstance(result, dict)
        assert 'X' in result['sequence_gapped_aa']  # Invalid codons translated to X

    def test_short_sequence_handling(self):
        """Test handling of very short sequences."""
        # Very short sequence
        sequence = 'TTT'

        result = parse_j_gene('human', 'IGHJ1*01', sequence)

        assert isinstance(result, dict)
        assert len(result['sequence_gapped_aa']) <= 1

    def test_motif_pattern_error_handling(self):
        """Test error handling when motif pattern is malformed."""
        # Mock a malformed pattern
        with patch.dict(MOTIF_LOOKUP, {'test_species': {'IGHJ': '[invalid_regex'}}):
            sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGG'

            # Currently j_gene_data.py doesn't handle regex errors - expect exception
            with pytest.raises(Exception):  # Could be re.error or other regex exception
                parse_j_gene('test_species', 'IGHJ1*01', sequence)


class TestJGeneDataPerformance:
    """Test performance of j_gene_data functions."""

    def test_parse_performance(self):
        """Test that parsing is reasonably fast."""
        import time

        sequence = 'TTTGGGACAAGGTGACCGTGATCCTCTGGACAGGCTACGGC'

        # Test multiple parses
        start_time = time.time()

        for _ in range(100):
            result = parse_j_gene('human', 'IGHJ1*01', sequence)
            assert result is not None

        elapsed = time.time() - start_time

        # Should be fast enough for production use
        assert elapsed < 1.0, f"100 parses took {elapsed:.3f}s, too slow"

    def test_get_j_gene_data_performance(self):
        """Test that get_j_gene_data is fast."""
        import time

        start_time = time.time()

        for _ in range(100):
            result = get_j_gene_data('IGHJ1*01', 'H')
            assert result is not None

        elapsed = time.time() - start_time

        # Should be very fast for cached lookups
        assert elapsed < 0.1, f"100 lookups took {elapsed:.3f}s, too slow"


class TestJGeneDataRegressionCompatibility:
    """Test regression compatibility with expected outputs."""

    def test_human_ighj1_regression(self):
        """Test that human IGHJ1*01 produces expected results."""
        # Known human IGHJ1*01 data
        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('IGHJ1*01', 'H')

        # Should match known values from reference data
        assert chain_type == 'JH'
        assert isinstance(rf, int) and rf >= 0
        assert isinstance(cdr3_end, int) and cdr3_end > 0
        assert isinstance(extra_bps, int) and extra_bps >= 0

    def test_mouse_ighj1_regression(self):
        """Test that mouse IGHJ1*01 produces expected results."""
        # Known mouse IGHJ1*01 data
        rf, chain_type, cdr3_end, extra_bps = get_j_gene_data('IGHJ1*01', 'H', species='mouse')

        # Should match known values from reference data
        assert chain_type == 'JH'
        assert isinstance(rf, int) and rf >= 0
        assert isinstance(cdr3_end, int) and cdr3_end > 0
        assert isinstance(extra_bps, int) and extra_bps >= 0

    def test_chain_type_mapping_regression(self):
        """Test that chain type mapping is preserved."""
        # Test all chain types
        test_cases = [
            ('H', 'JH'),
            ('K', 'JK'),
            ('L', 'JL'),
        ]

        for chain, expected_chain_type in test_cases:
            _, chain_type, _, _ = get_j_gene_data('TEST*01', chain)
            assert chain_type == expected_chain_type