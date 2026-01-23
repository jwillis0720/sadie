import pytest
from pathlib import Path


class TestCodonAwareGapping:
    @pytest.fixture
    def gapper_env(self, tmp_path):
        template_dir = tmp_path / "templates" / "human"
        template_dir.mkdir(parents=True)

        v_gapped = """>acc|IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG...TCTGGGGCT
"""
        (template_dir / "IGHV_gapped.fasta").write_text(v_gapped)

        return tmp_path / "templates"

    def test_translate_produces_amino_acids(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        aa = gapper._translate("CAGGTGCAGCTG")

        assert aa is not None
        assert len(aa) == 4
        assert aa == "QVQL"

    def test_translate_truncates_to_codon_boundary(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        aa = gapper._translate("CAGGTGCAGCT")

        assert aa is not None
        assert len(aa) == 3
        assert aa == "QVQ"

    def test_translate_handles_gaps(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        aa = gapper._translate("CAG...GTG...CAG")

        assert aa is not None
        assert "Q" in aa

    def test_gap_positions_extracted_correctly(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        positions = gapper._extract_gap_positions("CAG...GTG")

        assert len(positions) == 3
        assert positions == [3, 4, 5]

    def test_d_segment_not_gapped(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        result = gapper.gap_sequence(
            sequence="GTATTACTATGGTTCGGGGAGTT",
            segment="D",
            chain="H",
            gene_name="IGHD3-3*01"
        )

        assert result is None

    def test_unknown_segment_returns_none(self, gapper_env):
        from sadie.germlines.builders.gapper import GapperService

        gapper = GapperService(template_dir=gapper_env)
        result = gapper.gap_sequence(
            sequence="CAGGTGCAG",
            segment="X",
            chain="H",
            gene_name="TEST"
        )

        assert result is None

    def test_codon_table_complete(self, gapper_env):
        from sadie.germlines.builders.gapper import CODON_TABLE

        assert len(CODON_TABLE) == 64
        assert CODON_TABLE["ATG"] == "M"
        assert CODON_TABLE["TAA"] == "*"


class TestBatchGapping:
    @pytest.fixture
    def batch_env(self, tmp_path):
        template_dir = tmp_path / "templates" / "human"
        template_dir.mkdir(parents=True)

        v_gapped = """>acc|IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG
"""
        (template_dir / "IGHV_gapped.fasta").write_text(v_gapped)

        return tmp_path / "templates"

    def test_batch_processes_multiple(self, batch_env):
        from sadie.germlines.builders.gapper import gap_sequences_batch

        sequences = [
            ("IGHV1-69*01", "CAGGTGCAGCTGGTGCAG", "V", "H"),
            ("IGHD3-3*01", "GTATTACTAT", "D", "H"),
        ]

        results = gap_sequences_batch(sequences, batch_env)

        assert "IGHV1-69*01" in results
        assert "IGHD3-3*01" in results
        assert results["IGHD3-3*01"] is None
