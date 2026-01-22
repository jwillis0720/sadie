import pytest
from pathlib import Path


class TestHMMBuilder:
    @pytest.fixture
    def hmm_env(self, tmp_path):
        gapped_dir = tmp_path / "gapped"
        gapped_dir.mkdir()

        v_fasta = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG...TCTGGGGCT
>IGHV1-2*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG...TCTGGGGCT
>IGHV3-23*01|Homo sapiens|F
GAGGTGCAG...CTGGTGGAG...TCTGGGGGC
>IGHV4-34*01|Homo sapiens|F
CAGGTGCAG...CTGCAACAG...TCTGGGTCT
"""
        (gapped_dir / "IGHV.fasta").write_text(v_fasta)

        j_fasta = """>IGHJ4*01|Homo sapiens|F
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
>IGHJ6*01|Homo sapiens|F
ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG
>IGHJ5*01|Homo sapiens|F
ACAACTGGTTCGACTCCTGGGGCCAAGGCACCCTGGTCACCGTCTCCTCAG
"""
        (gapped_dir / "IGHJ.fasta").write_text(j_fasta)

        return tmp_path

    def test_build_creates_stockholm_file(self, hmm_env):
        from sadie.germlines.builders.hmm import HMMBuilder

        builder = HMMBuilder()
        output_dir = hmm_env / "output"

        builder.build_for_species(
            "human",
            source_dir=hmm_env / "gapped",
            output_dir=output_dir
        )

        sto_file = output_dir / "human_HV.sto"
        assert sto_file.exists()

    def test_stockholm_format_valid(self, hmm_env):
        from sadie.germlines.builders.hmm import HMMBuilder

        builder = HMMBuilder()
        output_dir = hmm_env / "output"

        builder.build_for_species(
            "human",
            source_dir=hmm_env / "gapped",
            output_dir=output_dir
        )

        sto_file = output_dir / "human_HV.sto"
        content = sto_file.read_text()

        assert "# STOCKHOLM 1.0" in content
        assert "#=GF ID human_HV" in content
        assert "//" in content

    def test_minimum_sequences_required(self, tmp_path):
        from sadie.germlines.builders.hmm import HMMBuilder

        gapped_dir = tmp_path / "gapped"
        gapped_dir.mkdir()

        v_fasta = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAG
>IGHV1-2*01|Homo sapiens|F
CAGGTGCAG
"""
        (gapped_dir / "IGHV.fasta").write_text(v_fasta)

        builder = HMMBuilder()
        output_dir = tmp_path / "output"

        builder.build_for_species(
            "human",
            source_dir=gapped_dir,
            output_dir=output_dir
        )

        sto_file = output_dir / "human_HV.sto"
        assert not sto_file.exists()

    def test_max_sequences_truncated(self, tmp_path):
        from sadie.germlines.builders.hmm import HMMBuilder

        gapped_dir = tmp_path / "gapped"
        gapped_dir.mkdir()

        lines = []
        for i in range(150):
            lines.append(f">IGHV{i}*01|Homo sapiens|F")
            lines.append("CAGGTGCAG")

        (gapped_dir / "IGHV.fasta").write_text("\n".join(lines))

        builder = HMMBuilder()
        output_dir = tmp_path / "output"

        builder.build_for_species(
            "human",
            source_dir=gapped_dir,
            output_dir=output_dir
        )

        sto_file = output_dir / "human_HV.sto"
        content = sto_file.read_text()

        seq_count = content.count("IGHV")
        assert seq_count == 100


class TestGetGappedSequences:
    def test_returns_tuples(self, tmp_path):
        from sadie.germlines.builders.hmm import get_gapped_sequences
        from sadie.germlines.manager import GermlineManager
        from sadie.germlines.providers.imgt import IMGTProvider

        (tmp_path / "imgt" / "human").mkdir(parents=True)

        v_fasta = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG
"""
        v_gapped = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG
"""
        (tmp_path / "imgt" / "human" / "IGHV.fasta").write_text(v_fasta)
        (tmp_path / "imgt" / "human" / "IGHV_gapped.fasta").write_text(v_gapped)

        provider = IMGTProvider(data_dir=tmp_path / "imgt")
        manager = GermlineManager(providers=[provider])

        result = get_gapped_sequences(manager, "human", "V")

        assert isinstance(result, list)
        if result:
            assert isinstance(result[0], tuple)
            assert len(result[0]) == 2
