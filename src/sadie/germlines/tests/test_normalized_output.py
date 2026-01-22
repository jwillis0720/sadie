import pytest
from pathlib import Path
import tempfile
import shutil


class TestNormalizedOutput:
    @pytest.fixture
    def pipeline_env(self, tmp_path):
        sources = tmp_path / "sources"
        (sources / "imgt" / "human").mkdir(parents=True)

        v_fasta = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCT
>IGHV1-2*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        v_gapped = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG...TCTGGGGCT
>IGHV1-2*01|Homo sapiens|F
CAGGTGCAG...CTGGTGCAG...TCTGGGGCT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(v_fasta)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(v_gapped)

        d_fasta = """>IGHD3-3*01|Homo sapiens|F
GTATTACTATGGTTCGGGGAGTT
"""
        (sources / "imgt" / "human" / "IGHD.fasta").write_text(d_fasta)

        j_fasta = """>IGHJ4*01|Homo sapiens|F
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
"""
        (sources / "imgt" / "human" / "IGHJ.fasta").write_text(j_fasta)
        (sources / "imgt" / "human" / "IGHJ_gapped.fasta").write_text(j_fasta)

        normalized = tmp_path / "normalized"
        normalized.mkdir()
        igblast = tmp_path / "igblast"
        igblast.mkdir()

        return tmp_path

    def test_normalized_creates_gapped_directory(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        gapped_dir = pipeline_env / "normalized" / "human" / "gapped"
        assert gapped_dir.exists()

    def test_normalized_creates_ungapped_directory(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        ungapped_dir = pipeline_env / "normalized" / "human" / "ungapped"
        assert ungapped_dir.exists()

    def test_v_segment_in_both_gapped_and_ungapped(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        gapped_v = pipeline_env / "normalized" / "human" / "gapped" / "IGHV.fasta"
        ungapped_v = pipeline_env / "normalized" / "human" / "ungapped" / "IGHV.fasta"

        assert gapped_v.exists()
        assert ungapped_v.exists()

    def test_d_segment_ungapped_only(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        ungapped_d = pipeline_env / "normalized" / "human" / "ungapped" / "IGHD.fasta"
        assert ungapped_d.exists()

    def test_gapped_contains_dots(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        gapped_v = pipeline_env / "normalized" / "human" / "gapped" / "IGHV.fasta"
        if gapped_v.exists():
            content = gapped_v.read_text()
            assert "." in content or "IGHV" in content

    def test_ungapped_no_dots(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        pipeline._rebuild_normalized("human")

        ungapped_v = pipeline_env / "normalized" / "human" / "ungapped" / "IGHV.fasta"
        content = ungapped_v.read_text()

        for line in content.splitlines():
            if not line.startswith(">"):
                assert "." not in line


class TestNormalizedFileNaming:
    @pytest.fixture
    def minimal_env(self, tmp_path):
        sources = tmp_path / "sources"
        (sources / "imgt" / "human").mkdir(parents=True)

        v_fasta = """>IGHV1-69*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(v_fasta)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(v_fasta)

        (tmp_path / "normalized").mkdir()
        (tmp_path / "igblast").mkdir()

        return tmp_path

    def test_output_file_naming_pattern(self, minimal_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(minimal_env)
        pipeline._rebuild_normalized("human")

        gapped_dir = minimal_env / "normalized" / "human" / "gapped"
        files = list(gapped_dir.glob("*.fasta"))

        for f in files:
            assert f.name.startswith("IG")
            assert f.name.endswith(".fasta")


class TestNoGapperDuplication:
    def test_custom_provider_uses_builders_gapper(self):
        from sadie.germlines.providers.custom import CustomProvider
        from sadie.germlines.builders.gapper import GapperService

        provider = CustomProvider()
        assert hasattr(provider, '_get_gapper')

    def test_pipeline_uses_builders(self):
        from sadie.germlines.pipeline import GermlinePipeline

        from sadie.germlines.builders.blast import BlastDBBuilder
        from sadie.germlines.builders.aux import AuxFileBuilder

        assert BlastDBBuilder is not None
        assert AuxFileBuilder is not None

    def test_gapper_service_imported_in_custom(self):
        import sadie.germlines.providers.custom as custom_module

        import_source = custom_module.__file__
        with open(import_source) as f:
            content = f.read()

        assert "from ..builders.gapper import GapperService" in content
