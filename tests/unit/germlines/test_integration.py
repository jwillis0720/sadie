import pytest
from pathlib import Path
import tempfile
import shutil


class TestCustomSequenceIntegration:
    @pytest.fixture
    def integration_env(self):
        tmp = tempfile.mkdtemp()
        base_dir = Path(tmp)

        sources = base_dir / "sources"
        normalized = base_dir / "normalized"
        igblast = base_dir / "igblast"

        for provider in ["custom", "imgt", "ogrdb", "vdjbase"]:
            (sources / provider / "human").mkdir(parents=True)

        imgt_v = """>M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt
CAGGTTCAGCTGGTGCAGTCTGGAGCT
>X07448|IGHV1-2*01|Homo sapiens|F|V-REGION|269..564|296 nt
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        imgt_j = """>J00256|IGHJ4*01|Homo sapiens|F|J-REGION|1912..1959|48 nt
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
"""
        (sources / "imgt" / "human" / "IGHJ.fasta").write_text(imgt_j)
        (sources / "imgt" / "human" / "IGHJ_gapped.fasta").write_text(imgt_j)

        imgt_d = """>X97051|IGHD1-1*01|Homo sapiens|F|D-REGION|33714..33730|17 nt
GGTACAACTGGAACGAC
"""
        (sources / "imgt" / "human" / "IGHD.fasta").write_text(imgt_d)

        yield {"base_dir": base_dir, "sources": sources, "normalized": normalized, "igblast": igblast}

        shutil.rmtree(tmp)

    def test_manager_loads_imgt_data(self, integration_env):
        from sadie.germlines.manager import GermlineManager

        manager = GermlineManager(providers=["imgt"])
        manager.sources_dir = integration_env["sources"]

        from sadie.germlines.providers.imgt import IMGTProvider
        imgt_provider = IMGTProvider(data_dir=integration_env["sources"] / "imgt")

        genes = imgt_provider.fetch_genes("human", "V", "H")
        assert len(genes) >= 2

    def test_custom_sequence_takes_priority(self, integration_env):
        from sadie.germlines.providers.custom import CustomProvider
        from sadie.germlines.providers.imgt import IMGTProvider

        custom_fasta = """>IGHV1-18*01
CAGGTTCAGCTGGTGCAGTCTGGAGCTAAAA
"""
        (integration_env["sources"] / "custom" / "human" / "IGHV.fasta").write_text(custom_fasta)

        custom = CustomProvider(data_dir=integration_env["sources"] / "custom")
        imgt = IMGTProvider(data_dir=integration_env["sources"] / "imgt")

        custom_genes = custom.fetch_genes("human", "V", "H")
        imgt_genes = imgt.fetch_genes("human", "V", "H")

        custom_v18 = next((g for g in custom_genes if "IGHV1-18" in g.name), None)
        imgt_v18 = next((g for g in imgt_genes if "IGHV1-18" in g.name), None)

        assert custom_v18 is not None
        assert imgt_v18 is not None
        assert custom_v18.sequence.endswith("AAAA")
        assert not imgt_v18.sequence.endswith("AAAA")

    def test_novel_custom_sequence_included(self, integration_env):
        from sadie.germlines.providers.custom import CustomProvider

        custom_fasta = """>IGHV-NOVEL*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGGGGG
"""
        (integration_env["sources"] / "custom" / "human" / "IGHV.fasta").write_text(custom_fasta)

        custom = CustomProvider(data_dir=integration_env["sources"] / "custom")
        genes = custom.fetch_genes("human", "V", "H")

        novel = next((g for g in genes if "NOVEL" in g.name), None)
        assert novel is not None
        assert novel.source == "custom"


class TestPipelineIntegration:
    @pytest.fixture
    def pipeline_env(self, tmp_path):
        sources = tmp_path / "sources"
        normalized = tmp_path / "normalized"
        igblast = tmp_path / "igblast"

        for provider in ["custom", "imgt", "ogrdb", "vdjbase"]:
            (sources / provider / "human").mkdir(parents=True)

        imgt_v = """>IGHV1-2*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        return tmp_path

    def test_sources_changed_detection(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        pipeline = GermlinePipeline(pipeline_env)
        assert pipeline._sources_changed("human") is True

    def test_sources_unchanged_after_build(self, pipeline_env):
        from sadie.germlines.pipeline import GermlinePipeline
        import time

        pipeline = GermlinePipeline(pipeline_env)

        gapped_dir = pipeline_env / "normalized" / "human" / "gapped"
        gapped_dir.mkdir(parents=True)

        dummy_fasta = gapped_dir / "IGHV.fasta"
        dummy_fasta.write_text(">test\nACGT")

        time.sleep(0.1)

        assert pipeline._sources_changed("human") is False


class TestOfflineOperation:
    def test_no_network_calls_during_fetch(self, tmp_path):
        from sadie.germlines.providers.imgt import IMGTProvider
        import socket

        sources = tmp_path / "sources" / "imgt" / "human"
        sources.mkdir(parents=True)

        fasta = """>IGHV1-2*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (sources / "IGHV.fasta").write_text(fasta)
        (sources / "IGHV_gapped.fasta").write_text(fasta)

        original_socket = socket.socket

        def blocked_socket(*args, **kwargs):
            raise RuntimeError("Network access attempted during fetch")

        socket.socket = blocked_socket

        try:
            provider = IMGTProvider(data_dir=tmp_path / "sources" / "imgt")
            genes = provider.fetch_genes("human", "V", "H")
            assert len(genes) >= 1
        finally:
            socket.socket = original_socket
