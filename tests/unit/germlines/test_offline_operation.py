import pytest
import socket
from pathlib import Path
import tempfile
import shutil


class TestOfflineOperation:
    @pytest.fixture
    def offline_env(self, tmp_path):
        sources = tmp_path / "sources"
        for provider in ["custom", "imgt", "ogrdb", "vdjbase"]:
            (sources / provider / "human").mkdir(parents=True)

        imgt_v = """>IGHV1-2*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        imgt_j = """>IGHJ4*01|Homo sapiens|F
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
"""
        (sources / "imgt" / "human" / "IGHJ.fasta").write_text(imgt_j)
        (sources / "imgt" / "human" / "IGHJ_gapped.fasta").write_text(imgt_j)

        imgt_d = """>IGHD1-1*01|Homo sapiens|F
GGTACAACTGGAACGAC
"""
        (sources / "imgt" / "human" / "IGHD.fasta").write_text(imgt_d)

        return tmp_path

    def test_pipeline_offline_no_network(self, offline_env):
        from sadie.germlines.pipeline import GermlinePipeline

        original_socket = socket.socket

        def blocked_socket(*args, **kwargs):
            raise RuntimeError("Network access attempted in offline mode")

        socket.socket = blocked_socket

        try:
            pipeline = GermlinePipeline(offline_env)
            assert pipeline._sources_changed("human") is True
        finally:
            socket.socket = original_socket

    def test_provider_fetch_offline(self, offline_env):
        from sadie.germlines.providers.imgt import IMGTProvider

        original_socket = socket.socket

        def blocked_socket(*args, **kwargs):
            raise RuntimeError("Network access attempted in offline mode")

        socket.socket = blocked_socket

        try:
            provider = IMGTProvider(data_dir=offline_env / "sources" / "imgt")
            genes = provider.fetch_genes("human", "V", "H")
            assert len(genes) >= 1
        finally:
            socket.socket = original_socket

    def test_manager_offline(self, offline_env):
        from sadie.germlines.providers.imgt import IMGTProvider

        original_socket = socket.socket

        def blocked_socket(*args, **kwargs):
            raise RuntimeError("Network access attempted in offline mode")

        socket.socket = blocked_socket

        try:
            provider = IMGTProvider(data_dir=offline_env / "sources" / "imgt")
            assert provider.is_available("human") is True
            genes = provider.fetch_genes("human", "V", "H")
            assert len(genes) >= 1
        finally:
            socket.socket = original_socket


class TestMissingDataErrors:
    def test_missing_species_error_message(self, tmp_path):
        from sadie.germlines.pipeline import GermlinePipeline, _check_offline_ready

        sources = tmp_path / "sources"
        sources.mkdir(parents=True)

        ready, error_msg = _check_offline_ready(sources, "nonexistent_species")
        assert ready is False
        assert "No germline data found" in error_msg
        assert "download_imgt.py" in error_msg
        assert "download_ogrdb.py" in error_msg

    def test_pipeline_raises_for_missing_species(self, tmp_path):
        from sadie.germlines.pipeline import GermlinePipeline

        sources = tmp_path / "sources"
        sources.mkdir(parents=True)

        pipeline = GermlinePipeline(tmp_path)

        with pytest.raises(RuntimeError) as excinfo:
            pipeline.update("nonexistent_species")

        assert "No germline data found" in str(excinfo.value)


class TestCachedDataUsage:
    def test_old_data_still_works(self, tmp_path):
        import os
        import time

        sources = tmp_path / "sources" / "imgt" / "human"
        sources.mkdir(parents=True)

        fasta = """>IGHV1-2*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        fasta_path = sources / "IGHV.fasta"
        fasta_path.write_text(fasta)

        old_time = time.time() - (180 * 24 * 60 * 60)
        os.utime(fasta_path, (old_time, old_time))

        from sadie.germlines.providers.imgt import IMGTProvider
        provider = IMGTProvider(data_dir=tmp_path / "sources" / "imgt")
        genes = provider.fetch_genes("human", "V", "H")
        assert len(genes) >= 1
