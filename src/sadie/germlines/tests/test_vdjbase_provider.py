import pytest
from pathlib import Path
import tempfile
import shutil


class TestVDJbaseProvider:
    @pytest.fixture
    def vdjbase_env(self, tmp_path):
        vdjbase_dir = tmp_path / "vdjbase" / "human"
        vdjbase_dir.mkdir(parents=True)

        fasta = """>IGHV1-69*01|VDJbase
CAGGTGCAGCTGGTGCAGTCTGGGGCTAAAA
>IGHV1-69*02|VDJbase
CAGGTGCAGCTGGTGCAGTCTGGGGCTAACC
>IGHV-VDJBASE-NOVEL*01|VDJbase
CAGGTGCAGCTGGTGCAGTCTGGGGCTCCCC
"""
        (vdjbase_dir / "IGHV.fasta").write_text(fasta)

        j_fasta = """>IGHJ4*01|VDJbase
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
"""
        (vdjbase_dir / "IGHJ.fasta").write_text(j_fasta)

        return tmp_path / "vdjbase"

    def test_fetch_genes_returns_genes(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        genes = provider.fetch_genes("human", "V", "H")

        assert len(genes) == 3
        names = [g.name for g in genes]
        assert any("IGHV1-69*01" in n for n in names)
        assert any("IGHV1-69*02" in n for n in names)
        assert any("NOVEL" in n for n in names)

    def test_genes_marked_vdjbase_source(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        genes = provider.fetch_genes("human", "V", "H")

        for gene in genes:
            assert gene.source == "vdjbase"

    def test_is_available(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        assert provider.is_available("human") is True
        assert provider.is_available("nonexistent") is False

    def test_fetch_gene_by_name(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        gene = provider.fetch_gene_by_name("IGHV1-69*01", "human")

        assert gene is not None
        assert "IGHV1-69*01" in gene.name
        assert gene.sequence.endswith("AAAA")

    def test_fetch_gene_by_name_not_found(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        gene = provider.fetch_gene_by_name("NONEXISTENT*01", "human")

        assert gene is None

    def test_empty_directory(self, tmp_path):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        empty_dir = tmp_path / "empty_vdjbase"
        empty_dir.mkdir(parents=True)

        provider = VDJbaseProvider(data_dir=empty_dir)
        genes = provider.fetch_genes("human", "V", "H")

        assert genes == []

    def test_get_metadata(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        metadata = provider.get_metadata()

        assert metadata.name == "vdjbase"
        assert "human" in metadata.species_available

    def test_j_genes_returned(self, vdjbase_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider

        provider = VDJbaseProvider(data_dir=vdjbase_env)
        genes = provider.fetch_genes("human", "J", "H")

        assert len(genes) >= 1
        assert any("IGHJ4" in g.name for g in genes)


class TestVDJbaseInPriority:
    @pytest.fixture
    def multi_provider_env(self, tmp_path):
        sources = tmp_path / "sources"

        (sources / "vdjbase" / "human").mkdir(parents=True)
        (sources / "imgt" / "human").mkdir(parents=True)

        vdjbase_v = """>IGHV1-69*01|VDJbase
CAGGTGCAGCTGGTGCAGTCTGGGGCTGGGG
"""
        (sources / "vdjbase" / "human" / "IGHV.fasta").write_text(vdjbase_v)

        imgt_v = """>X00001|IGHV1-69*01|Homo sapiens|F|V-REGION
CAGGTGCAGCTGGTGCAGTCTGGGGCTAAAA
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        return sources

    def test_vdjbase_takes_priority_when_first(self, multi_provider_env):
        from sadie.germlines.providers.vdjbase import VDJbaseProvider
        from sadie.germlines.providers.imgt import IMGTProvider

        vdjbase = VDJbaseProvider(data_dir=multi_provider_env / "vdjbase")
        imgt = IMGTProvider(data_dir=multi_provider_env / "imgt")

        vdjbase_genes = vdjbase.fetch_genes("human", "V", "H")
        imgt_genes = imgt.fetch_genes("human", "V", "H")

        vdjbase_69 = next((g for g in vdjbase_genes if "IGHV1-69" in g.name), None)
        imgt_69 = next((g for g in imgt_genes if "IGHV1-69" in g.name), None)

        assert vdjbase_69 is not None
        assert imgt_69 is not None
        assert vdjbase_69.sequence.endswith("GGGG")
        assert imgt_69.sequence.endswith("AAAA")
