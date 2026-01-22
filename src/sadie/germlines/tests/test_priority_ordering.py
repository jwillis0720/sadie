import pytest
from pathlib import Path
import tempfile
import shutil


class TestPriorityOrdering:
    @pytest.fixture
    def priority_env(self, tmp_path):
        sources = tmp_path / "sources"
        for provider in ["custom", "imgt", "ogrdb", "vdjbase"]:
            (sources / provider / "human").mkdir(parents=True)

        imgt_v = """>X00001|IGHV1-69*01|Homo sapiens|F|V-REGION
CAGGTGCAGCTGGTGCAGTCTGGGGCTAAAA
>X00002|IGHV1-2*01|Homo sapiens|F|V-REGION
CAGGTGCAGCTGGTGCAGTCTGGGGCTTTTT
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        ogrdb_v = """>IGHV1-69*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTCCCC
>IGHV-OGRDB-NOVEL*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTACAC
"""
        (sources / "ogrdb" / "human" / "IGHV.fasta").write_text(ogrdb_v)
        (sources / "ogrdb" / "human" / "IGHV_gapped.fasta").write_text(ogrdb_v)

        custom_v = """>IGHV1-69*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGGGG
"""
        (sources / "custom" / "human" / "IGHV.fasta").write_text(custom_v)

        return sources

    def test_default_priority_custom_over_ogrdb_over_imgt(self, priority_env):
        from sadie.germlines.providers.custom import CustomProvider
        from sadie.germlines.providers.imgt import IMGTProvider
        from sadie.germlines.providers.ogrdb import OGRDBProvider

        custom = CustomProvider(data_dir=priority_env / "custom")
        ogrdb = OGRDBProvider(data_dir=priority_env / "ogrdb")
        imgt = IMGTProvider(data_dir=priority_env / "imgt")

        custom_genes = custom.fetch_genes("human", "V", "H")
        ogrdb_genes = ogrdb.fetch_genes("human", "V", "H")
        imgt_genes = imgt.fetch_genes("human", "V", "H")

        custom_69 = next((g for g in custom_genes if "IGHV1-69" in g.name), None)
        ogrdb_69 = next((g for g in ogrdb_genes if "IGHV1-69" in g.name), None)
        imgt_69 = next((g for g in imgt_genes if "IGHV1-69" in g.name), None)

        assert custom_69 is not None
        assert ogrdb_69 is not None
        assert imgt_69 is not None

        assert custom_69.sequence.endswith("GGGG")
        assert ogrdb_69.sequence.endswith("CCCC")
        assert imgt_69.sequence.endswith("AAAA")

    def test_same_name_different_sequence_uses_higher_priority(self, priority_env):
        from sadie.germlines.providers.custom import CustomProvider
        from sadie.germlines.providers.ogrdb import OGRDBProvider

        custom = CustomProvider(data_dir=priority_env / "custom")
        ogrdb = OGRDBProvider(data_dir=priority_env / "ogrdb")

        custom_genes = custom.fetch_genes("human", "V", "H")
        ogrdb_genes = ogrdb.fetch_genes("human", "V", "H")

        custom_69 = next((g for g in custom_genes if "IGHV1-69" in g.name), None)
        ogrdb_69 = next((g for g in ogrdb_genes if "IGHV1-69" in g.name), None)

        assert custom_69.sequence != ogrdb_69.sequence

    def test_gene_in_imgt_not_ogrdb_included(self, priority_env):
        from sadie.germlines.providers.imgt import IMGTProvider
        from sadie.germlines.providers.ogrdb import OGRDBProvider

        imgt = IMGTProvider(data_dir=priority_env / "imgt")
        ogrdb = OGRDBProvider(data_dir=priority_env / "ogrdb")

        imgt_genes = imgt.fetch_genes("human", "V", "H")
        ogrdb_genes = ogrdb.fetch_genes("human", "V", "H")

        imgt_names = {g.name for g in imgt_genes}
        ogrdb_names = {g.name for g in ogrdb_genes}

        imgt_only_gene = next((n for n in imgt_names if "IGHV1-2" in n), None)
        assert imgt_only_gene is not None
        assert imgt_only_gene not in ogrdb_names

    def test_novel_genes_from_lower_priority_included(self, priority_env):
        from sadie.germlines.providers.ogrdb import OGRDBProvider
        from sadie.germlines.providers.imgt import IMGTProvider

        ogrdb = OGRDBProvider(data_dir=priority_env / "ogrdb")
        imgt = IMGTProvider(data_dir=priority_env / "imgt")

        ogrdb_genes = ogrdb.fetch_genes("human", "V", "H")
        imgt_genes = imgt.fetch_genes("human", "V", "H")

        ogrdb_names = {g.name for g in ogrdb_genes}
        imgt_names = {g.name for g in imgt_genes}

        ogrdb_novel = next((n for n in ogrdb_names if "NOVEL" in n), None)
        assert ogrdb_novel is not None
        assert ogrdb_novel not in imgt_names

    def test_same_name_same_sequence_keeps_one(self, tmp_path):
        sources = tmp_path / "sources"
        (sources / "imgt" / "human").mkdir(parents=True)
        (sources / "ogrdb" / "human").mkdir(parents=True)

        same_seq = "CAGGTGCAGCTGGTGCAGTCTGGGGCTAACC"

        imgt_v = f""">X00001|IGHV1-69*01|Homo sapiens|F|V-REGION
{same_seq}
"""
        (sources / "imgt" / "human" / "IGHV.fasta").write_text(imgt_v)
        (sources / "imgt" / "human" / "IGHV_gapped.fasta").write_text(imgt_v)

        ogrdb_v = f""">IGHV1-69*01
{same_seq}
"""
        (sources / "ogrdb" / "human" / "IGHV.fasta").write_text(ogrdb_v)
        (sources / "ogrdb" / "human" / "IGHV_gapped.fasta").write_text(ogrdb_v)

        from sadie.germlines.providers.imgt import IMGTProvider
        from sadie.germlines.providers.ogrdb import OGRDBProvider

        imgt = IMGTProvider(data_dir=sources / "imgt")
        ogrdb = OGRDBProvider(data_dir=sources / "ogrdb")

        imgt_genes = imgt.fetch_genes("human", "V", "H")
        ogrdb_genes = ogrdb.fetch_genes("human", "V", "H")

        imgt_69 = next((g for g in imgt_genes if "IGHV1-69" in g.name), None)
        ogrdb_69 = next((g for g in ogrdb_genes if "IGHV1-69" in g.name), None)

        assert imgt_69 is not None
        assert ogrdb_69 is not None
        assert imgt_69.sequence == ogrdb_69.sequence
