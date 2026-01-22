import pytest
from pathlib import Path
import tempfile
import shutil

from sadie.germlines.providers.custom import CustomProvider, _validate_sequence


class TestValidateSequence:
    def test_valid_nucleotides(self):
        valid, msg = _validate_sequence("ACGTACGT", "test_gene")
        assert valid is True
        assert msg == ""

    def test_valid_with_gaps(self):
        valid, msg = _validate_sequence("ACG...TGC", "test_gene")
        assert valid is True

    def test_valid_iupac_ambiguous(self):
        valid, msg = _validate_sequence("ACGTNRYSWKM", "test_gene")
        assert valid is True

    def test_empty_sequence(self):
        valid, msg = _validate_sequence("", "test_gene")
        assert valid is False
        assert "Empty" in msg

    def test_only_gaps(self):
        valid, msg = _validate_sequence("...", "test_gene")
        assert valid is False
        assert "gap characters" in msg

    def test_invalid_characters(self):
        valid, msg = _validate_sequence("ACGT123XYZ", "test_gene")
        assert valid is False
        assert "Invalid" in msg


class TestCustomProvider:
    @pytest.fixture
    def temp_dir(self):
        d = tempfile.mkdtemp()
        yield Path(d)
        shutil.rmtree(d)

    @pytest.fixture
    def provider_with_data(self, temp_dir):
        custom_dir = temp_dir / "custom" / "human"
        custom_dir.mkdir(parents=True)
        imgt_dir = temp_dir / "imgt" / "human"
        imgt_dir.mkdir(parents=True)

        fasta_content = """>IGHV1-TEST*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
>IGHV1-TEST*02
CAGGTGCAGCTGGTGCAGTCTGGGGCC
"""
        (custom_dir / "IGHV.fasta").write_text(fasta_content)

        return CustomProvider(data_dir=temp_dir / "custom", template_dir=imgt_dir.parent)

    def test_fetch_genes(self, provider_with_data):
        genes = provider_with_data.fetch_genes("human", "V", "H")
        assert len(genes) == 2
        assert genes[0].name == "IGHV1-TEST*01"
        assert genes[1].name == "IGHV1-TEST*02"

    def test_genes_marked_custom_source(self, provider_with_data):
        genes = provider_with_data.fetch_genes("human", "V", "H")
        for gene in genes:
            assert gene.source == "custom"

    def test_is_available(self, provider_with_data):
        assert provider_with_data.is_available("human") is True
        assert provider_with_data.is_available("mouse") is False

    def test_fetch_gene_by_name(self, provider_with_data):
        gene = provider_with_data.fetch_gene_by_name("IGHV1-TEST*01", "human")
        assert gene is not None
        assert gene.name == "IGHV1-TEST*01"

    def test_fetch_gene_by_name_not_found(self, provider_with_data):
        gene = provider_with_data.fetch_gene_by_name("NONEXISTENT", "human")
        assert gene is None

    def test_empty_directory(self, temp_dir):
        custom_dir = temp_dir / "custom"
        custom_dir.mkdir(parents=True)
        provider = CustomProvider(data_dir=custom_dir)
        genes = provider.fetch_genes("human", "V", "H")
        assert genes == []

    def test_invalid_sequence_skipped(self, temp_dir):
        custom_dir = temp_dir / "custom" / "human"
        custom_dir.mkdir(parents=True)

        fasta_content = """>VALID_GENE
ACGTACGT
>INVALID_GENE
123XYZ
"""
        (custom_dir / "IGHV.fasta").write_text(fasta_content)

        provider = CustomProvider(data_dir=temp_dir / "custom")
        genes = provider.fetch_genes("human", "V", "H")
        assert len(genes) == 1
        assert "VALID" in genes[0].name


class TestCustomProviderPriority:
    @pytest.fixture
    def setup_priority_test(self, tmp_path):
        from sadie.germlines.manager import GermlineManager

        custom_dir = tmp_path / "sources" / "custom" / "human"
        custom_dir.mkdir(parents=True)
        imgt_dir = tmp_path / "sources" / "imgt" / "human"
        imgt_dir.mkdir(parents=True)

        custom_fasta = """>IGHV1-SHARED*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTAAAA
"""
        (custom_dir / "IGHV.fasta").write_text(custom_fasta)

        imgt_fasta = """>IGHV1-SHARED*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTTTTT
>IGHV1-IMGTONLY*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (imgt_dir / "IGHV.fasta").write_text(imgt_fasta)

        return tmp_path / "sources"

    def test_custom_takes_priority_over_imgt(self, setup_priority_test):
        from sadie.germlines.providers.custom import CustomProvider
        from sadie.germlines.providers.imgt import IMGTProvider

        custom_provider = CustomProvider(data_dir=setup_priority_test / "custom")
        imgt_provider = IMGTProvider(data_dir=setup_priority_test / "imgt")

        custom_genes = custom_provider.fetch_genes("human", "V", "H")
        imgt_genes = imgt_provider.fetch_genes("human", "V", "H")

        custom_shared = next((g for g in custom_genes if "SHARED" in g.name), None)
        imgt_shared = next((g for g in imgt_genes if "SHARED" in g.name), None)

        assert custom_shared is not None
        assert imgt_shared is not None
        assert custom_shared.sequence.endswith("AAAA")
        assert imgt_shared.sequence.endswith("TTTT")
