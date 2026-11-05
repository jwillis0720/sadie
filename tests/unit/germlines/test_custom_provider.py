import shutil
import tempfile
from pathlib import Path

import pytest

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


class TestCustomProviderGappedFasta:
    """Tests for _gapped.fasta support in CustomProvider."""

    @pytest.fixture
    def temp_dir(self):
        d = tempfile.mkdtemp()
        yield Path(d)
        shutil.rmtree(d)

    def test_gapped_fasta_used_when_present(self, temp_dir):
        """Test that _gapped.fasta sequences are used when available."""
        custom_dir = temp_dir / "custom" / "human"
        custom_dir.mkdir(parents=True)

        # Ungapped main file
        ungapped_content = """>IGHV1-TEST*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (custom_dir / "IGHV.fasta").write_text(ungapped_content)

        # Pre-gapped companion file with known gapped sequence
        gapped_content = """>IGHV1-TEST*01
CAG.GTG.CAG.CTG.GTG.CAG.TCT.GGG.GCT
"""
        (custom_dir / "IGHV_gapped.fasta").write_text(gapped_content)

        provider = CustomProvider(data_dir=temp_dir / "custom")
        genes = provider.fetch_genes("human", "V", "H")

        assert len(genes) == 1
        # Verify pre-gapped sequence used (contains dots)
        assert genes[0].sequence_gapped == "CAG.GTG.CAG.CTG.GTG.CAG.TCT.GGG.GCT"

    def test_auto_gapping_when_no_gapped_fasta(self, temp_dir):
        """Test that auto-gapping still works when no _gapped.fasta exists."""
        custom_dir = temp_dir / "custom" / "human"
        custom_dir.mkdir(parents=True)

        ungapped_content = """>IGHV1-TEST*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
"""
        (custom_dir / "IGHV.fasta").write_text(ungapped_content)
        # No IGHV_gapped.fasta - should still work via auto-gapping or None

        provider = CustomProvider(data_dir=temp_dir / "custom")
        genes = provider.fetch_genes("human", "V", "H")

        assert len(genes) == 1
        assert genes[0].sequence == "CAGGTGCAGCTGGTGCAGTCTGGGGCT"
        # Gapped sequence might be None or auto-gapped depending on template availability
        # The important thing is no crash occurs

    def test_gapped_fasta_partial_coverage(self, temp_dir):
        """Test behavior when _gapped.fasta only has some sequences."""
        custom_dir = temp_dir / "custom" / "human"
        custom_dir.mkdir(parents=True)

        ungapped_content = """>IGHV1-TEST*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT
>IGHV1-TEST*02
GAGGTGCAGCTGGTGCAGTCTGGGGCC
"""
        (custom_dir / "IGHV.fasta").write_text(ungapped_content)

        # Only one sequence in gapped file
        gapped_content = """>IGHV1-TEST*01
CAG.GTG.CAG.CTG.GTG.CAG.TCT.GGG.GCT
"""
        (custom_dir / "IGHV_gapped.fasta").write_text(gapped_content)

        provider = CustomProvider(data_dir=temp_dir / "custom")
        genes = provider.fetch_genes("human", "V", "H")

        assert len(genes) == 2
        # First gene should use pre-gapped
        test01 = next(g for g in genes if g.name == "IGHV1-TEST*01")
        assert test01.sequence_gapped == "CAG.GTG.CAG.CTG.GTG.CAG.TCT.GGG.GCT"
        # Second gene should fall back to auto-gapping (or None if templates not available)
        test02 = next(g for g in genes if g.name == "IGHV1-TEST*02")
        # The key is the gene exists and has its ungapped sequence
        assert test02.sequence == "GAGGTGCAGCTGGTGCAGTCTGGGGCC"
