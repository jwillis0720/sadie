import pytest
from pathlib import Path

REGRESSION_DIR = Path(__file__).parent / "data" / "regression"

def test_regression_sequences_exist():
    fasta = REGRESSION_DIR / "FR025b_sequences.fasta"
    assert fasta.exists(), f"Regression test file missing: {fasta}"

def test_fr025b_ighv1_69_present():
    fasta = REGRESSION_DIR / "FR025b_sequences.fasta"
    content = fasta.read_text()
    assert "IGHV1-69*01" in content

def test_fr025b_ighv3_23_present():
    fasta = REGRESSION_DIR / "FR025b_sequences.fasta"
    content = fasta.read_text()
    assert "IGHV3-23*01" in content

def test_fr025b_ighd3_3_present():
    fasta = REGRESSION_DIR / "FR025b_sequences.fasta"
    content = fasta.read_text()
    assert "IGHD3-3*01" in content

def test_fr025b_ighj4_present():
    fasta = REGRESSION_DIR / "FR025b_sequences.fasta"
    content = fasta.read_text()
    assert "IGHJ4*01" in content

@pytest.fixture
def germline_manager():
    from sadie.germlines.manager import GermlineManager
    return GermlineManager(providers=["imgt"])

def test_germlines_returns_ighv1_69(germline_manager):
    genes = germline_manager.get_genes("human", "V", "H")
    gene_names = [g.name for g in genes]
    assert any("IGHV1-69" in name for name in gene_names)

def test_germlines_returns_ighv3_23(germline_manager):
    genes = germline_manager.get_genes("human", "V", "H")
    gene_names = [g.name for g in genes]
    assert any("IGHV3-23" in name for name in gene_names)

def test_germlines_returns_ighd3_3(germline_manager):
    genes = germline_manager.get_genes("human", "D", "H")
    gene_names = [g.name for g in genes]
    assert any("IGHD3-3" in name for name in gene_names)

def test_germlines_returns_ighj4(germline_manager):
    genes = germline_manager.get_genes("human", "J", "H")
    gene_names = [g.name for g in genes]
    assert any("IGHJ4" in name for name in gene_names)

def test_sequence_identity_ighv1_69(germline_manager):
    genes = germline_manager.get_genes("human", "V", "H")
    ighv1_69 = next((g for g in genes if "IGHV1-69*01" in g.name), None)
    assert ighv1_69 is not None
    expected_start = "caggtgcagctggtgcagtctggggct"
    actual = ighv1_69.sequence.replace(".", "").lower()
    assert actual.startswith(expected_start)

def test_sequence_identity_ighj4(germline_manager):
    genes = germline_manager.get_genes("human", "J", "H")
    ighj4 = next((g for g in genes if "IGHJ4*01" in g.name), None)
    assert ighj4 is not None
    expected = "actactttgactactggggccaaggaaccctggtcaccgtctcctcag"
    actual = ighj4.sequence.replace(".", "").lower()
    assert actual == expected

def test_cdr_fwr_boundaries_within_tolerance():
    pass
