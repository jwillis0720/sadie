import tempfile
import pytest
from sadie.anarci import Anarci, AnarciResult, AnarciResults
from sadie.antibody import exception
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def test_single_seq():
    anarci_api = Anarci(region_assign="imgt")
    result = anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    )
    # assert a bunch of numbering
    assert result.id == "MySweetAntibody"
    assert result.leader == "AAAADAFA"
    assert result.framework1_aa == "EVQLVESGG-GLEQPGGSLRLSCAGS"
    assert result.framework2_aa == "MTWVRQAPGKGLEWVSS"
    assert result.framework3_aa == "YYADSVK-GRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    assert result.framework4_aa == "WGQGTTVTVSS"
    assert result.cdr1_aa == "GFTF----RDYA"
    assert result.cdr2_aa == "ISGS--GGNT"
    assert result.cdr3_aa == "AKDRLSITIRPRYYGLDV"
    assert result.tail == "RRRESV"
    assert result.j_gene == "IGHJ6*01"
    assert result.v_gene == "IGHV3-23*04"
    assert result.scheme == "imgt"
    assert result.v_gene_identity == 0.93
    assert result.j_gene_identity == 0.93

    # get json string
    json_string = result.get_json()
    result_2 = AnarciResult.from_json(json_string)
    assert result == result_2

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json.gz") as tmpfile:
        result.to_json(tmpfile.name)
        result3 = AnarciResult.read_json(tmpfile.name)
        assert result == result3


def test_alternate_numbering():
    anarci_api = Anarci(scheme="chothia", region_assign="chothia")
    anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    )
    # chothia_result = result_chothia.to_antibody()
    # align = "IGHV3-23*04|IGHJ6*01  EVQLVESGGGLVQPGGSLRLSCAAS GFTFSSY AMSWVRQAPGKGLEWVSAI SGSGGS TYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC\nMySweetAntibody       ...........E...........G. ....RD. ..T..............S. .....N .......................................\n\nIGHV3-23*04|IGHJ6*01  A- K------YYYYYGMDV WGQGTTVTVSS\nMySweetAntibody       .K DRLSITIRPR...L.. ...........\n\n"
    # assert chothia_result.get_segmented_alignment_aa() == align
    anarci_api = Anarci(scheme="chothia", region_assign="abm")
    anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    )
    # align = "IGHV3-23*04|IGHJ6*01  EVQLVESGGGLVQPGGSLRLSCAAS GFTFSSYAMS WVRQAPGKGLEWVS AISGSGGSTY YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC\nMySweetAntibody       ...........E...........G. ....RD...T .............. S......N.. .....................................\n\nIGHV3-23*04|IGHJ6*01  A- K------YYYYYGMDV WGQGTTVTVSS\nMySweetAntibody       .K DRLSITIRPR...L.. ...........\n\n"
    # abm_result = result_abm.to_antibody()
    # assert abm_result.get_segmented_alignment_aa() == align


def test_long_hcdr3():
    anarci_api = Anarci(scheme="chothia", region_assign="chothia")
    with pytest.raises(exception.LongHCDR3Error) as e:
        anarci_api.run_single(
            "MySweetAntibody",
            "QVQLVESGGGVVQPGRSLRLSCAASGFTFNNYGMHWVRQAPGKGLEWVAVISYGGSDKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDGRGSLPRPKGGFINALSFHWPFGRWLGKSYGTYDSSEDSGGAFDIWGQGTLVTVSS",
        )
    assert e.value.hcdr3 == "ARDGRGSLPRPKGGFINALSFHWPFGRWLGKSYGTYDSSEDSGGAFDI"
    assert e.value.chosen_scheme == "chothia"
    assert e.value.acceptable_scheme == ["imgt", "aho"]

    anarci_api = Anarci(scheme="kabat", region_assign="chothia")
    with pytest.raises(exception.LongHCDR3Error) as e:
        anarci_api.run_single(
            "MySweetAntibody",
            "QVQLVESGGGVVQPGRSLRLSCAASGFTFNNYGMHWVRQAPGKGLEWVAVISYGGSDKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDGRGSLPRPKGGFINALSFHWPFGRWLGKSYGTYDSSEDSGGAFDIWGQGTLVTVSS",
        )
    assert e.value.hcdr3 == "ARDGRGSLPRPKGGFINALSFHWPFGRWLGKSYGTYDSSEDSGGAFDI"
    assert e.value.chosen_scheme == "kabat"
    assert e.value.acceptable_scheme == ["imgt", "aho"]
    assert e.value.sequence_name == "MySweetAntibody"


def test_anarci_multi_input():
    anarci_api = Anarci()

    seq_records = [
        SeqRecord(
            Seq(
                "EVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSS"
            ),
            id="DupulimabH",
        ),
        SeqRecord(
            Seq(
                "DIVMTQSPLSLPVTPGEPASISCRSSQSLLYSIGYNYLDWYLQKSGQSPQLLIYLGSNRASGVPDRFSGSGSGTDFTLKISRVEAEDVGFYYCMQALQTPYTFGQGTKLEIK"
            ),
            id="DupulimabL",
        ),
    ]
    results = anarci_api.run_multiple(seq_records)
    assert isinstance(results, AnarciResults)
    single_result_1 = results["DupulimabH"]
    single_result_2 = results["DupulimabL"]
    assert isinstance(single_result_1, AnarciResult)
    assert isinstance(single_result_2, AnarciResult)
    assert single_result_1.framework1_aa == "EVQLVESGG-GLEQPGGSLRLSCAGS"
    assert single_result_1.cdr1_aa == "GFTF----RDYA"
    assert single_result_1.framework2_aa == "MTWVRQAPGKGLEWVSS"
    assert single_result_1.cdr2_aa == "ISGS--GGNT"
    assert single_result_1.framework3_aa == "YYADSVK-GRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    assert single_result_1.cdr3_aa == "AKDRLSITIRPRYYGLDV"
    assert single_result_1.framework4_aa == "WGQGTTVTVSS"

    _segment_df = results.segment_table
    _summary_table = results.summary_table
    _alignment_table = results.alignment_table
    reconstructed_results = AnarciResults(
        summary_table=_summary_table,
        alignment_table=_alignment_table,
        segment_table=_segment_df,
    )
    assert reconstructed_results == results

    with tempfile.NamedTemporaryFile(suffix=".json.gz") as f:
        results.to_json(f.name)
        new_results = AnarciResults.read_json(f.name)
        assert new_results == results

    assert results == AnarciResults.from_json(results.get_json())


def test_io():
    """Test file io"""
    main_file = fixture_file("split_files/scfv_heavy_1.fasta")
    anarci_api = Anarci(allowed_species=["dog", "cat"])
    results = anarci_api.run_file(main_file)
    assert isinstance(results, AnarciResults)
    with tempfile.NamedTemporaryFile(suffix=".anarci.bz2") as temp:
        results.to_file(temp.name)
        other_results = AnarciResults.read_file(temp.name)
        assert other_results == results


def test_dog():
    anarci_api = Anarci()
    result = anarci_api.run_single(
        "V3-38_J4",
        "EVQLVESGGDLVKPGGTLRLSCVASGLSLTSNSMSWVRQSPGKGLQWVAVIWSNGGTYYADAVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCASIYYYDADYLHWGQGTLVTVSS",
    )
    assert result.framework1_aa == "EVQLVESGG-DLVKPGGTLRLSCVAS"
    assert result.cdr1_aa == "GLSL----TSNS"
    assert result.framework2_aa == "MSWVRQSPGKGLQWVAV"
    assert result.cdr2_aa == "IWSN---GGT"
    assert result.framework3_aa == "YYADAVK-GRFTISRDNAKNTLYLQMNSLRAEDTAVYYC"
    assert result.cdr3_aa == "ASIYYY-DADYLH"
    assert result.framework4_aa == "WGQGTLVTVSS"
    assert result.v_gene == "IGHV3-38*01"
    assert result.j_gene == "IGHJ2*01"
    assert result.v_gene_identity == 0.88
    assert result.j_gene_identity == 0.79

    result = anarci_api.run_single(
        "V3-38_J4",
        "EIVMTQSPASLSLSQEEKVTITCRASEGISNSLAWYQQKPGQAPKLLIYATSNRATGVPSRFSGSGSGTDFSFTISSLEPEDVAVYYCQQGYKFPLTFGAGTKVELK",
    )
    assert result.framework1_aa == "EIVMTQSPASLSLSQEEKVTITCRAS"
    assert result.cdr1_aa == "EGI------SNS"
    assert result.framework2_aa == "LAWYQQKPGQAPKLLIY"
    assert result.cdr2_aa == "AT-------S"
    assert result.framework3_aa == "NRATGVP-SRFSGSG--SGTDFSFTISSLEPEDVAVYYC"
    assert result.cdr3_aa == "QQGYK----FPLT"
    assert result.framework4_aa == "FGAGTKVELK"
    assert result.v_gene == "IGKV3S1*01"
    assert result.j_gene == "IGKJ1*01"
    assert result.v_gene_identity == 0.91
    assert result.j_gene_identity == 0.92


def test_cat():
    anarci_api = Anarci(allowed_species=["cat"])
    result = anarci_api.run_single(
        "CF-R01-D01",
        "DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF",
    )
    assert result.id == "CF-R01-D01"
    assert result.framework1_aa == "DVQLVESGG-DLAKPGGSLRLTCVAS"
    assert result.cdr1_aa == "GLSV----TSNS"
    assert result.framework2_aa == "MSWVRQAPGKGLRWVST"
    assert result.cdr2_aa == "IWSK---GGT"
    assert result.framework3_aa == "YYADSVK-GRFTVSRDSAKNTLYLQMDSLATEDTATYYC"
    assert result.framework4_aa == "WGQGALVTVSF"
    assert (
        result.vdj
        == "DVQLVESGG-DLAKPGGSLRLTCVASGLSV----TSNSMSWVRQAPGKGLRWVSTIWSK---GGTYYADSVK-GRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF"
    )
    assert result.v_gene == "IGHV17-1*01"
    assert result.j_gene == "IGHJ5*01"
    assert result.v_gene_identity == 0.84
    assert result.j_gene_identity == 0.79


def test_bad_sequences():
    bad_sequence = "SEQLTQPESLTLRPGQPLTIRCQVSYSVSSTGYATHWIRQPDGRGLEWIGGIRIGWKGAKDSLSSQFSLAVDGSSKTITLQGQNMQPGDSAVYYCAR"
    anarci_api = Anarci()
    result = anarci_api.run_single("bad_sequence", bad_sequence)
    assert result is None
