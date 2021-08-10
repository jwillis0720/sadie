import glob
import logging
import os
import tempfile
from itertools import product
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner
from pandas.testing import assert_frame_equal
from sadie.anarci import Anarci, AnarciDuplicateIdError, AnarciResults
from sadie.anarci.app import run_anarci

logger = logging.getLogger()


# def get_file(file):
#     """Helper method for test execution."""
#     _file = os.path.join(os.path.abspath(os.path.dirname(__file__)), f"fixtures/{file}")
#     if not os.path.exists(_file):
#         raise FileNotFoundError(_file)
#     return _file


# def fixture_file(file):
#     """Helper method for test execution."""
#     return resource_filename(__name__, "fixtures/{}".format(file))


def test_long_seq():
    anarci_api = Anarci(scheme="chothia", region_assign="imgt", allowed_species=["human"])
    anarci_api.run_single(
        "VRC26.27_KT371104_Homo_sapiens_anti-HIV-1_immunoglobulin",
        "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS",
    )
    anarci_api = Anarci(scheme="kabat", region_assign="imgt", allowed_species=["human"])
    anarci_api.run_single(
        "VRC26.27_KT371104_Homo_sapiens_anti-HIV-1_immunoglobulin",
        "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS",
    )


def test_no_j_gene():
    """no j gene found"""
    anarci_api = Anarci(scheme="chothia", region_assign="imgt", allowed_species=["rat"])
    anarci_api.run_single(
        "VRC26.27_KT371104_Homo_sapiens_anti-HIV-1_immunoglobulin",
        "QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS",
    )


def test_trouble_seqs():
    anarci = Anarci(scheme="kabat", region_assign="imgt")
    with pytest.warns(UserWarning):
        results = anarci.run_single(
            "POS",
            "DIQMTQSPSSLCASIGDRVTITCRASQSISSYLNQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPVTFGGGTKVEIK",
        )
        assert results.empty

    # check that we can still get out good sequence in presense of the bad
    seq_records = [
        SeqRecord(
            Seq("DIQMTQSPSSLCASIGDRVTITCRASQSISSYLNQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPVTFGGGTKVEIK"),
            id="POS",
        ),
        SeqRecord(
            Seq(
                "DIVMTQSPLSLPVTPGEPASISCRSSQSLLYSIGYNYLDWYLQKSGQSPQLLIYLGSNRASGVPDRFSGSGSGTDFTLKISRVEAEDVGFYYCMQALQTPYTFGQGTKLEIK"
            ),
            id="DupulimabL",
        ),
    ]
    results = anarci.run_multiple(seq_records)

    # can't capture this user warning with multiprocess
    anarci = Anarci(scheme="kabat", region_assign="imgt", run_multiproc=False)
    with pytest.warns(UserWarning):
        results = anarci.run_multiple(seq_records)
    assert len(results) == 1


def test_single_seq():
    anarci_api = Anarci(region_assign="imgt")
    result = anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    ).iloc[0]
    # assert a bunch of numbering

    assert result.Id == "MySweetAntibody"
    assert result.leader == "AAAADAFA"
    assert result.fwr1_aa_gaps == "EVQLVESGG-GLEQPGGSLRLSCAGS"
    assert result.fwr2_aa_gaps == "MTWVRQAPGKGLEWVSS"
    assert result.fwr3_aa_gaps == "YYADSVK-GRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    assert result.fwr4_aa_gaps == "WGQGTTVTVSS"
    assert result.cdr1_aa_gaps == "GFTF----RDYA"
    assert result.cdr2_aa_gaps == "ISGS--GGNT"
    assert result.cdr3_aa_gaps == "AKDRLSITIRPRYYGLDV"
    assert result.follow == "RRRESV"
    assert result.j_gene == "IGHJ6*01"
    assert result.v_gene == "IGHV3-23*04"
    assert result.scheme == "imgt"
    assert float(result.v_identity) == 0.93
    assert float(result.j_identity) == 0.93


def test_alternate_numbering():
    anarci_api = Anarci(scheme="chothia", region_assign="chothia")
    anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    )
    anarci_api = Anarci(scheme="chothia", region_assign="abm")
    anarci_api.run_single(
        "MySweetAntibody",
        "AAAADAFAEVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSSRRRESV",
    )


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
    single_result_1 = results.query("Id=='DupulimabH'").iloc[0]
    single_result_2 = results.query("Id=='DupulimabL'").iloc[0]
    assert single_result_1.fwr1_aa_gaps == "EVQLVESGG-GLEQPGGSLRLSCAGS"
    assert single_result_1.cdr1_aa_gaps == "GFTF----RDYA"
    assert single_result_1.fwr2_aa_gaps == "MTWVRQAPGKGLEWVSS"
    assert single_result_1.cdr2_aa_gaps == "ISGS--GGNT"
    assert single_result_1.fwr3_aa_gaps == "YYADSVK-GRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    assert single_result_1.cdr3_aa_gaps == "AKDRLSITIRPRYYGLDV"
    assert single_result_1.fwr4_aa_gaps == "WGQGTTVTVSS"

    assert single_result_2.fwr1_aa_gaps == "DIVMTQSPLSLPVTPGEPASISCRSS"
    assert single_result_2.cdr1_aa_gaps == "QSLLYS-IGYNY"
    assert single_result_2.fwr2_aa_gaps == "LDWYLQKSGQSPQLLIY"
    assert single_result_2.cdr2_aa_gaps == "LG-------S"
    assert single_result_2.fwr3_aa_gaps == "NRASGVP-DRFSGSG--SGTDFTLKISRVEAEDVGFYYC"
    assert single_result_2.cdr3_aa_gaps == "MQALQ----TPYT"
    assert single_result_2.fwr4_aa_gaps == "FGQGTKLEIK"


def test_io(fixture_setup):
    """Test file io"""
    main_file = fixture_setup.get_dog_aa_seqs()
    anarci_api = Anarci(allowed_species=["dog", "cat"])
    results = anarci_api.run_file(main_file)
    assert isinstance(results, AnarciResults)
    with tempfile.NamedTemporaryFile(suffix=".anarci.bz2") as temp:
        results.to_csv(temp.name)
        other_results = AnarciResults.read_csv(temp.name).fillna("")
        assert_frame_equal(other_results, results)


def test_dog():
    anarci_api = Anarci()
    result = anarci_api.run_single(
        "V3-38_J4",
        "EVQLVESGGDLVKPGGTLRLSCVASGLSLTSNSMSWVRQSPGKGLQWVAVIWSNGGTYYADAVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCASIYYYDADYLHWGQGTLVTVSS",
    ).iloc[0]
    assert result.fwr1_aa_gaps == "EVQLVESGG-DLVKPGGTLRLSCVAS"
    assert result.cdr1_aa_gaps == "GLSL----TSNS"
    assert result.fwr2_aa_gaps == "MSWVRQSPGKGLQWVAV"
    assert result.cdr2_aa_gaps == "IWSN---GGT"
    assert result.fwr3_aa_gaps == "YYADAVK-GRFTISRDNAKNTLYLQMNSLRAEDTAVYYC"
    assert result.cdr3_aa_gaps == "ASIYYY-DADYLH"
    assert result.fwr4_aa_gaps == "WGQGTLVTVSS"
    assert result.v_gene == "IGHV3-38*01"
    assert result.j_gene == "IGHJ2*01"
    assert float(result.v_identity) == 0.88
    assert float(result.j_identity) == 0.79

    result = anarci_api.run_single(
        "V3-38_J4",
        "EIVMTQSPASLSLSQEEKVTITCRASEGISNSLAWYQQKPGQAPKLLIYATSNRATGVPSRFSGSGSGTDFSFTISSLEPEDVAVYYCQQGYKFPLTFGAGTKVELK",
    ).iloc[0]
    assert result.fwr1_aa_gaps == "EIVMTQSPASLSLSQEEKVTITCRAS"
    assert result.cdr1_aa_gaps == "EGI------SNS"
    assert result.fwr2_aa_gaps == "LAWYQQKPGQAPKLLIY"
    assert result.cdr2_aa_gaps == "AT-------S"
    assert result.fwr3_aa_gaps == "NRATGVP-SRFSGSG--SGTDFSFTISSLEPEDVAVYYC"
    assert result.cdr3_aa_gaps == "QQGYK----FPLT"
    assert result.fwr4_aa_gaps == "FGAGTKVELK"
    assert result.v_gene == "IGKV3S1*01"
    assert result.j_gene == "IGKJ1*01"
    assert float(result.v_identity) == 0.91
    assert float(result.j_identity) == 0.92


def test_cat():
    anarci_api = Anarci(allowed_species=["cat"])
    result = anarci_api.run_single(
        "CF-R01-D01",
        "DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF",
    ).iloc[0]
    assert result.Id == "CF-R01-D01"
    assert result.fwr1_aa_gaps == "DVQLVESGG-DLAKPGGSLRLTCVAS"
    assert result.cdr1_aa_gaps == "GLSV----TSNS"
    assert result.fwr2_aa_gaps == "MSWVRQAPGKGLRWVST"
    assert result.cdr2_aa_gaps == "IWSK---GGT"
    assert result.fwr3_aa_gaps == "YYADSVK-GRFTVSRDSAKNTLYLQMDSLATEDTATYYC"
    assert result.fwr4_aa_gaps == "WGQGALVTVSF"
    assert result.v_gene == "IGHV17-1*01"
    assert result.j_gene == "IGHJ5*01"
    assert result.v_identity == 0.84
    assert result.j_identity == 0.79


def test_bad_sequences():
    bad_sequence = "SEQLTQPESLTLRPGQPLTIRCQVSYSVSSTGYATHWIRQPDGRGLEWIGGIRIGWKGAKDSLSSQFSLAVDGSSKTITLQGQNMQPGDSAVYYCAR"
    anarci_api = Anarci()
    result = anarci_api.run_single("bad_sequence", bad_sequence)
    assert result.empty


def test_duplicated_seq():
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
                "EVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSS"
            ),
            id="DupulimabH",
        ),
    ]
    with pytest.raises(AnarciDuplicateIdError):
        anarci_api.run_multiple(seq_records)


# cli helper
def _run_cli(p_tuple):
    with tempfile.NamedTemporaryFile() as tmpfile:
        cli_input = [
            "--query",
            str(p_tuple[0]),
            "-o",
            tmpfile.name,
            "-a",
            p_tuple[1],
            "-f",
            p_tuple[2],
            "-s",
            p_tuple[3],
            "-r",
            p_tuple[4],
            "-vvvvv",
        ]
        if not Anarci.check_combination(p_tuple[3], p_tuple[4]):
            logger.info(f"skipping {p_tuple[4]}-{p_tuple[3]}")
            return True

        logger.info(f"CLI input {' '.join(cli_input)}")
        runner = CliRunner()
        result = runner.invoke(run_anarci, cli_input)
        if result.exit_code != 0:
            logger.info(f"!!!STDERR {result}")

        assert result.exit_code == 0

        # anarci appends reults and alignment so glob
        out_file = Path(tmpfile.name).stem + "*"
    for f in glob.glob(out_file):
        os.remove(f)
    return True


def test_cli(fixture_setup):
    """Confirm the CLI works as expecte"""
    test_file_heavy = fixture_setup.get_catnap_heavy_aa()
    test_file_light = fixture_setup.get_catnap_light_aa()
    species = ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "dog", "cat"]
    species = ",".join(species)
    ft = ["csv", "json", "feather"]
    schemes = ["imgt", "kabat", "chothia"]
    regions = ["imgt", "kabat", "chothia", "abm", "contact", "scdr"]
    products = product([test_file_heavy, test_file_light], [species], ft, schemes, regions)

    # pool = Pool()
    results = list(map(_run_cli, products))
    print(results)


def test_df(fixture_setup):
    test_file_heavy = fixture_setup.get_catnap_heavy_aa()
    df = pd.DataFrame(
        [{"id": x.id, "seq": x.seq, "description": x.description} for x in SeqIO.parse(test_file_heavy, "fasta")]
    )
    anarci_obj = Anarci()
    anarci_results = anarci_obj.run_dataframe(df, "id", "seq")
    assert isinstance(anarci_results, (AnarciResults, pd.DataFrame))
