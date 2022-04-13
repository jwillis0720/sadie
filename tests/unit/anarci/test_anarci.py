import tempfile

from pyinstrument import Profiler
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.testing import assert_frame_equal
from sadie.anarci import Anarci, AnarciDuplicateIdError, AnarciResults
from sadie.anarci.aa._schemes import number_imgt
from sadie.anarci.aa._anarci import number_sequences_from_alignment, number_sequence_from_alignment


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
    anarci_api = Anarci(run_multiproc=True)
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


def test_df(fixture_setup):
    test_file_heavy = fixture_setup.get_catnap_heavy_aa()
    df = pd.DataFrame(
        [{"id": x.id, "seq": x.seq, "description": x.description} for x in SeqIO.parse(test_file_heavy, "fasta")]
    )
    anarci_obj = Anarci()
    anarci_results = anarci_obj.run_dataframe(df, "id", "seq")
    assert isinstance(anarci_results, (AnarciResults, pd.DataFrame))
    print("ehllo")


def test_numbering_seqs():
    sequences = [
        (
            "CF-R01-D01",
            "DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF",
        )
    ]
    _alignments = [
        (
            [
                ["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"],
                ["cat_H", "", 4.8e-57, 183.0, 1.8, 0, 122],
                ["dog_H", "", 2.4e-55, 177.5, 2.0, 0, 122],
                ["alpaca_H", "", 3.5e-54, 173.9, 2.0, 0, 122],
                ["human_H", "", 2.8e-51, 164.4, 1.6, 0, 122],
                ["pig_H", "", 3.7e-51, 164.2, 2.0, 1, 122],
                ["mouse_H", "", 1.3e-47, 152.4, 1.0, 0, 122],
                ["rabbit_H", "", 1.9e-46, 148.7, 1.2, 2, 122],
                ["cow_H", "", 2.8e-37, 118.9, 1.0, 1, 122],
                ["rhesus_H", "", 1.4e-34, 110.1, 2.1, 1, 122],
            ],
            [
                [
                    ((1, "m"), 0),
                    ((2, "m"), 1),
                    ((3, "m"), 2),
                    ((4, "m"), 3),
                    ((5, "m"), 4),
                    ((6, "m"), 5),
                    ((7, "m"), 6),
                    ((8, "m"), 7),
                    ((9, "m"), 8),
                    ((10, "d"), None),
                    ((11, "m"), 9),
                    ((12, "m"), 10),
                    ((13, "m"), 11),
                    ((14, "m"), 12),
                    ((15, "m"), 13),
                    ((16, "m"), 14),
                    ((17, "m"), 15),
                    ((18, "m"), 16),
                    ((19, "m"), 17),
                    ((20, "m"), 18),
                    ((21, "m"), 19),
                    ((22, "m"), 20),
                    ((23, "m"), 21),
                    ((24, "m"), 22),
                    ((25, "m"), 23),
                    ((26, "m"), 24),
                    ((27, "m"), 25),
                    ((28, "m"), 26),
                    ((29, "m"), 27),
                    ((30, "m"), 28),
                    ((31, "d"), None),
                    ((32, "d"), None),
                    ((33, "d"), None),
                    ((34, "d"), None),
                    ((35, "m"), 29),
                    ((36, "m"), 30),
                    ((37, "m"), 31),
                    ((38, "m"), 32),
                    ((39, "m"), 33),
                    ((40, "m"), 34),
                    ((41, "m"), 35),
                    ((42, "m"), 36),
                    ((43, "m"), 37),
                    ((44, "m"), 38),
                    ((45, "m"), 39),
                    ((46, "m"), 40),
                    ((47, "m"), 41),
                    ((48, "m"), 42),
                    ((49, "m"), 43),
                    ((50, "m"), 44),
                    ((51, "m"), 45),
                    ((52, "m"), 46),
                    ((53, "m"), 47),
                    ((54, "m"), 48),
                    ((55, "m"), 49),
                    ((56, "m"), 50),
                    ((57, "m"), 51),
                    ((58, "m"), 52),
                    ((59, "m"), 53),
                    ((60, "d"), None),
                    ((61, "d"), None),
                    ((62, "d"), None),
                    ((63, "m"), 54),
                    ((64, "m"), 55),
                    ((65, "m"), 56),
                    ((66, "m"), 57),
                    ((67, "m"), 58),
                    ((68, "m"), 59),
                    ((69, "m"), 60),
                    ((70, "m"), 61),
                    ((71, "m"), 62),
                    ((72, "m"), 63),
                    ((73, "d"), None),
                    ((74, "m"), 64),
                    ((75, "m"), 65),
                    ((76, "m"), 66),
                    ((77, "m"), 67),
                    ((78, "m"), 68),
                    ((79, "m"), 69),
                    ((80, "m"), 70),
                    ((81, "m"), 71),
                    ((82, "m"), 72),
                    ((83, "m"), 73),
                    ((84, "m"), 74),
                    ((85, "m"), 75),
                    ((86, "m"), 76),
                    ((87, "m"), 77),
                    ((88, "m"), 78),
                    ((89, "m"), 79),
                    ((90, "m"), 80),
                    ((91, "m"), 81),
                    ((92, "m"), 82),
                    ((93, "m"), 83),
                    ((94, "m"), 84),
                    ((95, "m"), 85),
                    ((96, "m"), 86),
                    ((97, "m"), 87),
                    ((98, "m"), 88),
                    ((99, "m"), 89),
                    ((100, "m"), 90),
                    ((101, "m"), 91),
                    ((102, "m"), 92),
                    ((103, "m"), 93),
                    ((104, "m"), 94),
                    ((105, "m"), 95),
                    ((106, "m"), 96),
                    ((107, "m"), 97),
                    ((108, "m"), 98),
                    ((109, "m"), 99),
                    ((110, "i"), 100),
                    ((110, "i"), 101),
                    ((110, "i"), 102),
                    ((110, "i"), 103),
                    ((110, "m"), 104),
                    ((111, "m"), 105),
                    ((112, "m"), 106),
                    ((113, "m"), 107),
                    ((114, "m"), 108),
                    ((115, "m"), 109),
                    ((116, "m"), 110),
                    ((117, "m"), 111),
                    ((118, "m"), 112),
                    ((119, "m"), 113),
                    ((120, "m"), 114),
                    ((121, "m"), 115),
                    ((122, "m"), 116),
                    ((123, "m"), 117),
                    ((124, "m"), 118),
                    ((125, "m"), 119),
                    ((126, "m"), 120),
                    ((127, "m"), 121),
                    ((128, "m"), 122),
                ]
            ],
            [
                {
                    "id": "cat_H",
                    "description": "",
                    "evalue": 4.8e-57,
                    "bitscore": 183.0,
                    "bias": 1.8,
                    "query_start": 0,
                    "query_end": 122,
                    "species": "cat",
                    "chain_type": "H",
                }
            ],
        )
    ]
    scheme = "imgt"
    allowed_chains = ["H", "K", "L"]
    assign_germline = True
    allowed_species = ["cat"]

    _numbered, _alignment_details, _hit_tables = number_sequences_from_alignment(
        sequences,
        _alignments,
        scheme=scheme,
        allow=allowed_chains,
        assign_germline=assign_germline,
        allowed_species=allowed_species,
    )


def test_numbering_seq():
    state_vector = [
        ((1, "m"), 0),
        ((2, "m"), 1),
        ((3, "m"), 2),
        ((4, "m"), 3),
        ((5, "m"), 4),
        ((6, "m"), 5),
        ((7, "m"), 6),
        ((8, "m"), 7),
        ((9, "m"), 8),
        ((10, "d"), None),
        ((11, "m"), 9),
        ((12, "m"), 10),
        ((13, "m"), 11),
        ((14, "m"), 12),
        ((15, "m"), 13),
        ((16, "m"), 14),
        ((17, "m"), 15),
        ((18, "m"), 16),
        ((19, "m"), 17),
        ((20, "m"), 18),
        ((21, "m"), 19),
        ((22, "m"), 20),
        ((23, "m"), 21),
        ((24, "m"), 22),
        ((25, "m"), 23),
        ((26, "m"), 24),
        ((27, "m"), 25),
        ((28, "m"), 26),
        ((29, "m"), 27),
        ((30, "m"), 28),
        ((31, "d"), None),
        ((32, "d"), None),
        ((33, "d"), None),
        ((34, "d"), None),
        ((35, "m"), 29),
        ((36, "m"), 30),
        ((37, "m"), 31),
        ((38, "m"), 32),
        ((39, "m"), 33),
        ((40, "m"), 34),
        ((41, "m"), 35),
        ((42, "m"), 36),
        ((43, "m"), 37),
        ((44, "m"), 38),
        ((45, "m"), 39),
        ((46, "m"), 40),
        ((47, "m"), 41),
        ((48, "m"), 42),
        ((49, "m"), 43),
        ((50, "m"), 44),
        ((51, "m"), 45),
        ((52, "m"), 46),
        ((53, "m"), 47),
        ((54, "m"), 48),
        ((55, "m"), 49),
        ((56, "m"), 50),
        ((57, "m"), 51),
        ((58, "m"), 52),
        ((59, "m"), 53),
        ((60, "d"), None),
        ((61, "d"), None),
        ((62, "d"), None),
        ((63, "m"), 54),
        ((64, "m"), 55),
        ((65, "m"), 56),
        ((66, "m"), 57),
        ((67, "m"), 58),
        ((68, "m"), 59),
        ((69, "m"), 60),
        ((70, "m"), 61),
        ((71, "m"), 62),
        ((72, "m"), 63),
        ((73, "d"), None),
        ((74, "m"), 64),
        ((75, "m"), 65),
        ((76, "m"), 66),
        ((77, "m"), 67),
        ((78, "m"), 68),
        ((79, "m"), 69),
        ((80, "m"), 70),
        ((81, "m"), 71),
        ((82, "m"), 72),
        ((83, "m"), 73),
        ((84, "m"), 74),
        ((85, "m"), 75),
        ((86, "m"), 76),
        ((87, "m"), 77),
        ((88, "m"), 78),
        ((89, "m"), 79),
        ((90, "m"), 80),
        ((91, "m"), 81),
        ((92, "m"), 82),
        ((93, "m"), 83),
        ((94, "m"), 84),
        ((95, "m"), 85),
        ((96, "m"), 86),
        ((97, "m"), 87),
        ((98, "m"), 88),
        ((99, "m"), 89),
        ((100, "m"), 90),
        ((101, "m"), 91),
        ((102, "m"), 92),
        ((103, "m"), 93),
        ((104, "m"), 94),
        ((105, "m"), 95),
        ((106, "m"), 96),
        ((107, "m"), 97),
        ((108, "m"), 98),
        ((109, "m"), 99),
        ((110, "i"), 100),
        ((110, "i"), 101),
        ((110, "i"), 102),
        ((110, "i"), 103),
        ((110, "m"), 104),
        ((111, "m"), 105),
        ((112, "m"), 106),
        ((113, "m"), 107),
        ((114, "m"), 108),
        ((115, "m"), 109),
        ((116, "m"), 110),
        ((117, "m"), 111),
        ((118, "m"), 112),
        ((119, "m"), 113),
        ((120, "m"), 114),
        ((121, "m"), 115),
        ((122, "m"), 116),
        ((123, "m"), 117),
        ((124, "m"), 118),
        ((125, "m"), 119),
        ((126, "m"), 120),
        ((127, "m"), 121),
        ((128, "m"), 122),
    ]
    seq = "DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF"
    scheme = "imgt"
    chain_type = "H"
    number_sequence_from_alignment(state_vector, seq, scheme=scheme, chain_type=chain_type)


def test_imgt():
    state_vector = [
        ((1, "m"), 0),
        ((2, "m"), 1),
        ((3, "m"), 2),
        ((4, "m"), 3),
        ((5, "m"), 4),
        ((6, "m"), 5),
        ((7, "m"), 6),
        ((8, "m"), 7),
        ((9, "m"), 8),
        ((10, "d"), None),
        ((11, "m"), 9),
        ((12, "m"), 10),
        ((13, "m"), 11),
        ((14, "m"), 12),
        ((15, "m"), 13),
        ((16, "m"), 14),
        ((17, "m"), 15),
        ((18, "m"), 16),
        ((19, "m"), 17),
        ((20, "m"), 18),
        ((21, "m"), 19),
        ((22, "m"), 20),
        ((23, "m"), 21),
        ((24, "m"), 22),
        ((25, "m"), 23),
        ((26, "m"), 24),
        ((27, "m"), 25),
        ((28, "m"), 26),
        ((29, "m"), 27),
        ((30, "m"), 28),
        ((31, "d"), None),
        ((32, "d"), None),
        ((33, "d"), None),
        ((34, "d"), None),
        ((35, "m"), 29),
        ((36, "m"), 30),
        ((37, "m"), 31),
        ((38, "m"), 32),
        ((39, "m"), 33),
        ((40, "m"), 34),
        ((41, "m"), 35),
        ((42, "m"), 36),
        ((43, "m"), 37),
        ((44, "m"), 38),
        ((45, "m"), 39),
        ((46, "m"), 40),
        ((47, "m"), 41),
        ((48, "m"), 42),
        ((49, "m"), 43),
        ((50, "m"), 44),
        ((51, "m"), 45),
        ((52, "m"), 46),
        ((53, "m"), 47),
        ((54, "m"), 48),
        ((55, "m"), 49),
        ((56, "m"), 50),
        ((57, "m"), 51),
        ((58, "m"), 52),
        ((59, "m"), 53),
        ((60, "d"), None),
        ((61, "d"), None),
        ((62, "d"), None),
        ((63, "m"), 54),
        ((64, "m"), 55),
        ((65, "m"), 56),
        ((66, "m"), 57),
        ((67, "m"), 58),
        ((68, "m"), 59),
        ((69, "m"), 60),
        ((70, "m"), 61),
        ((71, "m"), 62),
        ((72, "m"), 63),
        ((73, "d"), None),
        ((74, "m"), 64),
        ((75, "m"), 65),
        ((76, "m"), 66),
        ((77, "m"), 67),
        ((78, "m"), 68),
        ((79, "m"), 69),
        ((80, "m"), 70),
        ((81, "m"), 71),
        ((82, "m"), 72),
        ((83, "m"), 73),
        ((84, "m"), 74),
        ((85, "m"), 75),
        ((86, "m"), 76),
        ((87, "m"), 77),
        ((88, "m"), 78),
        ((89, "m"), 79),
        ((90, "m"), 80),
        ((91, "m"), 81),
        ((92, "m"), 82),
        ((93, "m"), 83),
        ((94, "m"), 84),
        ((95, "m"), 85),
        ((96, "m"), 86),
        ((97, "m"), 87),
        ((98, "m"), 88),
        ((99, "m"), 89),
        ((100, "m"), 90),
        ((101, "m"), 91),
        ((102, "m"), 92),
        ((103, "m"), 93),
        ((104, "m"), 94),
        ((105, "m"), 95),
        ((106, "m"), 96),
        ((107, "m"), 97),
        ((108, "m"), 98),
        ((109, "m"), 99),
        ((110, "i"), 100),
        ((110, "i"), 101),
        ((110, "i"), 102),
        ((110, "i"), 103),
        ((110, "m"), 104),
        ((111, "m"), 105),
        ((112, "m"), 106),
        ((113, "m"), 107),
        ((114, "m"), 108),
        ((115, "m"), 109),
        ((116, "m"), 110),
        ((117, "m"), 111),
        ((118, "m"), 112),
        ((119, "m"), 113),
        ((120, "m"), 114),
        ((121, "m"), 115),
        ((122, "m"), 116),
        ((123, "m"), 117),
        ((124, "m"), 118),
        ((125, "m"), 119),
        ((126, "m"), 120),
        ((127, "m"), 121),
        ((128, "m"), 122),
    ]
    sequence = "DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF"
    number_imgt(state_vector, sequence)


def test_anarci_multi_on():
    anarci_api = Anarci(run_multiproc=True, assign_germline=False)
    seq_records = [
        SeqRecord(
            Seq(
                "EVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSS"
            ),
            id=f"DupulimabH{i}",
        )
        for i in range(1000)
    ]
    _ = anarci_api.run_multiple(seq_records)
    
    
def test_anarci_multi_off():
    anarci_api = Anarci(run_multiproc=False, assign_germline=False)
    seq_records = [
        SeqRecord(
            Seq(
                "EVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSS"
            ),
            id=f"DupulimabH{i}",
        )
        for i in range(1000)
    ]
    _ = anarci_api.run_multiple(seq_records)
    
    
if __name__ == "__main__":
    for f in [test_anarci_multi_off]:
        profiler = Profiler()
        profiler.start()
        f()
        profiler.stop()
        profiler.print(show_all=True)
