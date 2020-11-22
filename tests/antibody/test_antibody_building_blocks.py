"""Unit tests for antibody objects lower level objects"""

import logging
import tempfile

import pytest
import pandas as pd
from Bio.Seq import Seq
from pkg_resources import resource_filename

from sadie.antibody import segment, exception, genetable

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def test_antibody_segment():
    """
    Test antibody segment base class, nt or AA
    """
    segment_aa = segment.AntibodySegment("DIQMTQSPASLSASLGETVSIECLAS")
    assert segment_aa.start == 1
    assert segment_aa.end == 26
    assert segment_aa.start_index == 0
    assert segment_aa.end_index == 25
    assert len(segment_aa) == len("DIQMTQSPASLSASLGETVSIECLAS")
    assert segment_aa.sequence == "DIQMTQSPASLSASLGETVSIECLAS"

    ##segmetn_nt
    segment_nt = segment.AntibodySegment("CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT")
    assert segment_nt.start == 1
    assert segment_nt.end == 75
    assert segment_nt.start_index == 0
    assert segment_nt.end_index == 74
    assert len(segment_nt) == len("CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT")
    assert segment_nt.sequence == "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT"
    assert not segment_nt.get_formatted_alignment()

    segment_nt.germline = "CAGGCGCAGCACGC"
    assert segment_nt.germline == "CAGGCGCAGCACGC"
    assert (
        segment_nt.get_formatted_alignment()
        == "germline    CAGGCGCAGCACG-------------------------------------------------------------C\ntarget      ....TT....TG.TGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT\n\n"
    )

    ##Explicit def
    segment_nt = segment.AntibodySegmentNT(
        "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCT"
    )
    assert segment_nt.aa == "QVQLVQSGAEVKKPGASVKVSCKAS"
    with pytest.raises(exception.BadNTSequenceError):
        segment.AntibodySegmentNT("ZTGCAGTC<TGGAGCT")

    segment_aa = segment.AntibodySegmentAA("DIQMTQSPASLSASLGETVSIECLAS")
    assert segment_aa == "DIQMTQSPASLSASLGETVSIECLAS"

    with pytest.raises(exception.BadAASequenceError):
        segment.AntibodySegmentAA("AVGQ<")


def test_framework_segment_nt():
    framework_1 = segment.FrameWork1NT("GACGTTCAGCTGGTGGAAAGTGGGGGTGACCTCTTGAAGCCGGGGGGCTCTCTCAGGTTAACATGCGTGGCTCCC")
    framework_2 = segment.FrameWork2NT("ATGAACTGGGTCTGTCAGGCACCTGGTAAGGGACTCCAGTGGGTAGCCTAT")
    framework_3 = segment.FrameWork3NT(
        "TACTACGCCGATAGCGTAAAGGGGAGATTTACCATCTCACGAGATAACGCCAAGAACACCTTGTACCTGCAGATGAACTCCTTGAAGGCCGAGGACACAGCCACACACTACTGT"
    )
    framework_4 = segment.FrameWork4NT("TGGGGCCACGGTACAATCGTAACCGTGAGCTCA")
    # fw1
    assert framework_1.start == 1
    assert framework_1.start_index == 0
    assert framework_1.end == 75
    assert framework_1.end_index == 74

    # fw2
    assert framework_2.start == 1
    assert framework_2.start_index == 0
    assert framework_2.end == 51
    assert framework_2.end_index == 50

    # fw3
    assert framework_3.start == 1
    assert framework_3.start_index == 0
    assert framework_3.end == 114
    assert framework_3.end_index == 113
    assert all(
        [
            isinstance(framework_1, segment.AntibodySegmentNT),
            isinstance(framework_2, segment.AntibodySegmentNT),
            isinstance(framework_3, segment.AntibodySegmentNT),
            isinstance(framework_4, segment.AntibodySegmentNT),
        ]
    )

    with pytest.raises(exception.BadNTSequenceError):
        framework_1 = segment.FrameWork1NT("DIQMTQSPASLSASLGET*VSIECLAS")
        framework_1 = segment.FrameWork1NT("ATCGN")
        framework_2 = segment.FrameWork2NT("LAxYQLKPGKSPQFLI")
        framework_3 = segment.FrameWork3NT("SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVxYC")
        framework_4 = segment.FrameWork4NT("FGKJ")


def test_framework_segment_aa():
    framework_1 = segment.FrameWork1AA("DIQMTQSPASLSASLGETVSIECLAS")
    framework_2 = segment.FrameWork2AA("LAWYQLKPGKSPQFLI")
    framework_3 = segment.FrameWork3AA("SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC")
    framework_4 = segment.FrameWork4AA("FGSGTKLKIK")
    assert framework_1.start == 1
    assert framework_1.start_index == 0
    assert framework_1.end == 26
    assert framework_1.end_index == 25

    # fw2
    assert framework_2.start == 1
    assert framework_2.start_index == 0
    assert framework_2.end == 16
    assert framework_2.end_index == 15

    # fw3
    assert framework_3.start == 1
    assert framework_3.start_index == 0
    assert framework_3.end == 36
    assert framework_3.end_index == 35

    # fw4
    assert framework_4.start == 1
    assert framework_4.start_index == 0
    assert framework_4.end == 10
    assert framework_4.end_index == 9
    assert all(
        [
            isinstance(framework_1, segment.AntibodySegmentAA),
            isinstance(framework_2, segment.AntibodySegmentAA),
            isinstance(framework_3, segment.AntibodySegmentAA),
            isinstance(framework_4, segment.AntibodySegmentAA),
        ]
    )

    with pytest.raises(exception.BadAASequenceError):
        framework_1 = segment.FrameWork1AA("DIQMTQSPASLSASLGET*VSIECLAS")
        framework_2 = segment.FrameWork2AA("LAxYQLKPGKSPQFLI")
        framework_3 = segment.FrameWork3AA("SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVxYC")
        framework_4 = segment.FrameWork4AA("FGKJ")


def test_cdr_segment_nt():
    cdr1 = segment.CDR1NT("GGCCTGTCCCTGACTTCTGGGTCT")
    cdr2 = segment.CDR2NT("ATTTACTCCAATGGTGGCACC")
    cdr3 = segment.CDR3NT("GCATCTATGTACTACTATGATGCCGACTACCTGCACTGGTATTTCGACTTC")
    # cdr1
    assert cdr1.start == 1
    assert cdr1.start_index == 0
    assert cdr1.end == 24
    assert cdr1.end_index == 23

    # cdr2
    assert cdr2.start == 1
    assert cdr2.start_index == 0
    assert cdr2.end == 21
    assert cdr2.end_index == 20

    # cdr3
    assert cdr3.start == 1
    assert cdr3.start_index == 0
    assert cdr3.end == 51
    assert cdr3.end_index == 50
    assert all(
        [
            isinstance(cdr1, segment.AntibodySegmentNT),
            isinstance(cdr2, segment.AntibodySegmentNT),
            isinstance(cdr3, segment.AntibodySegmentNT),
        ]
    )

    ###assign germlines
    cdr1.germline = "GGCCTGTCCCTGACTTCTGGGTCT"
    cdr2.germline = "GGCCTGTCCCTGACTTCTACT"
    cdr3.germline = ("CCGA", "GGCCTGTCC")
    assert cdr1.germline == "GGCCTGTCCCTGACTTCTGGGTCT"
    assert (
        cdr1.get_formatted_alignment()
        == "germline    GGCCTGTCCCTGACTTCTGGGTCT\ntarget      ........................\n\n"
    )
    assert cdr2.germline == "GGCCTGTCCCTGACTTCTACT"
    assert cdr2.get_formatted_alignment() == "germline    GGCCTGTCCCTGACTTCTACT\ntarget      ATTTAC...AATGG.GGC..C\n\n"

    assert cdr3.germline == "CCGA--------------------------------------GGCCTGTCC"
    assert (
        cdr3.get_formatted_alignment()
        == "germline    CCGA--------------------------------------GGCCTGTCC\ntarget      G.ATCTATGTACTACTATGATGCCGACTACCTGCACTGGTATTT.GAC.T.\n\n"
    )
    with pytest.raises(exception.BadNTSequenceError):
        cdr1 = segment.CDR1NT("DIQMTQSPASLSASLGET*VSIECLAS")
        cdr2 = segment.CDR2NT("AANC")
        cdr3 = segment.CDR3NT("CCACD")


def test_cdr_segment_aa():
    cdr1 = segment.CDR1AA("GLSLTSNS")
    cdr2 = segment.CDR2AA("IWSNGGT")
    cdr3 = segment.CDR3AA("ASIYYYDADYLHWYFDF")
    assert cdr1.start == 1
    assert cdr1.start_index == 0
    assert cdr1.end == 8
    assert cdr1.end_index == 7

    # fw2
    assert cdr2.start == 1
    assert cdr2.start_index == 0
    assert cdr2.end == 7
    assert cdr2.end_index == 6

    # fw3
    assert cdr3.start == 1
    assert cdr3.start_index == 0
    assert cdr3.end == 17
    assert cdr3.end_index == 16
    assert all(
        [
            isinstance(cdr3, segment.AntibodySegmentAA),
            isinstance(cdr2, segment.AntibodySegmentAA),
            isinstance(cdr3, segment.AntibodySegmentAA),
        ]
    )

    cdr1.germline = "GLSLGSNS"
    assert cdr1.get_formatted_alignment() == "germline    GLSLGSNS\ntarget      ....T...\n\n"
    cdr2.germline = "IWSNGAT"
    assert cdr2.get_formatted_alignment() == "germline    IWSNGAT\ntarget      .....G.\n\n"
    cdr3.germline = ("AT", "YV")
    assert cdr3.germline == "AT-------------YV"
    assert cdr3.get_formatted_alignment() == "germline    AT-------------YV\ntarget      .SIYYYDADYLHWYFDF\n\n"
    with pytest.raises(exception.BadAASequenceError):
        ##Z should cause error
        cdr1 = segment.CDR1AA("DIQMTQSPASLSASLGEZ")

    with pytest.warns(exception.BadAASequenceWarning):
        ##Z should cause error
        cdr1 = segment.CDR1AA("DIQMTQSPASLSASLGE*")
        cdr2 = segment.CDR2AA("DIQMTQSPASLSASLGEX")
