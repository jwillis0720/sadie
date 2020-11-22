"""Unit tests for antibody objects"""

import logging
import tempfile

import pytest
import pandas as pd
from Bio.Seq import Seq
from pkg_resources import resource_filename

from sadie.antibody import exception, segment
from sadie import antibody

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def test_antibody_chain_aa():
    heavy_chain_aa = antibody.AntibodyChainAA(
        fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
        cdr1_aa="GLSLTSNS",
        fwr2_aa="VSWIRQPPGKGLEWMGV",
        cdr2_aa="IWSNGGT",
        fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
        cdr3_aa="ASIYYYDADYLHWYFDF",
        fwr4_aa="WGPGTMVTVSS",
        v_gene="IGHV2-47",
        j_gene="IGHJ1",
        species="rat",
    )
    heavy_chain_aa.get_json()
    assert all(
        [
            isinstance(heavy_chain_aa, antibody.AntibodyChainAA),
            isinstance(heavy_chain_aa.cdr1_aa, segment.CDR1AA),
            isinstance(heavy_chain_aa.cdr2_aa, segment.CDR2AA),
            isinstance(heavy_chain_aa.cdr3_aa, segment.CDR3AA),
            isinstance(heavy_chain_aa.fwr1_aa, segment.FrameWork1AA),
            isinstance(heavy_chain_aa.fwr2_aa, segment.FrameWork2AA),
            isinstance(heavy_chain_aa.fwr3_aa, segment.FrameWork3AA),
            isinstance(heavy_chain_aa.fwr4_aa, segment.FrameWork4AA),
        ]
    )

    assert heavy_chain_aa.cdr1_aa == "GLSLTSNS"
    assert heavy_chain_aa.cdr2_aa == "IWSNGGT"
    assert heavy_chain_aa.cdr3_aa == "ASIYYYDADYLHWYFDF"
    assert heavy_chain_aa.fwr1_aa == "QVQLKESGPGLVQPSQTLSLTCTVS"
    assert heavy_chain_aa.fwr2_aa == "VSWIRQPPGKGLEWMGV"
    assert heavy_chain_aa.fwr3_aa == "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC"
    assert heavy_chain_aa.fwr1_aa_germline == "QVQLKESGPGLVQPSQTLSLTCTVS"
    assert heavy_chain_aa.cdr1_aa_germline == "GLSLTSNS"
    assert heavy_chain_aa.fwr2_aa_germline == "VSWIRQPPGKGLEWMGV"
    assert heavy_chain_aa.cdr2_aa_germline == "IWSNGGT"
    assert heavy_chain_aa.fwr3_aa_germline == "DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC"
    assert heavy_chain_aa.cdr3_aa_germline_v == "AR"
    assert heavy_chain_aa.cdr3_aa_germline_j == "YYWYFDF"

    segmented_aa = "QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS"
    assert heavy_chain_aa.get_segmented_vdj_aa() == segmented_aa
    assert heavy_chain_aa.__str__() == segmented_aa.replace(" ", "")
    segmented_alignment = "IGHV2-47|IGHJ1  QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC \nantibodychain   ......................... ........ ................. ....... ......E.....N................P........ \n\nIGHV2-47|IGHJ1  AR--------YYWYFDF WGPGTMVTVSS\nantibodychain   .SIYYYDADYLH..... ...........\n\n"
    assert heavy_chain_aa.get_segmented_alignment_aa() == segmented_alignment
    assert (
        heavy_chain_aa.__repr__()
        == "antibodychain\nFrameWork1AA 1-25:QVQLKESGPGLVQPSQTLSLTCTVS\nCDR1AA 26-33:GLSLTSNS\nFrameWork2AA 34-50:VSWIRQPPGKGLEWMGV\nCDR2AA 51-57:IWSNGGT\nFrameWork3AA 58-95:DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC\nCDR3AA 96-112:ASIYYYDADYLHWYFDF\nFrameWork4AA 113-123:WGPGTMVTVSS"
    )


def test_antibody_chain_nt():
    chain_nt = antibody.AntibodyChainNT(
        fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
        cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
        fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
        cdr2_nt="ATTTGGTCCAACGGTGGCACC",
        fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
        cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
        fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
        v_gene="IGHV2-47",
        j_gene="IGHJ3",
        species="rat",
    )
    assert all(
        [
            isinstance(chain_nt, antibody.AntibodyChainNT),
            isinstance(chain_nt.cdr1_nt, segment.CDR1NT),
            isinstance(chain_nt.cdr2_nt, segment.CDR2NT),
            isinstance(chain_nt.cdr3_nt, segment.CDR3NT),
            isinstance(chain_nt.fwr1_nt, segment.FrameWork1NT),
            isinstance(chain_nt.fwr2_nt, segment.FrameWork2NT),
            isinstance(chain_nt.fwr3_nt, segment.FrameWork3NT),
            isinstance(chain_nt.fwr4_nt, segment.FrameWork4NT),
        ]
    )

    ##Get nucletodies
    assert chain_nt.fwr1_nt == "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA"
    assert chain_nt.cdr1_nt == "GGCCTATCGCTCACCAGCAACTCC"
    assert chain_nt.fwr2_nt == "GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG"
    assert chain_nt.cdr2_nt == "ATTTGGTCCAACGGTGGCACC"
    assert (
        chain_nt.fwr3_nt
        == "GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC"
    )
    assert chain_nt.cdr3_nt == "GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC"
    assert chain_nt.fwr4_nt == "TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC"

    ##Get nucletodie germline
    assert chain_nt.fwr1_nt_germline == "CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT"
    assert chain_nt.cdr1_nt_germline == "GGGTTATCATTAACCAGCAATAGT"

    assert chain_nt.fwr2_nt_germline == "GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA"
    assert chain_nt.cdr2_nt_germline == "ATATGGAGTAATGGAGGCACA"
    assert (
        chain_nt.fwr3_nt_germline
        == "GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT"
    )
    assert chain_nt.fwr4_nt_germline == "TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG"

    ##Get amino acid sequences
    assert chain_nt.fwr1_aa == "QVQLKESGPGLVQPSQTLSLTCTVS"
    assert chain_nt.cdr1_aa == "GLSLTSNS"
    assert chain_nt.fwr2_aa == "VSWIRQPPGKGLEWMGV"
    assert chain_nt.cdr2_aa == "IWSNGGT"
    assert chain_nt.fwr3_aa == "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC"
    assert chain_nt.cdr3_aa == "ASIYYYDADYLHWYFDF"
    assert chain_nt.fwr4_aa == "WGPGTMVTVSS"

    assert (
        chain_nt.get_segmented_vdj_nt()
        == "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA GGCCTATCGCTCACCAGCAACTCC GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG ATTTGGTCCAACGGTGGCACC GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC"
    )
    assert (
        chain_nt.get_segmented_vdj_aa()
        == "QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS"
    )

    assert (
        chain_nt.get_segmented_alignment_nt()
        == "IGHV2-47|IGHJ3  CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT\nantibodychain   ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC\n\nIGHV2-47|IGHJ3   GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCC\nantibodychain    ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG..\n\nIGHV2-47|IGHJ3  GACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--\nantibodychain   .CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTT\n\nIGHV2-47|IGHJ3  ------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG\nantibodychain   ATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C\n\n"
    )

    assert (
        chain_nt.get_segmented_alignment_aa()
        == "IGHV2-47|IGHJ3  QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC \nantibodychain   ......................... ........ ................. ....... ......E.....N................P........ \n\nIGHV2-47|IGHJ3  AR----------NWFAY WGQGTLVTVSS\nantibodychain   .SIYYYDADYLHWY.DF ..P..M.....\n\n"
    )


##Test heavy chain, kappa and lambda chain AA
def test_antibody_chain_heavy():
    heavy_chain_aa = antibody.HeavyChainAA(
        fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
        cdr1_aa="GLSLTSNS",
        fwr2_aa="VSWIRQPPGKGLEWMGV",
        cdr2_aa="IWSNGGT",
        fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
        cdr3_aa="ASIYYYDADYLHWYFDF",
        fwr4_aa="WGPGTMVTVSS",
        v_gene="IGHV2-47",
        j_gene="IGHJ1",
        species="rat",
    )
    assert all(
        [
            isinstance(heavy_chain_aa, antibody.HeavyChainAA),
            isinstance(heavy_chain_aa.cdr1_aa, segment.CDR1AA),
            isinstance(heavy_chain_aa.cdr2_aa, segment.CDR2AA),
            isinstance(heavy_chain_aa.cdr3_aa, segment.CDR3AA),
            isinstance(heavy_chain_aa.fwr1_aa, segment.FrameWork1AA),
            isinstance(heavy_chain_aa.fwr2_aa, segment.FrameWork2AA),
            isinstance(heavy_chain_aa.fwr3_aa, segment.FrameWork3AA),
            isinstance(heavy_chain_aa.fwr4_aa, segment.FrameWork4AA),
        ]
    )

    assert heavy_chain_aa.cdr1_aa == "GLSLTSNS"
    assert heavy_chain_aa.cdr2_aa == "IWSNGGT"
    assert heavy_chain_aa.cdr3_aa == "ASIYYYDADYLHWYFDF"
    assert heavy_chain_aa.fwr1_aa == "QVQLKESGPGLVQPSQTLSLTCTVS"
    assert heavy_chain_aa.fwr2_aa == "VSWIRQPPGKGLEWMGV"
    assert heavy_chain_aa.fwr3_aa == "DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC"
    assert heavy_chain_aa.fwr1_aa_germline == "QVQLKESGPGLVQPSQTLSLTCTVS"
    assert heavy_chain_aa.cdr1_aa_germline == "GLSLTSNS"
    assert heavy_chain_aa.fwr2_aa_germline == "VSWIRQPPGKGLEWMGV"
    assert heavy_chain_aa.cdr2_aa_germline == "IWSNGGT"
    assert heavy_chain_aa.fwr3_aa_germline == "DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC"
    assert heavy_chain_aa.cdr3_aa_germline_v == "AR"
    assert heavy_chain_aa.cdr3_aa_germline_j == "YYWYFDF"
    segmented_aa = "QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC ASIYYYDADYLHWYFDF WGPGTMVTVSS"
    assert heavy_chain_aa.get_segmented_vdj_aa() == segmented_aa
    assert heavy_chain_aa.__str__() == segmented_aa.replace(" ", "")
    segmented_alignment = "IGHV2-47|IGHJ1  QVQLKESGPGLVQPSQTLSLTCTVS GLSLTSNS VSWIRQPPGKGLEWMGV IWSNGGT DYNSAIKSRLSISRDTSKSQVFLKMNSLQTEDTAMYFC \nheavy_chain     ......................... ........ ................. ....... ......E.....N................P........ \n\nIGHV2-47|IGHJ1  AR--------YYWYFDF WGPGTMVTVSS\nheavy_chain     .SIYYYDADYLH..... ...........\n\n"
    assert heavy_chain_aa.get_segmented_alignment_aa() == segmented_alignment
    assert (
        heavy_chain_aa.__repr__()
        == "heavy_chain\nFrameWork1AA 1-25:QVQLKESGPGLVQPSQTLSLTCTVS\nCDR1AA 26-33:GLSLTSNS\nFrameWork2AA 34-50:VSWIRQPPGKGLEWMGV\nCDR2AA 51-57:IWSNGGT\nFrameWork3AA 58-95:DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC\nCDR3AA 96-112:ASIYYYDADYLHWYFDF\nFrameWork4AA 113-123:WGPGTMVTVSS"
    )
    assert heavy_chain_aa.locus == "IGH"


def test_antibody_chain_kappa():
    kappa_chain_aa = antibody.KappaChainAA(
        fwr1_aa="DIQMTQSPASLSASLGETVSIECLAS",
        cdr1_aa="EGISNS",
        fwr2_aa="LAWYQLKPGKSPQFLI",
        cdr2_aa="YATS",
        fwr3_aa="SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC",
        cdr3_aa="QQGYKFPLT",
        fwr4_aa="FGSGTKLKIK",
        v_gene="IGKV12S11",
        j_gene="IGKJ5*01",
        species="rat",
    )

    assert all(
        [
            isinstance(kappa_chain_aa, antibody.KappaChainAA),
            isinstance(kappa_chain_aa.cdr1_aa, segment.CDR1AA),
            isinstance(kappa_chain_aa.cdr2_aa, segment.CDR2AA),
            isinstance(kappa_chain_aa.cdr3_aa, segment.CDR3AA),
            isinstance(kappa_chain_aa.fwr1_aa, segment.FrameWork1AA),
            isinstance(kappa_chain_aa.fwr2_aa, segment.FrameWork2AA),
            isinstance(kappa_chain_aa.fwr3_aa, segment.FrameWork3AA),
            isinstance(kappa_chain_aa.fwr4_aa, segment.FrameWork4AA),
        ]
    )
    assert kappa_chain_aa.cdr1_aa == "EGISNS"
    assert kappa_chain_aa.cdr2_aa == "YATS"
    assert kappa_chain_aa.cdr3_aa == "QQGYKFPLT"
    assert kappa_chain_aa.fwr1_aa == "DIQMTQSPASLSASLGETVSIECLAS"
    assert kappa_chain_aa.fwr2_aa == "LAWYQLKPGKSPQFLI"
    assert kappa_chain_aa.fwr3_aa == "SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC"
    assert kappa_chain_aa.fwr4_aa == "FGSGTKLKIK"

    assert kappa_chain_aa.fwr1_aa_germline == "DIQMTQSPHSLSASLGETVSIECLAS"
    assert kappa_chain_aa.cdr1_aa_germline == "EGISNY"
    assert kappa_chain_aa.fwr2_aa_germline == "LAWYQQKPGKSPQLLIY"
    assert kappa_chain_aa.cdr2_aa_germline == "YAS"
    assert kappa_chain_aa.fwr3_aa_germline == "SLQDGVPSRFSGSGSGTQYSLKISNMQPEDEGVYYC"
    assert kappa_chain_aa.fwr4_aa_germline == "FGSGTKLEIK"

    assert (
        kappa_chain_aa.get_segmented_vdj_aa()
        == "DIQMTQSPASLSASLGETVSIECLAS EGISNS LAWYQLKPGKSPQFLI YATS SLQDGVPSRFSGSGSGTQYSLKISGMQPEDEGVYYC QQGYKFPLT FGSGTKLKIK"
    )
    segmented_align = "kappa_chain         DIQMTQSPHSLSASLGETVSIECLAS EGISNY LAWYQQKPGKSPQLLIY YA-S SLQDGVPSRFSGSGSGTQYSLKISNMQPEDEGVYYC QQGYKY\nIGKV12S11|IGKJ5*01  ........A................. .....S .....L.......F..- ..T. ........................G........... .....F\n\nkappa_chain         PLT FGSGTKLEIK\nIGKV12S11|IGKJ5*01  ... .......K..\n\n"
    segmented_align = "IGKV12S11|IGKJ5*01  DIQMTQSPHSLSASLGETVSIECLAS EGISNY LAWYQQKPGKSPQLLIY YA-S SLQDGVPSRFSGSGSGTQYSLKISNMQPEDEGVYYC QQGYKY\nkappa_chain         ........A................. .....S .....L.......F..- ..T. ........................G........... .....F\n\nIGKV12S11|IGKJ5*01  PLT FGSGTKLEIK\nkappa_chain         ... .......K..\n\n"
    assert kappa_chain_aa.get_segmented_alignment_aa() == segmented_align


def test_antibody_chain_lambda():
    lambda_chain_aa = antibody.LambdaChainAA(
        name="fezakinumab",
        fwr1_aa="QAVLTQPPSVSGAPGQRVTISCTGS",
        cdr1_aa="SSNIGAGYG",
        fwr2_aa="VHWYQQLPGTAPKLLIY",
        cdr2_aa="GDS",
        fwr3_aa="NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC",
        cdr3_aa="QSYDNSLSGYV",
        fwr4_aa="FGGGTQLTVL",
        v_gene="IGLV1-40*01",
        j_gene="IGLJ7*01",
        species="human",
    )

    assert all(
        [
            isinstance(lambda_chain_aa, antibody.LambdaChainAA),
            isinstance(lambda_chain_aa.cdr1_aa, segment.CDR1AA),
            isinstance(lambda_chain_aa.cdr2_aa, segment.CDR2AA),
            isinstance(lambda_chain_aa.cdr3_aa, segment.CDR3AA),
            isinstance(lambda_chain_aa.fwr1_aa, segment.FrameWork1AA),
            isinstance(lambda_chain_aa.fwr2_aa, segment.FrameWork2AA),
            isinstance(lambda_chain_aa.fwr3_aa, segment.FrameWork3AA),
            isinstance(lambda_chain_aa.fwr4_aa, segment.FrameWork4AA),
        ]
    )
    assert lambda_chain_aa.cdr1_aa == "SSNIGAGYG"
    assert lambda_chain_aa.cdr2_aa == "GDS"
    assert lambda_chain_aa.cdr3_aa == "QSYDNSLSGYV"
    assert lambda_chain_aa.fwr1_aa == "QAVLTQPPSVSGAPGQRVTISCTGS"
    assert lambda_chain_aa.fwr2_aa == "VHWYQQLPGTAPKLLIY"
    assert lambda_chain_aa.fwr3_aa == "NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC"
    assert lambda_chain_aa.fwr4_aa == "FGGGTQLTVL"

    assert lambda_chain_aa.fwr1_aa_germline == "QSVLTQPPSVSGAPGQRVTISCTGS"
    assert lambda_chain_aa.cdr1_aa_germline == "SSNIGAGYD"
    assert lambda_chain_aa.fwr2_aa_germline == "VHWYQQLPGTAPKLLIY"
    assert lambda_chain_aa.cdr2_aa_germline == "GNS"
    assert lambda_chain_aa.fwr3_aa_germline == "NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC"
    assert lambda_chain_aa.fwr4_aa_germline == "FGGGTQLTVL"

    assert (
        lambda_chain_aa.get_segmented_vdj_aa()
        == "QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSGYV FGGGTQLTVL"
    )
    segmented_aling = "IGLV1-40*01|IGLJ7*01  QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDS\nfezakinumab           .A....................... ........G ................. .D. .................................... ....N\n\nIGLV1-40*01|IGLJ7*01  SLSGAV FGGGTQLTVL\nfezakinumab           ....Y. ..........\n\n"
    assert lambda_chain_aa.get_segmented_alignment_aa() == segmented_aling


##Test nt objects
def test_antibody_chain_heavy_nt():
    heavy_chain_nt = antibody.HeavyChainNT(
        fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
        cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
        fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
        cdr2_nt="ATTTGGTCCAACGGTGGCACC",
        fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
        cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
        fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
        v_gene="IGHV2-47",
        j_gene="IGHJ3",
        species="rat",
    )
    assert all(
        [
            isinstance(heavy_chain_nt, antibody.HeavyChainNT),
            isinstance(heavy_chain_nt.cdr1_nt, segment.CDR1NT),
            isinstance(heavy_chain_nt.cdr2_nt, segment.CDR2NT),
            isinstance(heavy_chain_nt.cdr3_nt, segment.CDR3NT),
            isinstance(heavy_chain_nt.fwr1_nt, segment.FrameWork1NT),
            isinstance(heavy_chain_nt.fwr2_nt, segment.FrameWork2NT),
            isinstance(heavy_chain_nt.fwr3_nt, segment.FrameWork3NT),
            isinstance(heavy_chain_nt.fwr4_nt, segment.FrameWork4NT),
        ]
    )

    assert heavy_chain_nt.cdr1_nt == "GGCCTATCGCTCACCAGCAACTCC"
    assert heavy_chain_nt.cdr2_nt == "ATTTGGTCCAACGGTGGCACC"
    assert heavy_chain_nt.cdr3_nt == "GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC"
    assert heavy_chain_nt.fwr1_nt == "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA"
    assert heavy_chain_nt.fwr2_nt == "GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG"
    assert (
        heavy_chain_nt.fwr3_nt
        == "GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC"
    )
    assert (
        heavy_chain_nt.fwr1_nt_germline == "CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT"
    )
    assert heavy_chain_nt.cdr1_nt_germline == "GGGTTATCATTAACCAGCAATAGT"
    assert heavy_chain_nt.fwr2_nt_germline == "GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA"
    assert heavy_chain_nt.cdr2_nt_germline == "ATATGGAGTAATGGAGGCACA"
    assert (
        heavy_chain_nt.fwr3_nt_germline
        == "GATTATAATTCAGCTATCAAATCCCGACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT"
    )
    assert heavy_chain_nt.cdr3_nt_germline_v == "GCCAGAAA"
    segmented_nt = "CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA GGCCTATCGCTCACCAGCAACTCC GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG ATTTGGTCCAACGGTGGCACC GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC"
    assert heavy_chain_nt.get_segmented_vdj_nt() == segmented_nt
    assert heavy_chain_nt.__str__() == segmented_nt.replace(" ", "")
    segmented_alignment = "IGHV2-47|IGHJ3  CAAGTGCAACTAAAGGAGTCAGGACCTGGTCTGGTACAGCCATCACAGACCCTGTCTCTCACCTGCACTGTCTCT GGGTTATCATTAACCAGCAATAGT\nheavy_chain     ..G.....G..G......AGC..C......T....G...........A..T..T.....G..A.....C..G..A ..CC....GC.C........CTCC\n\nIGHV2-47|IGHJ3   GTAAGCTGGATTCGGCAGCCTCCAGGAAAGGGTCTGGAGTGGATGGGAGTA ATATGGAGTAATGGAGGCACA GATTATAATTCAGCTATCAAATCCC\nheavy_chain      ..C........A..T.....G.....C..A.................T..G ..T...TCC..C..T.....C ..C..C..C..C......G.GAG..\n\nIGHV2-47|IGHJ3  GACTGAGCATCAGCAGGGACACCTCGAAGAGCCAAGTTTTCTTAAAGATGAACAGTCTGCAAACTGAAGACACAGCCATGTACTTCTGT GCCAGAAA--\nheavy_chain     .CT..TCT....A.C.C........T..ATC...G........G.........TCG..T...C....G..T..G..T........T..C ...TCC.TTT\n\nIGHV2-47|IGHJ3  ------------------------ACAATTGGTTTGCTTAC TGGGGCCAAGGCACTCTGGTCACTGTCTCTTCAG\nheavy_chain     ATTACTATGACGCTGACTACCTCC..TGG.AC..C.A..T. .......CC......A....G..C..GAGC..-C\n\n"
    assert heavy_chain_nt.get_segmented_alignment_nt() == segmented_alignment
    assert heavy_chain_nt.locus == "IGH"


def test_antibody_chain_kappa_nt():
    kappa_chain_nt = antibody.KappaChainNT(
        fwr1_nt="GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
        cdr1_nt="CACTCAATCTACCGGAAT",
        fwr2_nt="CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
        cdr2_nt="AACGCCAAC",
        fwr3_nt="TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
        cdr3_nt="CAACAGTACTATCGAGGATGGACG",
        fwr4_nt="TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
        v_gene="IGKV22S4",
        j_gene="IGKJ1",
        species="rat",
    )
    assert all(
        [
            isinstance(kappa_chain_nt, antibody.KappaChainNT),
            isinstance(kappa_chain_nt.cdr1_nt, segment.CDR1NT),
            isinstance(kappa_chain_nt.cdr2_nt, segment.CDR2NT),
            isinstance(kappa_chain_nt.cdr3_nt, segment.CDR3NT),
            isinstance(kappa_chain_nt.fwr1_nt, segment.FrameWork1NT),
            isinstance(kappa_chain_nt.fwr2_nt, segment.FrameWork2NT),
            isinstance(kappa_chain_nt.fwr3_nt, segment.FrameWork3NT),
            isinstance(kappa_chain_nt.fwr4_nt, segment.FrameWork4NT),
        ]
    )

    assert kappa_chain_nt.cdr1_nt == "CACTCAATCTACCGGAAT"
    assert kappa_chain_nt.cdr2_nt == "AACGCCAAC"
    assert kappa_chain_nt.cdr3_nt == "CAACAGTACTATCGAGGATGGACG"
    assert kappa_chain_nt.fwr1_nt == "GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC"
    assert kappa_chain_nt.fwr2_nt == "CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC"
    assert (
        kappa_chain_nt.fwr3_nt
        == "TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC"
    )
    assert (
        kappa_chain_nt.fwr1_nt_germline
        == "GACATCCAGATGACCCAGTCTCCTTCATTCCTGTCTGCATCTGTGGGAGACAGAGTCACTATCAACTGCAAAGCAAGT"
    )
    assert kappa_chain_nt.cdr1_nt_germline == "CAGAATATTAACAGGTAC"
    assert kappa_chain_nt.fwr2_nt_germline == "TTAAACTGGTACCAGCAAAAGCTTGGAGAAGCTCCCAAACTCCTGATATAT"
    assert kappa_chain_nt.cdr2_nt_germline == "AATGCAAAC"
    assert (
        kappa_chain_nt.fwr3_nt_germline
        == "AGTTTGCAAACGGGCATCCCATCAAGGTTCAGTGGCAGTGGATCTGGTACTGATTTCACACTCACCATCAGCAGCCTGCAGCCTGAAGATGTTGCCACATATTTCTGC"
    )
    assert kappa_chain_nt.cdr3_nt_germline_v == "TTGCAGCATAATAGTTGGCCG"
    segmented_nt = "GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC CACTCAATCTACCGGAAT CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC AACGCCAAC TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC CAACAGTACTATCGAGGATGGACG TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG"
    assert kappa_chain_nt.get_segmented_vdj_nt() == segmented_nt
    assert kappa_chain_nt.__str__() == segmented_nt.replace(" ", "")
    segmented_alignment = "IGKV22S4|IGKJ1  GACATCCAGATGACCCAGTCTCCTTCATTCCTGTCTGCATCTGTGGGAGACAGAGTCACTATCAACTGCAAAGCAAGT CAGAATATTAACAGGTAC TT\nkappa_chain     ........A.....A..T..G........G...AG...G........T...C.C....G.C.G........G..CTCC ..CTCA..CT..C..A.T C.\n\nIGKV22S4|IGKJ1  AAACTGGTACCAGCAAAAGCTTGGAGAAGCTCCCAAACTCCTGATATAT AATGCAAAC AGTTTGCAAACGGGCATCCCATCAAGGTTCAGTGGCAGTG\nkappa_chain     GGC.........A..G..A..C..T..G.....A.....A..C..C..C ..C..C... TC.C....G..A..A.....G..T..A..T..C..ATCC.\n\nIGKV22S4|IGKJ1  GATCTGGTACTGATTTCACACTCACCATCAGCAGCCTGCAGCCTGAAGATGTTGCCACATATTTCTGC TTGCAGCATAATAGTTGGCCGGTGGACG T-\nkappa_chain     .C..C.....C..C.....C..G.....T...TC.........C..G.....G..G..C......... CAA...T.CT.----.C.AG.A...... .G\n\nIGKV22S4|IGKJ1  -TCGGTGGAGGCACCAAGCTGGAATTGAAAC\nkappa_chain     GA..T.C.GT.G.GGT.CAAA.CTGGAGCTG\n\n"
    assert kappa_chain_nt.get_segmented_alignment_nt() == segmented_alignment
    assert kappa_chain_nt.locus == "IGK"

    assert kappa_chain_nt.cdr1_aa == "HSIYRN"
    assert kappa_chain_nt.cdr2_aa == "NAN"
    assert kappa_chain_nt.cdr3_aa == "QQYYRGWT"
    assert kappa_chain_nt.fwr1_aa == "DIQMTHSPSLLSASVGDRVSLNCKAS"
    assert kappa_chain_nt.fwr2_aa == "LAWYQQKLGEAPKLLIY"
    assert kappa_chain_nt.fwr3_aa == "SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC"
    assert kappa_chain_nt.fwr4_aa == "WTFGGGTKLEL"

    assert kappa_chain_nt.fwr1_aa_germline == "DIQMTQSPSFLSASVGDRVTINCKAS"
    assert kappa_chain_nt.cdr1_aa_germline == "QNINRY"
    assert kappa_chain_nt.fwr2_aa_germline == "LNWYQQKLGEAPKLLIY"
    assert kappa_chain_nt.cdr2_aa_germline == "NAN"
    assert kappa_chain_nt.fwr3_aa_germline == "SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC"
    assert kappa_chain_nt.fwr4_aa_germline == "FGGGTKLELK"

    assert (
        kappa_chain_nt.get_segmented_vdj_aa()
        == "DIQMTHSPSLLSASVGDRVSLNCKAS HSIYRN LAWYQQKLGEAPKLLIY NAN SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC QQYYRGWT WTFGGGTKLEL"
    )
    segmented_alignment = "IGKV22S4|IGKJ1  DIQMTQSPSFLSASVGDRVTINCKAS QNINRY LNWYQQKLGEAPKLLIY NAN SLQTGIPSRFSGSGSGTDFTLTISSLQPEDVATYFC LQHNSWP\nkappa_chain     .....H...L.........SL..... HS.Y.N .A............... ... .................................... -.QYYRG\n\nIGKV22S4|IGKJ1  WT FGGGTKLELK-\nkappa_chain     .. WTF.GGTK.EL\n\n"
    assert kappa_chain_nt.get_segmented_alignment_aa() == segmented_alignment
    assert (
        kappa_chain_nt.fwr2_nt.get_formatted_alignment()
        == "germline    TTAAACTGGTACCAGCAAAAGCTTGGAGAAGCTCCCAAACTCCTGATATAT\ntarget      C.GGC.........A..G..A..C..T..G.....A.....A..C..C..C\n\n"
    )


def test_antibody_chain_lamba_nt():
    lambda_chain_nt = antibody.LambdaChainNT(
        name="fezakinumab",
        fwr1_nt="CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC",
        cdr1_nt="TCTTCCAACATCGGCGCCGGTTACGGC",
        fwr2_nt="GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC",
        cdr2_nt="GGCGACAGC",
        fwr3_nt="AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC",
        cdr3_nt="CAGTCTTATGACAATTCCTTGAGTGGC",
        fwr4_nt="TTCGGGGGAGGGACCCAGTTGACTGTTCTT",
        v_gene="IGLV1-40*01",
        j_gene="IGLJ7*01",
        species="human",
    )
    assert all(
        [
            isinstance(lambda_chain_nt, antibody.LambdaChainNT),
            isinstance(lambda_chain_nt.cdr1_nt, segment.CDR1NT),
            isinstance(lambda_chain_nt.cdr2_nt, segment.CDR2NT),
            isinstance(lambda_chain_nt.cdr3_nt, segment.CDR3NT),
            isinstance(lambda_chain_nt.fwr1_nt, segment.FrameWork1NT),
            isinstance(lambda_chain_nt.fwr2_nt, segment.FrameWork2NT),
            isinstance(lambda_chain_nt.fwr3_nt, segment.FrameWork3NT),
            isinstance(lambda_chain_nt.fwr4_nt, segment.FrameWork4NT),
        ]
    )
    assert lambda_chain_nt.cdr1_nt == "TCTTCCAACATCGGCGCCGGTTACGGC"
    assert lambda_chain_nt.cdr2_nt == "GGCGACAGC"
    assert lambda_chain_nt.cdr3_nt == "CAGTCTTATGACAATTCCTTGAGTGGC"
    assert lambda_chain_nt.fwr1_nt == "CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC"
    assert lambda_chain_nt.fwr2_nt == "GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC"
    assert (
        lambda_chain_nt.fwr3_nt
        == "AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC"
    )
    assert (
        lambda_chain_nt.fwr1_nt_germline
        == "CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGC"
    )
    assert lambda_chain_nt.cdr1_nt_germline == "AGCTCCAACATCGGGGCAGGTTATGAT"
    assert lambda_chain_nt.fwr2_nt_germline == "GTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTAT"
    assert lambda_chain_nt.cdr2_nt_germline == "GGTAACAGC"
    assert (
        lambda_chain_nt.fwr3_nt_germline
        == "AATCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGC"
    )
    assert lambda_chain_nt.cdr3_nt_germline_v == "CAGTCCTATGACAGCAGCCTGAGTGGTTC"
    segmented_nt = "CAGGCGGTGCTCACCCAGCCACCTAGTGTGAGCGGTGCACCTGGGCAGCGTGTGACCATCTCTTGCACTGGGTCC TCTTCCAACATCGGCGCCGGTTACGGC GTGCACTGGTACCAACAGCTTCCGGGCACCGCCCCCAAGCTGCTCATCTAC GGCGACAGC AATCGTCCATCAGGGGTTCCGGATCGCTTTAGCGGGTCTAAGTCAGGGACCTCAGCCTCCCTGGCGATCACTGGGCTGCAGGCGGAGGACGAGGCAGACTATTACTGC CAGTCTTATGACAATTCCTTGAGTGGC TTCGGGGGAGGGACCCAGTTGACTGTTCTT"
    assert lambda_chain_nt.get_segmented_vdj_nt() == segmented_nt
    assert lambda_chain_nt.__str__() == segmented_nt.replace(" ", "")
    segmented_alignment = "IGLV1-40*01|IGLJ7*01  CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGC AGCTCCAACATCGGGGCAGGTTAT\nfezakinumab           ...G.G.....C..C.....A..TAGT...AGC..T..A..T......C.T..G........T.........TC. TCT...........C..C.....C\n\nIGLV1-40*01|IGLJ7*01  GAT GTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTAT GGTAACAGC AATCGGCCCTCAGGGGTCCCTGACCGATTCTCTG\nfezakinumab           .GC ..G...........A........G..C..C........G..G........C ..CG..... .....T..A........T..G..T..C..TAGC.\n\nIGLV1-40*01|IGLJ7*01  GCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGC CAGTCCTATGACAGCAGCCTGAGTG\nfezakinumab           .G..T.....A..G.................G...........G.....G.....C.....A..C......... .....T.......ATTC.T......\n\nIGLV1-40*01|IGLJ7*01  GTTCTGCTGTG TTCGGAGGAGGCACCCAGCTGACCGTCCTCG\nfezakinumab           .---------C .....G.....G......T....T..T..-T\n\n"
    assert lambda_chain_nt.get_segmented_alignment_nt() == segmented_alignment
    assert lambda_chain_nt.locus == "IGL"

    assert lambda_chain_nt.cdr1_aa == "SSNIGAGYG"
    assert lambda_chain_nt.cdr2_aa == "GDS"
    assert lambda_chain_nt.cdr3_aa == "QSYDNSLSG"
    assert lambda_chain_nt.fwr1_aa == "QAVLTQPPSVSGAPGQRVTISCTGS"
    assert lambda_chain_nt.fwr2_aa == "VHWYQQLPGTAPKLLIY"
    assert lambda_chain_nt.fwr3_aa == "NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC"
    assert lambda_chain_nt.fwr4_aa == "FGGGTQLTVL"

    assert lambda_chain_nt.fwr1_aa_germline == "QSVLTQPPSVSGAPGQRVTISCTGS"
    assert lambda_chain_nt.cdr1_aa_germline == "SSNIGAGYD"
    assert lambda_chain_nt.fwr2_aa_germline == "VHWYQQLPGTAPKLLIY"
    assert lambda_chain_nt.cdr2_aa_germline == "GNS"
    assert lambda_chain_nt.fwr3_aa_germline == "NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC"
    assert lambda_chain_nt.fwr4_aa_germline == "FGGGTQLTVL"

    assert (
        lambda_chain_nt.get_segmented_vdj_aa()
        == "QAVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYG VHWYQQLPGTAPKLLIY GDS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDNSLSG FGGGTQLTVL"
    )
    segmented_alignment = "IGLV1-40*01|IGLJ7*01  QSVLTQPPSVSGAPGQRVTISCTGS SSNIGAGYD VHWYQQLPGTAPKLLIY GNS NRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC QSYDS\nfezakinumab           .A....................... ........G ................. .D. .................................... ....N\n\nIGLV1-40*01|IGLJ7*01  SLSGAV FGGGTQLTVL\nfezakinumab           ....-- ..........\n\n"
    assert lambda_chain_nt.get_segmented_alignment_aa() == segmented_alignment


def test_chain_io():
    chain_aa = antibody.AntibodyChainAA(
        fwr1_aa="QVQLKESGPGLVQPSQTLSLTCTVS",
        cdr1_aa="GLSLTSNS",
        fwr2_aa="VSWIRQPPGKGLEWMGV",
        cdr2_aa="IWSNGGT",
        fwr3_aa="DYNSAIESRLSINRDTSKSQVFLKMNSLQPEDTAMYFC",
        cdr3_aa="ASIYYYDADYLHWYFDF",
        fwr4_aa="WGPGTMVTVSS",
        v_gene="IGHV2-47",
        j_gene="IGHJ3",
        species="rat",
    )
    first_class_json = chain_aa.get_json()
    second_class = antibody.AntibodyChainAA.from_json(first_class_json)
    assert chain_aa == second_class

    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        chain_aa.to_json(t_file.name)
        second_class = antibody.AntibodyChainAA.read_json(t_file.name)
        assert chain_aa == second_class

    chain_nt = antibody.AntibodyChainNT(
        fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
        cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
        fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
        cdr2_nt="ATTTGGTCCAACGGTGGCACC",
        fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
        cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
        fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
        v_gene="IGHV2-47",
        j_gene="IGHJ3",
        species="rat",
    )

    first_class_json = chain_nt.get_json()
    second_class = antibody.AntibodyChainNT.from_json(first_class_json)
    assert chain_nt == second_class

    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        chain_nt.to_json(t_file.name)
        second_class = antibody.AntibodyChainNT.read_json(t_file.name)
        assert chain_nt == second_class

    kappa_chain_nt = antibody.KappaChainNT(
        fwr1_nt="GACATCCAAATGACACATTCGCCTTCATTGCTGAGTGCGTCTGTGGGTGACCGCGTCAGTCTGAACTGCAAGGCCTCC",
        cdr1_nt="CACTCAATCTACCGGAAT",
        fwr2_nt="CTGGCCTGGTACCAACAGAAACTCGGTGAGGCTCCAAAACTACTCATCTAC",
        cdr2_nt="AACGCCAAC",
        fwr3_nt="TCTCTGCAGACAGGAATCCCGTCTAGATTTAGCGGATCCGGCTCCGGTACCGACTTCACCCTGACCATTAGCTCCCTGCAGCCCGAGGATGTGGCGACCTATTTCTGC",
        cdr3_nt="CAACAGTACTATCGAGGATGGACG",
        fwr4_nt="TGGACGTTCGGTGGAGGTACAAAGCTGGAGCTG",
        v_gene="IGKV22S4",
        j_gene="IGKJ1",
        species="rat",
    )

    first_class_json = kappa_chain_nt.get_json()
    second_class = antibody.KappaChainNT.from_json(first_class_json)
    assert kappa_chain_nt == second_class

    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        kappa_chain_nt.to_json(t_file.name)
        second_class = antibody.KappaChainNT.read_json(t_file.name)
        assert kappa_chain_nt == second_class

    heavy_chain_nt = antibody.HeavyChainNT(
        fwr1_nt="CAGGTGCAGCTGAAGGAGAGCGGCCCTGGTTTGGTGCAGCCATCACAAACTCTTTCTCTGACATGCACCGTGTCA",
        cdr1_nt="GGCCTATCGCTCACCAGCAACTCC",
        fwr2_nt="GTCAGCTGGATACGTCAGCCGCCAGGCAAAGGTCTGGAGTGGATGGGTGTG",
        cdr2_nt="ATTTGGTCCAACGGTGGCACC",
        fwr3_nt="GACTACAACTCCGCTATCGAGAGCCGCTTGTCTATCAACCGCGACACCTCTAAATCCCAGGTTTTCTTGAAGATGAACTCGCTTCAACCTGAGGATACGGCTATGTACTTTTGC",
        cdr3_nt="GCCTCCATTTATTACTATGACGCTGACTACCTCCACTGGTACTTCGATTTC",
        fwr4_nt="TGGGGCCCCGGCACTATGGTGACCGTGAGCTCC",
        v_gene="IGHV2-47",
        j_gene="IGHJ3",
        species="rat",
    )

    first_class_json = heavy_chain_nt.get_json()
    second_class = antibody.HeavyChainNT.from_json(first_class_json)
    assert heavy_chain_nt == second_class

    with tempfile.NamedTemporaryFile(suffix="json") as t_file:
        heavy_chain_nt.to_json(t_file.name)
        second_class = antibody.HeavyChainNT.read_json(t_file.name)
        assert heavy_chain_nt == second_class
