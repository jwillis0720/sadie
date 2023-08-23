"""Unit tests for antibody."""

import logging

import pandas as pd
import pytest
from pkg_resources import resource_filename

from sadie import antibody
from sadie.antibody import exception, genetable, segment

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def test_gene_table():
    """
    testing the gene table class
    """
    fixture_v_gene = (
        pd.read_csv(fixture_file("VSEGMENT.csv.gz"), index_col=0)
        .fillna("")
        .sort_values(["species", "full"])
        .reset_index(drop=True)
    )
    v_gene_table = genetable.VGeneTable()

    # Assert everything is equal between our fixture and test frame
    pd._testing.assert_frame_equal(fixture_v_gene, v_gene_table.gene_table)

    assert list(v_gene_table.available_species) == [
        "alpaca",
        "boar",
        "camel",
        "cat",
        "catfish",
        "cod",
        "cow",
        "crabmacaque",
        "dog",
        "dolphin",
        "ferret",
        "goat",
        "horse",
        "human",
        "junglefowl",
        "macaque",
        "mouse",
        "nhp",
        "night_monkey",
        "platypus",
        "rabbit",
        "rat",
        "salmon",
        "sharks",
        "sheep",
        "teleosts",
        "trout",
        "zebrafish",
    ]

    assert list(v_gene_table.get_available_genes("human", locus="IGL")) == [
        "IGLV(I)-68*01",
        "IGLV(I)-70*01",
        "IGLV(III)-59-1*01",
        "IGLV(IV)-66-1*01",
        "IGLV(V)-66*01",
        "IGLV1-36*01",
        "IGLV1-40*01",
        "IGLV1-41*01",
        "IGLV1-44*01",
        "IGLV1-47*01",
        "IGLV1-50*01",
        "IGLV1-51*01",
        "IGLV1-62*01",
        "IGLV10-54*01",
        "IGLV10-67*01",
        "IGLV11-55*01",
        "IGLV2-11*01",
        "IGLV2-14*01",
        "IGLV2-18*01",
        "IGLV2-23*01",
        "IGLV2-28*01",
        "IGLV2-33*01",
        "IGLV2-34*01",
        "IGLV2-5*01",
        "IGLV2-8*01",
        "IGLV3-1*01",
        "IGLV3-10*01",
        "IGLV3-12*01",
        "IGLV3-13*01",
        "IGLV3-16*01",
        "IGLV3-17*01",
        "IGLV3-19*01",
        "IGLV3-2*01",
        "IGLV3-21*01",
        "IGLV3-22*01",
        "IGLV3-24*01",
        "IGLV3-25*01",
        "IGLV3-26*01",
        "IGLV3-27*01",
        "IGLV3-29*01",
        "IGLV3-30*01",
        "IGLV3-31*01",
        "IGLV3-32*01",
        "IGLV3-6*01",
        "IGLV3-7*01",
        "IGLV3-9*01",
        "IGLV4-3*01",
        "IGLV4-60*01",
        "IGLV4-69*01",
        "IGLV5-37*01",
        "IGLV5-39*01",
        "IGLV5-45*01",
        "IGLV5-48*01",
        "IGLV5-52*01",
        "IGLV6-57*01",
        "IGLV7-35*01",
        "IGLV7-43*01",
        "IGLV7-46*01",
        "IGLV8-61*01",
        "IGLV8/OR8-1*01",
        "IGLV9-49*01",
    ]

    v323_human = v_gene_table.get_gene("IGHV3-23*01", "human")
    assert isinstance(v323_human, pd.Series)
    assert list(v323_human.index) == [
        "gene",
        "fwr1_nt",
        "cdr1_nt",
        "fwr2_nt",
        "cdr2_nt",
        "fwr3_nt",
        "cdr3_nt",
        "fwr1_aa",
        "cdr1_aa",
        "fwr2_aa",
        "cdr2_aa",
        "fwr3_aa",
        "cdr3_aa",
    ]
    with pytest.raises(exception.AmbiguousGene):
        v_gene_table.get_gene("IGHV3-23", "human")
        v_gene_table.get_gene("IGHV1-69", "human")
    with pytest.raises(exception.BadGene):
        v_gene_table.get_gene("IGHV3-700", "human")
        v_gene_table.get_gene("IGHZ3-33", "human")
    with pytest.raises(exception.NoSpecies):
        v_gene_table.get_gene("IGHV3-23", "robot")
        v_gene_table.get_gene("IGHV1-69", "elon")

    fixture_j_gene = (
        pd.read_csv(fixture_file("JSEGMENT.csv.gz"), index_col=0)
        .fillna("")
        .sort_values(["species", "full"])
        .reset_index(drop=True)
    )
    j_gene_table = genetable.JGeneTable()

    # Assert everything is equal between our fixture and test frame
    pd._testing.assert_frame_equal(fixture_j_gene, j_gene_table.gene_table)

    j6_human = j_gene_table.get_gene("IGHJ6*01", "human")
    assert isinstance(j6_human, pd.Series)
    with pytest.raises(exception.AmbiguousGene):
        j_gene_table.get_gene("IGHJ4", "human")
        j_gene_table.get_gene("IGHJ6", "human")
    with pytest.raises(exception.NoSpecies):
        j_gene_table.get_gene("IGHJ4", "robot")

    with pytest.raises(exception.BadGene):
        j_gene_table.get_gene("IGHJ12", "dog")


def test_v_gene():
    # test we can create v gene objects
    v_segments_table_manual = pd.read_csv(fixture_file("VSEGMENT.csv.gz"), index_col=0).fillna("")
    for gb_index, v_gene_df in v_segments_table_manual.groupby(["species", "gene"]):
        species = gb_index[0]
        v_gene = gb_index[1]
        if len(v_gene_df) == 1:
            lookup = v_gene_df.iloc[0]
            # we have a single allelic match
            cdr1_aa = lookup["cdr1_aa"]
            cdr2_aa = lookup["cdr2_aa"]
            cdr3_aa = lookup["cdr3_aa"]
            fwr1_aa = lookup["fwr1_aa"]
            fwr2_aa = lookup["fwr2_aa"]
            fwr3_aa = lookup["fwr3_aa"]
            #  cdr nt
            cdr1_nt = lookup["cdr1_nt"]
            cdr2_nt = lookup["cdr2_nt"]
            cdr3_nt = lookup["cdr3_nt"]
            fwr1_nt = lookup["fwr1_nt"]
            fwr2_nt = lookup["fwr2_nt"]
            fwr3_nt = lookup["fwr3_nt"]
            v_gene_object = antibody.VGene(v_gene, species)

            # Make sure we are getting the right subclass
            assert isinstance(v_gene_object, antibody.VGene)
            assert isinstance(v_gene_object.cdr1_aa, segment.CDR1AA)
            assert isinstance(v_gene_object.cdr2_aa, segment.CDR2AA)
            assert isinstance(v_gene_object.cdr3_aa, segment.CDR3AA)
            assert isinstance(v_gene_object.cdr1_nt, segment.CDR1NT)
            assert isinstance(v_gene_object.cdr2_nt, segment.CDR2NT)
            assert isinstance(v_gene_object.cdr3_nt, segment.CDR3NT)
            assert isinstance(v_gene_object.fwr1_aa, segment.FrameWork1AA)
            assert isinstance(v_gene_object.fwr2_aa, segment.FrameWork2AA)
            assert isinstance(v_gene_object.fwr3_aa, segment.FrameWork3AA)
            assert isinstance(v_gene_object.fwr1_nt, segment.FrameWork1NT)
            assert isinstance(v_gene_object.fwr2_nt, segment.FrameWork2NT)
            assert isinstance(v_gene_object.fwr3_nt, segment.FrameWork3NT)

            # assert we have the same objects strings
            assert v_gene_object.cdr1_aa == cdr1_aa
            assert v_gene_object.cdr2_aa == cdr2_aa
            assert v_gene_object.cdr3_aa == cdr3_aa
            assert v_gene_object.fwr1_aa == fwr1_aa
            assert v_gene_object.fwr2_aa == fwr2_aa
            assert v_gene_object.fwr3_aa == fwr3_aa
            assert v_gene_object.cdr1_nt == cdr1_nt
            assert v_gene_object.cdr2_nt == cdr2_nt
            assert v_gene_object.cdr3_nt == cdr3_nt
            assert v_gene_object.fwr1_nt == fwr1_nt
            assert v_gene_object.fwr2_nt == fwr2_nt
            assert v_gene_object.fwr3_nt == fwr3_nt
        else:
            # else test that we can raise an ambigous gene for those that have more than one allele
            with pytest.raises(exception.AmbiguousGene):
                v_gene_object = antibody.VGene(v_gene, species)
            for index in range(0, len(v_gene_df)):
                lookup = v_gene_df.iloc[index]
                # we have to use the full allel
                v_gene = lookup["full"]
                # we have a single allelic match
                cdr1_aa = lookup["cdr1_aa"]
                cdr2_aa = lookup["cdr2_aa"]
                cdr3_aa = lookup["cdr3_aa"]
                fwr1_aa = lookup["fwr1_aa"]
                fwr2_aa = lookup["fwr2_aa"]
                fwr3_aa = lookup["fwr3_aa"]
                # cdr nt
                cdr1_nt = lookup["cdr1_nt"]
                cdr2_nt = lookup["cdr2_nt"]
                cdr3_nt = lookup["cdr3_nt"]
                fwr1_nt = lookup["fwr1_nt"]
                fwr2_nt = lookup["fwr2_nt"]
                fwr3_nt = lookup["fwr3_nt"]
                v_gene_object = antibody.VGene(v_gene, species)
                # Make sure we are getting the right subclass
                assert isinstance(v_gene_object, antibody.VGene)
                assert isinstance(v_gene_object.cdr1_aa, segment.CDR1AA)
                assert isinstance(v_gene_object.cdr2_aa, segment.CDR2AA)
                assert isinstance(v_gene_object.cdr3_aa, segment.CDR3AA)
                assert isinstance(v_gene_object.cdr1_nt, segment.CDR1NT)
                assert isinstance(v_gene_object.cdr2_nt, segment.CDR2NT)
                assert isinstance(v_gene_object.cdr3_nt, segment.CDR3NT)
                assert isinstance(v_gene_object.fwr1_aa, segment.FrameWork1AA)
                assert isinstance(v_gene_object.fwr2_aa, segment.FrameWork2AA)
                assert isinstance(v_gene_object.fwr3_aa, segment.FrameWork3AA)
                assert isinstance(v_gene_object.fwr1_nt, segment.FrameWork1NT)
                assert isinstance(v_gene_object.fwr2_nt, segment.FrameWork2NT)
                assert isinstance(v_gene_object.fwr3_nt, segment.FrameWork3NT)

                # assert we have the same objects strings
                assert v_gene_object.cdr1_aa == cdr1_aa
                assert v_gene_object.cdr2_aa == cdr2_aa
                assert v_gene_object.cdr3_aa == cdr3_aa
                assert v_gene_object.fwr1_aa == fwr1_aa
                assert v_gene_object.fwr2_aa == fwr2_aa
                assert v_gene_object.fwr3_aa == fwr3_aa
                assert v_gene_object.cdr1_nt == cdr1_nt
                assert v_gene_object.cdr2_nt == cdr2_nt
                assert v_gene_object.cdr3_nt == cdr3_nt
                assert v_gene_object.fwr1_nt == fwr1_nt
                assert v_gene_object.fwr2_nt == fwr2_nt
                assert v_gene_object.fwr3_nt == fwr3_nt


def test_v_gene_numbering():
    v_gene_object = antibody.VGene("IGHV3-23*01", "human")
    assert v_gene_object.fwr1_aa.start_index == 0
    assert v_gene_object.fwr1_aa.start == 1
    assert v_gene_object.fwr1_aa.end == 25
    assert v_gene_object.fwr1_aa.end_index == 24

    assert v_gene_object.cdr1_aa.start_index == 25
    assert v_gene_object.cdr1_aa.start == 26
    assert v_gene_object.cdr1_aa.end == 33
    assert v_gene_object.cdr1_aa.end_index == 32

    assert v_gene_object.fwr2_aa.start_index == 33
    assert v_gene_object.fwr2_aa.start == 34
    assert v_gene_object.fwr2_aa.end == 50
    assert v_gene_object.fwr2_aa.end_index == 49

    assert v_gene_object.cdr2_aa.start_index == 50
    assert v_gene_object.cdr2_aa.start == 51
    assert v_gene_object.cdr2_aa.end == 58
    assert v_gene_object.cdr2_aa.end_index == 57

    assert v_gene_object.fwr3_aa.start_index == 58
    assert v_gene_object.fwr3_aa.start == 59
    assert v_gene_object.fwr3_aa.end == 96
    assert v_gene_object.fwr3_aa.end_index == 95

    assert v_gene_object.fwr1_nt.start_index == 0
    assert v_gene_object.fwr1_nt.start == 1
    assert v_gene_object.fwr1_nt.end == 75
    assert v_gene_object.fwr1_nt.end_index == 74

    assert v_gene_object.cdr1_nt.start_index == 75
    assert v_gene_object.cdr1_nt.start == 76
    assert v_gene_object.cdr1_nt.end == 99
    assert v_gene_object.cdr1_nt.end_index == 98

    assert v_gene_object.fwr2_nt.start_index == 99
    assert v_gene_object.fwr2_nt.start == 100
    assert v_gene_object.fwr2_nt.end == 150
    assert v_gene_object.fwr2_nt.end_index == 149

    assert v_gene_object.cdr2_nt.start_index == 150
    assert v_gene_object.cdr2_nt.start == 151
    assert v_gene_object.cdr2_nt.end == 174
    assert v_gene_object.cdr2_nt.end_index == 173

    assert v_gene_object.fwr3_nt.start_index == 174
    assert v_gene_object.fwr3_nt.start == 175
    assert v_gene_object.fwr3_nt.end == 288
    assert v_gene_object.fwr3_nt.end_index == 287

    rep = "<VGene>human->IGHV3-23*01\n        fwr1_nt-GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT\n        fwr1_aa-EVQLLESGGGLVQPGGSLRLSCAAS\n        cdr1_nt-GGATTCACCTTTAGCAGCTATGCC\n        cdr1_aa-GFTFSSYA\n        fwr2_nt-ATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCT\n        fwr2_aa-MSWVRQAPGKGLEWVSA\n        cdr2_nt-ATTAGTGGTAGTGGTGGTAGCACA\n        cdr2_aa-ISGSGGST\n        fwr3_nt-TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGT\n        fwr3_aa-YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC\n        cdr3_nt-GCGAAAGA\n        cdr3_aa-AK\n        "
    assert v_gene_object.__repr__() == rep

    # alternative numbering
    with pytest.raises(NotImplementedError):
        antibody.VGene("IGHV3-23*01", "human", "kabat")
    # assert v_gene_kabat.cdr1_aa == "SYAMS"
    # assert v_gene_kabat.cdr2_aa == "AISGSGGSTYYADSVKG"
    # assert v_gene_kabat.cdr3_aa == "K"


def test_j_gene():
    # test we can create j gene objects
    j_segments_table_manual = pd.read_csv(fixture_file("JSEGMENT.csv.gz"), index_col=0).fillna("")
    for gb_index, j_gene_df in j_segments_table_manual.groupby(["species", "gene"]):
        species = gb_index[0]
        j_gene = gb_index[1]
        if len(j_gene_df) == 1:
            lookup = j_gene_df.iloc[0]
            # we have a single allelic match
            cdr3_nt = lookup["cdr3_nt"]
            cdr3_aa = lookup["cdr3_aa"]

            # FW4
            fwr4_nt = lookup["fwr4_nt"]
            fwr4_aa = lookup["fwr4_aa"]

            # J Gene Object
            j_gene_object = antibody.JGene(j_gene, species)

            # Make sure we are getting the right subclass
            assert isinstance(j_gene_object, antibody.JGene)
            assert isinstance(j_gene_object.cdr3_nt, segment.CDR3NT)
            assert isinstance(j_gene_object.cdr3_aa, segment.CDR3AA)
            assert isinstance(j_gene_object.fwr4_nt, segment.FrameWork4NT)
            assert isinstance(j_gene_object.fwr4_aa, segment.FrameWork4AA)

            # assert we have the same objects strings
            assert j_gene_object.cdr3_nt == cdr3_nt
            assert j_gene_object.cdr3_aa == cdr3_aa
            assert j_gene_object.fwr4_nt == fwr4_nt
            assert j_gene_object.fwr4_aa == fwr4_aa
    j_gene_object = antibody.JGene("IGHJ6*01", "human")
    rep = "<JGene>human->IGHJ6*01\n        cdr3_nt-ATTACTACTACTACTACGGTATGGACGTC\n        cdr3_aa-YYYYYGMDV\n        fwr4_nt-TGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG\n        fwr4_aa-WGQGTTVTVSS\n        "
    assert j_gene_object.__repr__() == rep
