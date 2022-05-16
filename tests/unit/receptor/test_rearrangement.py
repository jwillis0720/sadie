from pprint import pprint
from sadie.receptor.rearrangment import (
    InputSequence,
    PrimaryAnnotations,
    AlignmentAnnotations,
    AlignmentPositions,
    RegionSequences,
    RegionPositions,
    JunctionLengths,
    ReceptorChain,
)


def test_init_models():
    vrc01_heavy_seqeunce = "AAAAAACAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCTGGATATGAATTTATTGATTGTACGCTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGGCTGAAGCCTCGGGGGGGGGCCGTCAACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGGGGCCGGGGCACCCCGGTCATCGTCTCATCAGGGGGG"
    leader_and_vrc01_and_constant = "AAACTCA" + vrc01_heavy_seqeunce + "CACACACA"
    vrc01_aa_heavy = "QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS"

    # Model 1 - Input Dict
    input_dict = dict(
        sequence_id="test_sequence",
        sequence=vrc01_heavy_seqeunce,
        raw_sequence=leader_and_vrc01_and_constant,
        sequnce_aa=vrc01_aa_heavy,
    )
    input_sequence_model = InputSequence(**input_dict)
    pprint(input_sequence_model.__dict__)

    # Model 2 - Primary Annotations
    primary_sequence_annotation_model = PrimaryAnnotations(
        rev_comp=False,
        productive=True,
        vj_in_frame=True,
        stop_codon=False,
        complete_vdj=True,
        locus="IGH",
        v_call="IGHV1-2*02",
        d_call=["IGHD3-16*01", "IGHD3-16*02"],
        j_call="IGHJ1*01",
        v_call_top="IGHV1-2*02",
        d_call_top="IGHD3-16*01",
        j_call_top="IGHJ1*01",
        c_call="IGHG1*01",
    )
    pprint(primary_sequence_annotation_model.__dict__)

    # Model 3 - Alignment Annotations
    alignment_annotations_model = AlignmentAnnotations(
        sequence_alignment="CAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCTGGATATGAATTTATTGATTGTACGCTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGGCTGAAGCCTCGGGGGGGGGCCGTCAACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGGGGCCGGGGCACCCCGGTCATCGTCTCATCAG",
        sequence_alignment_aa="QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS",
        germline_alignment="CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGNNNNNNNNNNNNTGATTACGTTTGGGACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
        germline_alignment_aa="QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAXXXXXDYVWDFQHWGQGTLVTVSS",
        v_score=168.2,
        d_score=17.8,
        j_score=52.6,
        v_identity=0.6825,
        d_identity=0.85,
        j_identity=0.86,
        v_cigar="6S293M76S3N",
        d_cigar="311S6N14M50S17N",
        j_cigar="325S7N45M5S",
        v_support=6.796e-44,
        d_support=0.5755,
        j_support=5.727e-11,
        junction="TGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGG",
        junction_aa="CTRGKNCDYNWDFEHW",
        np1="GGGAAAAAACTG",
        c_score=100,
        c_identity=1,
        c_support=1e-44,
        c_cigar="6S293M76S3N",
    )
    # alignment_sequence_annotation_model = AlignmentAnnotations(**alignment_dict)
    pprint(alignment_annotations_model.__dict__)

    # Model 4 - Alignment Positions
    alignment_positions_dict = dict(
        v_sequence_start=7,
        v_sequence_end=299,
        v_germline_start=1,
        v_germline_end=293,
        v_alignment_start=1,
        v_alignment_end=293,
        d_sequence_start=312,
        d_sequence_end=325,
        d_germline_start=7,
        d_germline_end=20,
        d_alignment_start=306,
        d_alignment_end=319,
        j_sequence_start=326,
        j_sequence_end=370,
        j_germline_start=8,
        j_germline_end=52,
        j_alignment_start=320,
        j_alignment_end=364,
    )
    alignment_positions_model = AlignmentPositions(**alignment_positions_dict)
    pprint(alignment_positions_model.__dict__)

    # Model 5 - RegionSequences
    region_sequence_dict = dict(
        fwr="CAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCT",
        fwr1_aa="QVQLVQSGGQMKKPGESMRISCRAS",
        cdr1="GGATATGAATTTATTGATTGTACG",
        cdr1_aa="GYEFIDCT",
        fwr2="CTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGG",
        fwr2_aa="LNWIRLAPGKRPEWMGW",
        cdr2="CTGAAGCCTCGGGGGGGGGCCGTC",
        cdr2_aa="LKPRGGAV",
        fwr3="AACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGT",
        fwr3_aa="NYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFC",
        cdr3="ACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACAC",
        cdr3_aa="TRGKNCDYNWDFEH",
        fwr4="TGGGGCCGGGGCACCCCGGTCATCGTCTCATCA",
        fwr4_aa="WGRGTPVIVSS",
    )
    region_sequence_model = RegionSequences(**region_sequence_dict)
    pprint(region_sequence_model.__dict__)

    # Model 6 - RegionPositions
    region_positions_dict = dict(
        fwr1_start=7,
        fwr1_end=81,
        cdr1_start=82,
        cdr1_end=105,
        fwr2_start=106,
        fwr2_end=156,
        cdr2_start=157,
        cdr2_end=180,
        fwr3_start=181,
        fwr3_end=294,
        cdr3_start=295,
        cdr3_end=336,
        fwr4_start=337,
        fwr4_end=369,
    )
    region_position_model = RegionPositions(**region_positions_dict)
    pprint(region_position_model.__dict__)

    # Model 7 - JunctionLenghts
    junction_length_dict = dict(
        junction_length=48,
        junction_aa_length=None,
        np1_length=None,
        np2_length=None,
        np3_length=None,
        n1_length=None,
        n2_length=None,
        n3_length=None,
        p3v_length=None,
        p5d_length=None,
        p3d_length=None,
        p5d2_length=None,
        p3d2_length=None,
        p5j_length=None,
    )
    junction_length_model = JunctionLengths(**junction_length_dict)
    pprint(junction_length_model.__dict__)

    # Model 8 - ReceptorChain - required
    receptor_chain_model = ReceptorChain(
        input_sequence=input_sequence_model,
        primary_annotations=primary_sequence_annotation_model,
        alignment_annotations=alignment_annotations_model,
    )
    pprint(receptor_chain_model.__dict__)

    receptor_chain_model = ReceptorChain(
        input_sequence=input_sequence_model,
        primary_annotations=primary_sequence_annotation_model,
        alignment_annotations=alignment_annotations_model,
        alignment_positions=alignment_positions_model,
        region_sequences=region_sequence_model,
        region_positions=region_position_model,
        junction_lengths=junction_length_model,
    )
    print(receptor_chain_model)
