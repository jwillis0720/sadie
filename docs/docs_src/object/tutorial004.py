from pprint import pprint

from sadie.receptor.rearrangment import AlignmentPositions

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
# pretty print dictonary
pprint(alignment_positions_model.__dict__)
