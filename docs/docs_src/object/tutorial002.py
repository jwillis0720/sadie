from pprint import pprint

from sadie.receptor.rearrangment import PrimaryAnnotations

# make a primary sequence model
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
# pretty print the dictionary attribute
pprint(primary_sequence_annotation_model.__dict__)
