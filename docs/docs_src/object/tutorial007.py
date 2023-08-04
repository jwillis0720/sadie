from pprint import pprint

from sadie.receptor.rearrangment import JunctionLengths

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
