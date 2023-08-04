from pprint import pprint

from sadie.receptor.rearrangment import RegionPositions

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
