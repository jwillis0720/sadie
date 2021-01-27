"""Unit tests for airrtable."""

import logging
import pytest
from pkg_resources import resource_filename


import pandas as pd

from sadie.airr import AirrTable, Airr, ScfvAirrTable
from sadie.airr.airrtable import MissingAirrColumns
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger()


def fixture_file(file):
    """Helper method for test execution."""
    return resource_filename(__name__, "fixtures/{}".format(file))


def test_airrtable_init():
    test_csv = fixture_file("airr_tables/dog_igh.csv.gz")
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_csv(test_csv)
    assert not airr_table.empty
    assert sorted(airr_table.non_airr_columns) == [
        "d_call_top",
        "j_call_top",
        "note",
        "species",
        "v_call_top",
        "vdj_aa",
        "vdj_nt",
    ]
    # Tes twe can init with a direct pandas call
    airr_table = AirrTable(pd.read_csv(test_csv))
    assert not airr_table.empty
    assert isinstance(airr_table, AirrTable)

    # test we can read in json too
    test_json = fixture_file("airr_tables/heavy_sample.json.gz")
    airr_table = AirrTable.read_json(test_json)
    assert sorted(airr_table.non_airr_columns) == [
        "cdr3_aa_length",
        "chain",
        "d_call_top",
        "d_family",
        "j_call_top",
        "j_family",
        "note",
        "species",
        "v_call_top",
        "v_family",
        "vdj_aa",
        "vdj_nt",
    ]
    assert not airr_table.empty

    # gen bank
    genbanks = airr_table.get_genbank()
    assert all([isinstance(i, SeqRecord) for i in genbanks])

    # I will not accept a busted table sam I am
    busted_table = fixture_file("airr_tables/busted.csv.gz")
    with pytest.raises(MissingAirrColumns) as e:
        AirrTable.read_csv(busted_table)
    assert e.value.__str__()


def test_scfv_airrtable():
    airr_api = Airr("human")
    scfv_airr_table = airr_api.run_file(fixture_file("fasta_inputs/scfv.fasta"), scfv=True)
    heavy, light = ScfvAirrTable.deconstruct_scfv(scfv_airr_table)
    assert type(heavy) == AirrTable
    assert type(light) == AirrTable
