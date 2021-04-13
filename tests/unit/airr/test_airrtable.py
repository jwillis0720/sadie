"""Unit tests for airrtable."""

import logging
import pytest
from pkg_resources import resource_filename


import pandas as pd

from sadie.airr import AirrTable, LinkedAirrTable
from sadie.airr.exceptions import MissingAirrColumns
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
    non_airr_columns = [
        "d_call_top",
        "d_mutation",
        "d_mutation_aa",
        "j_call_top",
        "j_mutation",
        "j_mutation_aa",
        "liable",
        "note",
        "species",
        "v_call_top",
        "v_mutation",
        "v_mutation_aa",
        "vdj_aa",
        "vdj_nt",
    ]
    assert sorted(airr_table.non_airr_columns) == non_airr_columns
    # Tes twe can init with a direct pandas call
    airr_table = AirrTable(pd.read_csv(test_csv))
    assert not airr_table.empty
    assert isinstance(airr_table, AirrTable)

    # test we can read in json too
    test_json = fixture_file("airr_tables/heavy_sample.json.gz")
    airr_table = AirrTable.read_json(test_json)
    assert not airr_table.empty

    # gen bank
    genbanks = airr_table.get_genbank()
    assert all([isinstance(i, SeqRecord) for i in genbanks])

    # I will not accept a busted table sam I am
    busted_table = fixture_file("airr_tables/busted.csv.gz")
    with pytest.raises(MissingAirrColumns) as e:
        AirrTable.read_csv(busted_table)
    assert e.value.__str__()


def test_indel_correction():
    test_csv = fixture_file("airr_tables/dog_igh.csv.gz")
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_csv(test_csv)
    airr_table_indel = AirrTable.correct_indel(airr_table)
    assert "germline_alignment_aa_corrected" in airr_table_indel.columns
    assert "v_germline_alignment_aa_corrected" in airr_table_indel.columns
    assert isinstance(airr_table_indel, AirrTable)


def test_scfv_airrtable():
    file_path = fixture_file("airr_tables/linked_airr_table_dummy.csv.gz")
    dummy_scfv_table = pd.read_csv(file_path, index_col=0)
    linked_table = LinkedAirrTable(dummy_scfv_table)
    # test if we can split
    heavy_table, light_table = linked_table.get_split_table()
    assert isinstance(heavy_table, AirrTable)
    assert isinstance(light_table, AirrTable)
    rebuild_table = LinkedAirrTable(heavy_table.merge(light_table, on="sequence_id", suffixes=["_heavy", "_light"]))
    assert rebuild_table == rebuild_table
    assert rebuild_table == linked_table
