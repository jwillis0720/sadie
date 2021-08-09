"""Unit tests for airrtable."""

import pytest


import pandas as pd

from sadie.airr import AirrTable, LinkedAirrTable
from sadie.airr.exceptions import MissingAirrColumns
from Bio.SeqRecord import SeqRecord


def test_airrtable_init(fixture_setup):
    test_csv = fixture_setup.get_dog_airrdataframe_file()
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
    test_json = fixture_setup.get_json_as_dataframe()
    airr_table = AirrTable.read_json(test_json)
    assert not airr_table.empty

    # gen bank
    genbanks = airr_table.get_genbank()
    assert all([isinstance(i, SeqRecord) for i in genbanks])

    # I will not accept a busted table sam I am
    busted_table = fixture_setup.get_busted_airrtable()
    with pytest.raises(MissingAirrColumns) as e:
        AirrTable.read_csv(busted_table)
    assert e.value.__str__()


def test_indel_correction(fixture_setup):
    test_csv = fixture_setup.get_dog_airrdataframe_file()
    # Test we can initllize with staic meathod
    airr_table = AirrTable.read_csv(test_csv)
    airr_table_indel = AirrTable.correct_indel(airr_table)
    assert "germline_alignment_aa_corrected" in airr_table_indel.columns
    assert "v_germline_alignment_aa_corrected" in airr_table_indel.columns
    assert isinstance(airr_table_indel, AirrTable)


def test_scfv_airrtable(fixture_setup):
    file_path = fixture_setup.get_dummy_scfv_table()
    dummy_scfv_table = pd.read_csv(file_path, index_col=0)
    linked_table = LinkedAirrTable(dummy_scfv_table)
    # test if we can split
    heavy_table, light_table = linked_table.get_split_table()
    assert isinstance(heavy_table, AirrTable)
    assert isinstance(light_table, AirrTable)
    rebuild_table = LinkedAirrTable(heavy_table.merge(light_table, on="sequence_id", suffixes=["_heavy", "_light"]))
    assert rebuild_table == rebuild_table
    assert rebuild_table == linked_table

    heavy_table["cell_id"] = heavy_table["sequence_id"]
    light_table["cell_id"] = light_table["sequence_id"]
    with pytest.raises(MissingAirrColumns):
        LinkedAirrTable(heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]))
    rebuild_data = LinkedAirrTable(
        heavy_table.merge(light_table, on="cell_id", suffixes=["_heavy", "_light"]), key_column="cell_id"
    )
    heavy_table_split, light_table_split = rebuild_data.get_split_table()
    assert heavy_table_split.columns.difference(heavy_table.columns).empty
    assert light_table_split.columns.difference(light_table.columns).empty
    assert heavy_table == heavy_table_split[heavy_table.columns]
    assert light_table == light_table_split[light_table.columns]
    assert rebuild_data.key_column == "cell_id"
    assert rebuild_data.suffixes == ["_heavy", "_light"]

    rebuild_data = LinkedAirrTable(
        heavy_table.merge(light_table, on="cell_id", suffixes=["_h", "_l"]), suffixes=["_h", "_l"], key_column="cell_id"
    )
    assert rebuild_data.suffixes == ["_h", "_l"]
