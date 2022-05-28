from distutils.version import StrictVersion
from math import nan

import pandas as pd
from sadie.airr import Airr, AirrTable
from sadie.airr.airrtable import constants
from tests.conftest import SadieFixture


def fillna(df, fill_value=""):
    """
    Replace null values with `fill_value`.

    Also replaces in categorical columns.
    """
    for col in df.dtypes[df.dtypes == "category"].index:
        if fill_value not in df[col].cat.categories:
            df[col].cat.add_categories([fill_value], inplace=True)
    # Known bug https://github.com/pandas-dev/pandas/issues/25472
    if StrictVersion(pd.__version__) >= StrictVersion("1.0"):
        for col in df.dtypes[df.dtypes.apply(lambda x: x in ["float64", "Int16", "Int64"])].index:
            df[col] = df[col].astype("float")
    return df.fillna(fill_value)


# straight ignore these
ignore = [
    "v_score",
    "d_score",
    "j_score",
    "v_support",
    "d_support",
    "j_support",
]

# cast these to integers
starts_and_ends = [
    "cdr1_end",
    "cdr1_start",
    "cdr2_end",
    "cdr2_start",
    "cdr3_end",
    "cdr3_start",
    "d_alignment_end",
    "d_alignment_start",
    "d_germline_end",
    "d_germline_start",
    "d_sequence_end",
    "d_sequence_start",
    "fwr1_end",
    "fwr1_start",
    "fwr2_end",
    "fwr2_start",
    "fwr4_end",
    "fwr4_start",
    "j_alignment_end",
    "j_alignment_start",
    "j_germline_end",
    "j_germline_start",
    "j_sequence_end",
    "j_sequence_start",
]

# check these in integration
check_these = [
    "sequence",
    "locus",
    "stop_codon",
    "vj_in_frame",
    "productive",
    "rev_comp",
    "complete_vdj",
    "fwr1",
    "cdr1",
    "fwr2",
    "cdr2",
    "fwr3",
    "cdr3",
    "fwr1_aa",
    "cdr1_aa",
    "fwr2_aa",
    "cdr2_aa",
    "fwr3_aa",
    "cdr3_aa",
    "fwr4_aa",
]


def _make_sadie_comparable(df):
    """Takes sadie df and makes it comparable wiht IMGT

    Parameters
    ----------
    df : AirrTable

    Returns
    -------
    pd.DataFrame
        returns a pandas dataframe for comparison
    """

    # comparison keys between imgt and sadie
    compare_key = list(constants.IGBLAST_AIRR.keys()) + [
        "v_call_top",
        "d_call_top",
        "j_call_top",
    ]

    # Drop frameshift since imgt does not have it
    compare_key.remove("v_frameshift")

    # Just get compare keys
    df = pd.DataFrame(df)
    df = df[compare_key].drop(ignore, axis=1)
    df.loc[:, starts_and_ends] = df[starts_and_ends].astype("Int64")
    # Just get the gene top call IGHV1-2*01 -> IGHV1-2
    df.insert(
        df.columns.get_loc("v_call_top"),
        "v_gene_top",
        df["v_call_top"].str.split("*").str.get(0),
    )
    df.insert(
        df.columns.get_loc("d_call_top"),
        "d_gene_top",
        df["d_call_top"].str.split("*").str.get(0),
    )
    df.insert(
        df.columns.get_loc("j_call_top"),
        "j_gene_top",
        df["j_call_top"].str.split("*").str.get(0),
    )
    return df


def _make_imgt_comparable(df: pd.DataFrame) -> pd.DataFrame:

    """Takes Hi-Vquest and return a compariable dataframe

    Parameters
    ----------
    df : pd.Dataframe

    Returns
    -------
    pd.Dataframe
        A comparable dataframe
    """

    def map_bool(mapp_df, col):
        return mapp_df.loc[:, col].map(
            {"T": True, "F": False, True: True, False: False, nan: False},
            na_action="ignore",
        )

    compare_key = list(constants.IGBLAST_AIRR.keys())
    compare_key.remove("v_frameshift")
    df = df[compare_key].copy()
    bool_cols = ["vj_in_frame", "productive", "rev_comp", "complete_vdj", "stop_codon"]
    for col in bool_cols:
        df[col] = map_bool(df, col)

    upper_columns = [
        "fwr1",
        "cdr1",
        "fwr2",
        "cdr2",
        "fwr3",
        "cdr3",
        "fwr4",
        "germline_alignment",
        "germline_alignment_aa",
        "sequence_alignment",
        "sequence_alignment_aa",
        "v_germline_alignment",
        "v_germline_alignment_aa",
        "v_sequence_alignment",
        "v_sequence_alignment_aa",
        "d_germline_alignment",
        "d_germline_alignment_aa",
        "d_sequence_alignment",
        "d_sequence_alignment_aa",
        "j_germline_alignment",
        "j_germline_alignment_aa",
        "j_sequence_alignment",
        "j_sequence_alignment_aa",
        "np1",
        "np2",
        "junction",
        "sequence",
    ]
    for col in upper_columns:
        df[col] = df[col].str.upper()

    # strip the "." on left and right by imgt
    df["fwr1"] = df["fwr1"].str.lstrip(".")
    df["fwr4"] = df["fwr4"].str.rstrip(".")

    # Get top gens
    df["v_call_top"] = df["v_call"].str.split(",").str.get(0).str.split().str.get(1)
    df["v_gene_top"] = df["v_call_top"].str.split("*").str.get(0)
    df["d_call_top"] = df["d_call"].str.split(",").str.get(0).str.split("_").str.get(-1)
    df["d_gene_top"] = df["d_call_top"].str.split("*").str.get(0)
    df["j_call_top"] = df["j_call"].str.split(",").str.get(0).str.split().str.get(1)
    df["j_gene_top"] = df["j_call_top"].str.split("*").str.get(0)

    # Drop the ignore
    df = df.drop(ignore, axis=1)

    # Convert the integer types
    df.loc[:, starts_and_ends] = df[starts_and_ends].astype("Int64")
    return df


def test_imgt_integration(fixture_setup: SadieFixture) -> None:
    # 1 of 2 major airr integerations. Checks that we match imgt on the ones we care about
    ignore = [
        48,
        67,
        143,
        149,
        193,
        213,
        239,
        286,
        305,
        306,
        364,
        367,
        383,
        419,
        436,
        450,
        457,
        486,
        490,
        520,
        521,
        590,
        606,
        612,
        631,
        698,
        715,
        720,
        732,
        760,
        767,
        827,
        839,
        888,
        899,
        983,
        989,
    ]

    file = fixture_setup.get_oas_fasta()
    airr_api = Airr("human", adaptable=True)
    sadie_airr = airr_api.run_fasta(file)
    sadie_comparable = _make_sadie_comparable(sadie_airr)[check_these]
    sadie_comparable = sadie_comparable.drop(pd.Index(ignore))

    imgt_df = pd.read_csv(fixture_setup.get_imgt_airrtable(), low_memory=False)
    imgt_comparable = _make_imgt_comparable(imgt_df)[check_these]
    imgt_comparable = imgt_comparable.drop(pd.Index(ignore))
    for x in check_these:
        if not (sadie_comparable[x] == imgt_comparable[x]).all():
            not_index = sadie_comparable[~(sadie_comparable[x] == imgt_comparable[x])].index
            if not not_index.empty:
                raise AssertionError(f"{not_index} does not match {x}")


# @pytest.mark.skip(reason="integration tests will change under this active development")
def test_catnap_integration(fixture_setup: SadieFixture) -> None:
    # 2 of 2 major integrations. Make sure our catnap calls dont change from last time we ran it
    heavy_file = fixture_setup.get_catnap_heavy_nt()
    light_file = fixture_setup.get_catnap_light_nt()
    airr_api = Airr("human", adaptable=True)
    heavy_at = AirrTable(pd.read_feather(fixture_setup.get_catnap_heavy_airrtable()))
    light_at = AirrTable(pd.read_feather(fixture_setup.get_catnap_light_airrtable()))
    catnap_heavy = airr_api.run_fasta(heavy_file)
    catnap_light = airr_api.run_fasta(light_file)
    diffs = catnap_light.columns.symmetric_difference(light_at.columns)
    if not diffs.empty:
        raise AssertionError(f"Light table has the following different columns {diffs}")
    diffs = catnap_heavy.columns.symmetric_difference(heavy_at.columns)
    if not diffs.empty:
        raise AssertionError(f"Heavy table has the following different columns {diffs}")
    pd.testing.assert_frame_equal(light_at, catnap_light, check_dtype=False, rtol=0.15)
    pd.testing.assert_frame_equal(heavy_at, catnap_heavy, check_dtype=False, rtol=0.15)
