from distutils.version import StrictVersion
from math import nan
from pathlib import Path

import pandas as pd
from pkg_resources import resource_filename
from sadie.airr import constants, AirrTable, Airr


def fixture_file(file):
    """Helper method for test execution."""
    _file = Path(resource_filename(__name__, f"fixtures/{file}"))
    if not _file.exists():
        raise FileExistsError(f"Fixutre file not found {_file}")
    return _file


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


def make_sadie_comparable(df):
    """Takes sadie df and makes it comparable wiht IMGT

    Parameters
    ----------
    df : AirrTable

    Returns
    -------
    pd.DataFrame
        returns a pandas dataframe for comparison
    """
    if isinstance(df, AirrTable):
        df = df.table

    # comparison keys between imgt and sadie
    compare_key = list(constants.IGBLAST_AIRR.keys()) + [
        "v_call_top",
        "d_call_top",
        "j_call_top",
    ]

    # Drop frameshift since imgt does not have it
    compare_key.remove("v_frameshift")

    # Just get compare keys
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


def make_imgt_comparable(df: pd.DataFrame) -> pd.DataFrame:

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


def test_imgt_integration():
    file = fixture_file("OAS_subsample_good_anarci.fasta.gz")
    airr_api = Airr(species="human", database="imgt", functional="all")
    # try 1 with -1, -2
    # These are the default
    airr_api.igblast.v_penalty = -1
    airr_api.igblast.j_penalty = -2
    print(file)


def test_catnap_heavy_integration():
    file = fixture_file("catnap_nt_heavy.fasta.gz")
    airr_api = Airr(species="human", database="imgt", functional="all")
    # try 1 with -1, -2
    # These are the default
    airr_api.igblast.v_penalty = -1
    airr_api.igblast.j_penalty = -2
    result_airr = airr_api.run_file(file)
    print(result_airr)


# def test_vs_oas_sampling():
#     """
#     Test airr vs the OAS airr datasets. I have chosen 15K heavy and 10K light to compare against
#     """
#     end_point = "https://sadie.s3.us-east-2.amazonaws.com/integration/OAS_sample_subsample.bz2"

#     # uzing bzip for max compression to minimize egress S3 expense
#     df_bz2 = pd.read_csv(end_point, index_col=0).reset_index()

#     # Airr api
#     airr_api = Airr(species="human", database="imgt", functional="all")
#     airr_api.run_dataframe(
#         df_bz2,
#     )


# def test_cli():
#     """Confirm the CLI works as expecte"""
#     runner = CliRunner()
#     fastas = get_fasta_inputs("*")
#     for query_file in fastas:
#         with tempfile.NamedTemporaryFile() as tmpfile:
#             logger.debug(tmpfile.name)
#             result = runner.invoke(app.run_airr, ["--query", query_file, "-o", tmpfile.name])
#             assert result.exit_code == 0
#             assert os.path.exists(tmpfile.name)
#     assert not os.path.exists(tmpfile.name)
#     for query_file in fastas:
#         with tempfile.NamedTemporaryFile() as tmpfile:
#             logger.debug(tmpfile.name)
#             result = runner.invoke(app.run_airr, ["--query", query_file, "-s", "dog", "-o", tmpfile.name])
#             assert result.exit_code == 0
#             assert os.path.exists(tmpfile.name)
#         assert not os.path.exists(tmpfile.name)
