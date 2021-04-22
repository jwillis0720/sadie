from typing import Union
import logging

# third party
import pandas as pd

# module/package level
from sadie.airr.airrtable import AirrTable, LinkedAirrTable
from sadie.anarci import Anarci

logger = logging.getLogger("AirrMethod")


def run_mutational_analysis(
    airrtable: Union[AirrTable, LinkedAirrTable], scheme: str
) -> Union[AirrTable, LinkedAirrTable]:
    """Run a mutational analysis given a numbering scheme. Returns an AirrTable with added mutational analysis columns

    This method is computationally expensive.

    Parameters
    ----------
    airrtable : Union[AirrTable,LinkedAirrTable]
        An AirrTable or LinkedAirrTable
    scheme : str
        the numbering scheme: ex, 'martin','kabat','imgt','chothia'

    Returns
    -------
    AirrTable
        returns an airrtable with mutation and scheme fields containing the germline mutations

    Raises
    ------
    TypeError
        if input is not an instance of airrtable
    """
    if not isinstance(airrtable, AirrTable):
        raise TypeError(f"{type(airrtable)} must be an instance of AirrTable")

    key = airrtable.key_column
    if airrtable.__class__ == LinkedAirrTable:
        # split table into left and right (heavy and light) tables
        left_table, right_table = airrtable.get_split_table()
        left_table = run_mutational_analysis(left_table, scheme)
        right_table = run_mutational_analysis(right_table, scheme)
        key = airrtable.key_column
        return LinkedAirrTable(left_table.merge(right_table, on=key, suffixes=airrtable.suffixes), key_column=key)

    if not airrtable.index.is_monotonic_increasing:
        raise IndexError(f"{airrtable.index} must be monotonic increasing")

    # create anarci api
    logger.info("Running ANARCI on germline alignment")
    anarci_api = Anarci(scheme=scheme, allowed_chain=["H", "K", "L"])
    germline_results_anarci = anarci_api.run_dataframe(
        airrtable["germline_alignment_aa"].str.replace("-", "").to_frame().join(airrtable[key]),
        key,
        "germline_alignment_aa",
    )
    logger.info("Running ANARCI on mature alignment")
    mature_results_anarci = anarci_api.run_dataframe(
        airrtable["sequence_alignment_aa"].str.replace("-", "").to_frame().join(airrtable[key]),
        key,
        "sequence_alignment_aa",
    )
    logger.info("Getting ANARCI on alignment tables")
    sets_of_lists = [
        set(germline_results_anarci["Id"]),
        set(mature_results_anarci["Id"]),
        set(airrtable[key].astype("str")),
    ]
    sets_of_lists = sorted(sets_of_lists, key=lambda x: len(x))
    common_results = set.intersection(*sets_of_lists)

    logger.info(f"Can run mutational analysis on {len(common_results)} out of {len(airrtable)} results")
    germline_results_anarci = germline_results_anarci.loc[germline_results_anarci["Id"].isin(common_results), :]
    germline_results_anarci_at = germline_results_anarci.get_alignment_table()
    mature_results_anarci = mature_results_anarci.loc[mature_results_anarci["Id"].isin(common_results), :]
    mature_results_anarci_at = mature_results_anarci.get_alignment_table()
    lookup_dataframe = (
        mature_results_anarci_at.drop(["chain_type", "scheme"], axis=1)
        .set_index("Id")
        .transpose()
        .join(
            germline_results_anarci_at.drop(["chain_type", "scheme"], axis=1).set_index("Id").transpose(),
            lsuffix="_mature",
            rsuffix="_germ",
        )
    )
    lookup_dataframe = lookup_dataframe[sorted(lookup_dataframe.columns)].fillna("-")
    mutation_arrays = []
    logger.info(f"Finding mutations on {len(mature_results_anarci)} sequences")
    for x in mature_results_anarci["Id"]:
        germ_tag = x + "_germ"
        mat_tag = x + "_mature"

        # get section of dataframe for only the two we are interested in
        lookup_specific = lookup_dataframe[[germ_tag, mat_tag]]

        # mutation array are all the mutations in a list
        mutation_array = lookup_specific[lookup_specific.apply(lambda x: x[0] != x[1] and x[0] != "X", axis=1)].apply(
            lambda x: x[0] + x.name + x[1], axis=1
        )
        if mutation_array.empty:
            mutation_array = []
        else:
            mutation_array = mutation_array.to_list()
        mutation_arrays.append(mutation_array)

    mature_results_anarci["mutations"] = mutation_arrays
    return AirrTable(
        airrtable.astype({key: "str"}).merge(
            mature_results_anarci.rename({"Id": key}, axis=1)[[key, "scheme", "mutations"]], on=key
        ),
        key_column=key,
    )


def _get_igl(row: pd.Series) -> str:
    """Get infered germline sequences from Airr Row

    Parameters
    ----------
    row : pd.Series
        A row from the airr table
    Returns
    -------
    str
        the igl sequecne
    """

    # get germline components
    v_germline = row.v_germline_alignment_aa
    full_germline = row.germline_alignment_aa
    if isinstance(v_germline, float):
        return
    cdr3_j_germline = full_germline[len(v_germline) :]

    # get mature components
    v_mature = row.v_sequence_alignment_aa
    full_mature = row.sequence_alignment_aa
    cdr3_j_mature = full_mature[len(v_mature) :]

    # if the mature and cdr3 are not the same size
    # this will happen on non-productive
    if len(cdr3_j_mature) != len(cdr3_j_germline):
        logger.debug(f"{row.name} - strange iGL")
        return

    iGL_cdr3 = ""
    for mature, germline in zip(cdr3_j_mature, cdr3_j_germline):
        if germline == "X" or germline == "-":
            iGL_cdr3 += mature
            continue
        iGL_cdr3 += germline

    full_igl = v_germline.replace("-", "") + iGL_cdr3.replace("-", "")
    return full_igl


def run_igl_assignment(airrtable: Union[AirrTable, LinkedAirrTable]) -> Union[AirrTable, LinkedAirrTable]:
    if not isinstance(airrtable, AirrTable):
        raise TypeError(f"{type(airrtable)} must be an instance of AirrTable")

    key = airrtable.key_column
    if airrtable.__class__ == LinkedAirrTable:
        # split table into left and right (heavy and light) tables
        left_table, right_table = airrtable.get_split_table()
        left_table = run_igl_assignment(left_table)
        right_table = run_igl_assignment(right_table)
        key = airrtable.key_column
        return LinkedAirrTable(left_table.merge(right_table, on=key, suffixes=airrtable.suffixes), key_column=key)

    airrtable["iGL"] = airrtable.apply(_get_igl, axis=1)
    return airrtable
