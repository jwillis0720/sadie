import logging
from typing import Any, List, Union

import numpy as np

# third party
import pandas as pd
from Bio.Seq import Seq
from Levenshtein import distance

# module/package level
from sadie.airr.airrtable import AirrTable, LinkedAirrTable
from sadie.renumbering import Renumbering

logger = logging.getLogger("AirrMethod")


def run_mutational_analysis(
    airrtable: Union[AirrTable, LinkedAirrTable], scheme: str, run_multiproc: bool = True
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
        left_table = run_mutational_analysis(left_table, scheme, run_multiproc)
        right_table = run_mutational_analysis(right_table, scheme, run_multiproc)
        key = airrtable.key_column

        l_suffix = airrtable.suffixes[0]
        r_suffix = airrtable.suffixes[1]

        return LinkedAirrTable(left_table.merge(right_table, on=key, suffixes=(l_suffix, r_suffix)), key_column=key)

    if not airrtable.index.is_monotonic_increasing:
        raise IndexError(f"{airrtable.index} must be monotonic increasing")

    # create Renumbering api
    logger.info("Running Renumbering on germline alignment")
    renumbering_api = Renumbering(scheme=scheme, allowed_chain=["H", "K", "L"], run_multiproc=run_multiproc)
    airrtable = pd.DataFrame(airrtable)  # type: ignore[assignment]
    germline_results_renumbering = renumbering_api.run_dataframe(
        airrtable["germline_alignment_aa"].str.replace("-", "").to_frame().join(airrtable[key]),
        key,
        "germline_alignment_aa",
    )
    logger.info("Running Renumbering on mature alignment")

    mature_results_renumbering = renumbering_api.run_dataframe(
        airrtable["sequence_alignment_aa"].str.replace("-", "").to_frame().join(airrtable[key]),
        key,
        "sequence_alignment_aa",
    )
    logger.info("Getting Renumbering on alignment tables")
    sets_of_lists = [
        set(germline_results_renumbering["Id"]),
        set(mature_results_renumbering["Id"]),
        set(airrtable[key].astype("str")),
    ]
    sets_of_lists = sorted(sets_of_lists, key=lambda x: len(x))
    common_results = set.intersection(*sets_of_lists)

    logger.info(f"Can run mutational analysis on {len(common_results)} out of {len(airrtable)} results")
    germline_results_renumbering = germline_results_renumbering.loc[
        germline_results_renumbering["Id"].isin(common_results), :  # type: ignore
    ]
    germline_results_renumbering_at = germline_results_renumbering.get_alignment_table()
    mature_results_renumbering = mature_results_renumbering.loc[  # type ignore
        mature_results_renumbering["Id"].isin(common_results), :  # type: ignore
    ]
    mature_results_renumbering_at = mature_results_renumbering.get_alignment_table()
    lookup_dataframe = (
        mature_results_renumbering_at.drop(["chain_type", "scheme"], axis=1)
        .set_index("Id")
        .transpose()
        .join(
            germline_results_renumbering_at.drop(["chain_type", "scheme"], axis=1).set_index("Id").transpose(),
            lsuffix="_mature",
            rsuffix="_germ",
        )
    )
    lookup_dataframe = lookup_dataframe[sorted(lookup_dataframe.columns)].fillna("-")
    mutation_arrays = []
    logger.info(f"Finding mutations on {len(mature_results_renumbering)} sequences")
    for x in mature_results_renumbering["Id"]:
        germ_tag = x + "_germ"
        mat_tag = x + "_mature"

        # get section of dataframe for only the two we are interested in
        lookup_specific = lookup_dataframe[[germ_tag, mat_tag]]

        # mutation array are all the mutations in a list
        mutation_array = lookup_specific[lookup_specific.apply(lambda x: x[0] != x[1] and x[0] != "X", axis=1)].apply(
            lambda x: x[0] + x.name + x[1], axis=1
        )
        if mutation_array.empty:
            mutation_array = []  # type: ignore[assignment]
        else:
            mutation_array = mutation_array.to_list()
        mutation_arrays.append(mutation_array)

    mature_results_renumbering["mutations"] = mutation_arrays
    return AirrTable(
        pd.DataFrame(airrtable)
        .astype({key: "str"})
        .merge(
            pd.DataFrame(mature_results_renumbering).rename({"Id": key}, axis=1)[[key, "scheme", "mutations"]], on=key
        ),
        key_column=key,
    )


def _get_igl_aa(row: pd.Series) -> Union[str, float]:  # type: ignore
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
        return np.nan
    cdr3_j_germline = full_germline[len(v_germline) :]

    # get mature components
    v_mature = row.v_sequence_alignment_aa
    full_mature = row.sequence_alignment_aa

    cdr3_j_mature = full_mature[len(v_mature) :]

    # if the mature and cdr3 are not the same size
    # this will happen on non-productive
    if len(cdr3_j_mature) != len(cdr3_j_germline):
        logger.debug(f"{row.name} - strange iGL")
        return np.nan

    iGL_cdr3 = ""
    for mature, germline in zip(cdr3_j_mature, cdr3_j_germline):
        if germline == "X" or germline == "-" or germline == "*":
            iGL_cdr3 += mature
            continue
        iGL_cdr3 += germline

    # remove insertions from seq
    full_igl: str = v_germline.replace("-", "") + iGL_cdr3.replace("-", "")
    if "*" in full_igl:
        raise ValueError(f"{row.name} - strange iGL has * in it")
    return full_igl


codon_table = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "I": ("ATT", "ATC", "ATA"),
    "H": ("CAT", "CAC"),
    "K": ("AAA", "AAG"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "*": ("TAA", "TAG", "TGA"),
}


def _get_3mers(string: str) -> List[str]:
    parts = [string[i : i + 3] for i in range(0, len(string), 3)]
    return parts
    # return list(filter(lambda x: len(x) == 3, parts))


def _find_best_codon(codon: str, aa: str) -> str:
    # first arg is a partial codon, second arg is an aa.
    if len(codon) == 3:
        return codon
    candidates = codon_table[aa]

    def best_dist(x: str) -> Any:
        return distance(x, codon)

    best_codon: str = sorted(candidates, key=best_dist)[0]
    return best_codon


def _get_igl_nt(row: pd.Series) -> Union[str, float]:  # type: ignore
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

    # this will be our germline iGL nt
    germline_igl = ""

    # grab the alignment that has hopefully been corrected
    germline_alignment_aa = row["germline_alignment_aa"]
    sequence_alignment_aa = row["sequence_alignment_aa"]

    # make sure to take out the indels so we can just get the codond
    germline_alignment_codons = _get_3mers(row["germline_alignment"].replace("-", ""))
    sequence_alignment_codons = _get_3mers(row["sequence_alignment"].replace("-", ""))

    if isinstance(germline_alignment_aa, float) or isinstance(sequence_alignment_aa, float):
        return np.nan

    if len(germline_alignment_aa) != len(sequence_alignment_aa):
        raise ValueError(
            f"{row.index} - germline aa alignment is not the same length as the sequence aa alignment {len(germline_alignment_aa)} != {len(sequence_alignment_aa)}"
        )

    # indexer
    germline_index, sequence_index = 0, 0

    # iterate over bot alignments
    for (germline_aa, sequence_aa) in zip(germline_alignment_aa, sequence_alignment_aa):

        # we can't have both - in germ and sequence
        if germline_aa == "-" and sequence_aa == "-":
            raise Exception("this is not possible")

        # if germline == sequence, take the germline codon
        if germline_aa == sequence_aa:
            codon = sequence_alignment_codons[sequence_index]
            germline_igl += codon
            germline_index += 1
            sequence_index += 1

        # if germline is -, only increment sequence but not take a "codon"
        elif germline_aa == "-":
            if germline_index == len(germline_alignment_aa) - 1:
                # we are at the end so pad it with sequence
                # print('here',germline_igl)
                partial_codon = sequence_alignment_codons[sequence_index]
                best_codon = _find_best_codon(partial_codon, sequence_aa)
                logger.debug(f"Partial codon:{partial_codon} for mature {sequence_aa} choosing {best_codon}")
                germline_igl += best_codon

            sequence_index += 1

        # if sequence is - (a deletion), increment the germline codon and take it
        elif sequence_aa == "-":
            germline_igl += germline_alignment_codons[germline_index]
            germline_index += 1
        else:
            # finally, deal with the cases where the germline and sequence do not equal each other but are not idnesl
            # if gerline is X or *, take mature codon
            if germline_aa == "X" or germline_aa == "*":
                germline_igl += sequence_alignment_codons[sequence_index]

            # else take germline
            else:
                codon = germline_alignment_codons[germline_index]
                germline_igl += codon
            sequence_index += 1
            germline_index += 1
    real_igl = row["iGL_aa"]
    contrived_igl = str(Seq(germline_igl).translate())
    if real_igl != contrived_igl:
        if real_igl[:-1] == contrived_igl:
            ending_codon = _find_best_codon(sequence_alignment_codons[-1], real_igl[-1])
            germline_igl = germline_igl[: -(len(germline_igl) % 3)] + ending_codon
            # we are trying to fix it
            if real_igl != str(Seq(germline_igl).translate()):
                if row["productive"]:
                    raise ValueError(
                        f"{row.name} - iGL_aa {row['iGL_aa']} != contrived {Seq(germline_igl).translate()}"
                    )
        else:
            if row["productive"]:
                raise ValueError(f"{row.name} - iGL_aa {row['iGL_aa']} != contrived {Seq(germline_igl).translate()}")
    return germline_igl


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
        l_suffix = airrtable.suffixes[0]
        r_suffix = airrtable.suffixes[1]
        return LinkedAirrTable(left_table.merge(right_table, on=key, suffixes=(l_suffix, r_suffix)), key_column=key)

    airrtable["iGL_aa"] = airrtable.apply(_get_igl_aa, axis=1)
    airrtable["iGL"] = airrtable.apply(_get_igl_nt, axis=1)
    return airrtable
