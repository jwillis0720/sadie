"""[summary]"""
import logging
import re
import tempfile

import pandas as pd
from sadie.anarci import Anarci
from Bio import Seq

from .imgt import assign_index_position, multiindex_pivot
from .settings import MOTIF_LOOKUP, RENAME_DICT_TRANSLATE

logger = logging.getLogger(__name__)
# from gspread_pandas import Spread, Client

# spread = Spread("Germline Reference", create_spread=True)


def make_segments(row):
    vdj_nt = row["v_gene_nt"]
    fw1_aa = row["fwr1_aa"]
    fw1_aa_end = len(fw1_aa)
    fw1_nt_start = 0
    fw1_nt_end = fw1_aa_end * 3
    fw1_nt = vdj_nt[fw1_nt_start:fw1_nt_end]
    assert str(Seq.Seq(fw1_nt).translate()) == fw1_aa

    cdr1_aa = row["cdr1_aa"]
    cdr1_aa_end = len(cdr1_aa)
    cdr1_nt_start = fw1_nt_end
    cdr1_nt_end = fw1_nt_end + (cdr1_aa_end * 3)
    cdr1_nt = vdj_nt[cdr1_nt_start:cdr1_nt_end]
    assert str(Seq.Seq(cdr1_nt).translate()) == cdr1_aa

    fw2_aa = row["fwr2_aa"]
    fw2_aa_end = len(fw2_aa)
    fw2_nt_start = cdr1_nt_end
    fw2_nt_end = cdr1_nt_end + (fw2_aa_end * 3)
    fw2_nt = vdj_nt[fw2_nt_start:fw2_nt_end]
    assert str(Seq.Seq(fw2_nt).translate()) == fw2_aa

    cdr2_aa = row["cdr2_aa"]
    cdr2_aa_end = len(cdr2_aa)
    cdr2_nt_start = fw2_nt_end
    cdr2_nt_end = fw2_nt_end + (cdr2_aa_end * 3)
    cdr2_nt = vdj_nt[cdr2_nt_start:cdr2_nt_end]
    assert str(Seq.Seq(cdr2_nt).translate()) == cdr2_aa

    fw3_aa = row["fwr3_aa"]
    fw3_aa_end = len(fw3_aa)
    fw3_nt_start = cdr2_nt_end
    fw3_nt_end = cdr2_nt_end + (fw3_aa_end * 3)
    fw3_nt = vdj_nt[fw3_nt_start:fw3_nt_end]
    assert str(Seq.Seq(fw3_nt).translate()) == fw3_aa

    cdr3_nt = vdj_nt[fw3_nt_end:]

    return pd.Series(
        {
            "fwr1_nt": fw1_nt,
            "cdr1_nt": cdr1_nt,
            "fwr2_nt": fw2_nt,
            "cdr2_nt": cdr2_nt,
            "fwr3_nt": fw3_nt,
            "cdr3_nt": cdr3_nt,
        }
    )


def generate_v_segment_data(imgt_db_engine):
    """[summary]

    Args:
        imgt_db_path ([type]): sqlalchemy connection

    Returns:
        [type]: [description]
    """
    imgt_db_unique = pd.read_sql("imgt_unique", con=imgt_db_engine, index_col="index")
    custom = pd.read_sql("custom", con=imgt_db_engine, index_col="index")
    v_segment_only_clean_imgt = imgt_db_unique.loc[
        (imgt_db_unique["gene_segment"] == "V")
        & (imgt_db_unique["functional"] == "F")
        & (imgt_db_unique["label"] != "v_gene_nt"),
        :,
    ]

    v_fit_index = (v_segment_only_clean_imgt["nt_sequence_no_gaps"].str.len() % 3 == 0) | (
        v_segment_only_clean_imgt["label"] == "cdr3_nt"
    )
    v_segment_only_clean_imgt.loc[v_fit_index, "aa_sequence_no_gaps"] = v_segment_only_clean_imgt.loc[
        v_fit_index, "nt_sequence_no_gaps"
    ].apply(lambda x: Seq.Seq(x).translate().__str__())

    # Now do the same for custom
    v_segment_only_clean_custom = custom.loc[(custom["gene_segment"] == "V") & (custom["label"] != "v_gene_nt"), :]
    v_fit_index = (v_segment_only_clean_custom["nt_sequence_no_gaps"].str.len() % 3 == 0) | (
        v_segment_only_clean_custom["label"] == "cdr3_nt"
    )
    v_segment_only_clean_custom.loc[v_fit_index, "aa_sequence_no_gaps"] = v_segment_only_clean_custom.loc[
        v_fit_index, "nt_sequence_no_gaps"
    ].apply(lambda x: Seq.Seq(x).translate().__str__())

    v_segment_combined = pd.concat([v_segment_only_clean_custom, v_segment_only_clean_imgt]).reset_index(drop=True)
    v_segment_combined_pivot_nt = multiindex_pivot(
        v_segment_combined.set_index(["gene", "latin", "common", "functional", "gene_segment", "source"]),
        columns="label",
        values="nt_sequence_no_gaps",
    )
    v_segment_only_clean_pivot_aa = multiindex_pivot(
        v_segment_combined.set_index(["gene", "latin", "common", "functional", "gene_segment", "source"]),
        columns="label",
        values="aa_sequence_no_gaps",
    ).rename(RENAME_DICT_TRANSLATE, axis=1)
    v_segment_joined = v_segment_combined_pivot_nt.join(v_segment_only_clean_pivot_aa).fillna("")
    dropped_v_genes = []

    for key in [
        "fwr1_nt",
        "cdr1_nt",
        "fwr2_nt",
        "cdr2_nt",
        "fwr3_nt",
        "fwr1_aa",
        "cdr1_aa",
        "fwr2_aa",
        "cdr2_aa",
        "fwr3_aa",
    ]:
        dropped_index = v_segment_joined[v_segment_joined[key] == ""].index
        logger.debug(f"Dropping {len(dropped_index)} indexes because of missing {key} region - {list(dropped_index)}")
        dropped_v_genes.append(v_segment_joined.loc[dropped_index])
        v_segment_joined.drop(v_segment_joined[v_segment_joined[key] == ""].index, inplace=True)
    dropped_v_genes = pd.concat(dropped_v_genes)
    preferred_order = [
        "common",
        "latin",
        "gene",
        "functional",
        "gene_segment",
        "source",
        "fwr1_nt",
        "cdr1_nt",
        "fwr2_nt",
        "cdr2_nt",
        "fwr3_nt",
        "cdr3_nt",
        "fwr1_aa",
        "cdr1_aa",
        "fwr2_aa",
        "cdr2_aa",
        "fwr3_aa",
        "cdr3_aa",
    ]
    v_segment_joined = (
        v_segment_joined.reset_index()[preferred_order]
        .sort_values(["common", "latin", "gene", "source"])
        .reset_index(drop=True)
    )
    dropped_v_genes = (
        dropped_v_genes.reset_index()[preferred_order]
        .sort_values(["common", "latin", "gene", "source"])
        .reset_index(drop=True)
    )

    concat = pd.concat([imgt_db_unique, custom]).reset_index()
    # Add all v gene nt as defined by imgt
    v_segment_joined = v_segment_joined.merge(
        concat[concat["label"] == "v_gene_nt"][
            [
                "common",
                "latin",
                "gene",
                "source",
                "nt_sequence_no_gaps",
                "nt_sequence_gaps",
            ]
        ].rename(
            {
                "nt_sequence_no_gaps": "v_gene_nt_imgt",
                "nt_sequence_gaps": "v_gene_nt_imgt_gaps",
            },
            axis=1,
        ),
        how="left",
    )

    # but lets keep track of the dropped genes
    dropped_v_genes = dropped_v_genes.merge(
        concat[concat["label"] == "v_gene_nt"][
            [
                "common",
                "latin",
                "gene",
                "source",
                "nt_sequence_no_gaps",
                "nt_sequence_gaps",
            ]
        ].rename(
            {
                "nt_sequence_no_gaps": "v_gene_nt_imgt",
                "nt_sequence_gaps": "v_gene_nt_imgt_gaps",
            },
            axis=1,
        ),
        how="left",
    )
    # take our segmented definition and add it backtogether
    v_segment_joined.loc[:, "v_gene_nt"] = v_segment_joined[
        ["fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt"]
    ].apply(lambda x: "".join(x), axis=1)
    dropped_v_genes.loc[:, "v_gene_nt"] = dropped_v_genes[
        ["fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt"]
    ].apply(lambda x: "".join(x), axis=1)
    v_segment_joined.loc[:, "v_gene_aa"] = v_segment_joined[
        ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa"]
    ].apply(lambda x: "".join(x), axis=1)
    dropped_v_genes.loc[:, "v_gene_aa"] = dropped_v_genes[
        ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa"]
    ].apply(lambda x: "".join(x), axis=1)

    # add a flag to see if if imgt v gene nt is the same as our piecemeal definition
    v_segment_joined.loc[:, "v_gene_imgt_is_trans"] = list(
        v_segment_joined["v_gene_nt"] == v_segment_joined["v_gene_nt_imgt"]
    )
    dropped_v_genes.loc[:, "v_gene_imgt_is_trans"] = list(
        dropped_v_genes["v_gene_nt"] == dropped_v_genes["v_gene_nt_imgt"]
    )
    v_segment_joined.loc[:, "region_definiton"] = "imgt"
    v_segment_joined.loc[:, "assigned_by"] = "imgt"

    # This will break down the v segment into its positions in teh sequence
    v_segment_joined = assign_index_position(v_segment_joined, molecule="aa")
    v_segment_joined = assign_index_position(v_segment_joined, molecule="nt")
    dropped_v_genes = assign_index_position(dropped_v_genes, molecule="aa")
    dropped_v_genes = assign_index_position(dropped_v_genes, molecule="nt")
    logger.info(f"Writing IMGT V segment to {imgt_db_engine}")
    v_segment_joined.to_sql("v_segment_imgt", con=imgt_db_engine, if_exists="replace")
    dropped_v_genes.to_sql("v_segment_imgt_dropped", con=imgt_db_engine, if_exists="replace")
    with tempfile.NamedTemporaryFile("w", dir=".", suffix=".fasta") as tempfile_fasta:
        for group, group_df in v_segment_joined.groupby(["latin", "common", "source"]):
            _group_df = group_df
            # for each sequence lets add some x's and add it to a temperofy file
            for sequence_name, sequence in zip(_group_df["gene"], _group_df["v_gene_aa"] + "X" * 14):
                tempfile_fasta.write(f">{group[0]}|{group[1]}|{group[2]}|{sequence_name}\n{sequence}\n")

        for definition in ["scdr", "kabat", "chothia", "abm", "contact"]:
            logger.info(f"Resegmenting database with {definition} boundaries")
            # Taking imgt definition and resegmenting it as a new definitnon
            anarci_api = Anarci(region_assign=definition)
            assignments = anarci_api.run_file(tempfile_fasta.name, multi=True).segment_table_no_gaps.reset_index()
            assignments = (
                assignments["Id"]
                .str.split("|")
                .apply(lambda x: pd.Series({"latin": x[0], "common": x[1], "gene": x[-1], "source": x[2]}))
                .join(assignments)
            )
            assignments_clean = assignments[assignments["leader"] == ""]

            # sometimes the cdr3 portion gets assigned as a tail
            assignments_clean.loc[:, "cdr3_aa"] = assignments_clean["cdr3_aa"] + assignments_clean["tail"]

            # remove the x we used as a buffer
            assignments_clean.loc[:, "cdr3_aa"] = assignments_clean["cdr3_aa"].str.rstrip("X")
            assignments_clean.loc[:, "fwr3_aa"] = assignments_clean["fwr3_aa"].str.rstrip("X")
            assignments_clean.loc[:, "vdj"] = assignments_clean["vdj"].str.rstrip("X")
            assignments_clean = assignments_clean[
                [
                    "common",
                    "latin",
                    "gene",
                    "source",
                    "fwr1_aa",
                    "cdr1_aa",
                    "fwr2_aa",
                    "cdr2_aa",
                    "fwr3_aa",
                    "cdr3_aa",
                ]
            ]
            assignments_clean.loc[:, "v_gene_aa"] = assignments_clean[
                ["fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa"]
            ].apply(lambda x: "".join(x), axis=1)
            assignments_clean = assignments_clean.merge(
                v_segment_joined[
                    [
                        "common",
                        "latin",
                        "gene",
                        "source",
                        "functional",
                        "v_gene_nt_imgt",
                        "v_gene_nt",
                        "v_gene_aa",
                    ]
                ],
                on=["common", "latin", "gene", "source", "v_gene_aa"],
                how="inner",
            )
            translations = assignments_clean.apply(make_segments, axis=1)
            assignments_clean = assignments_clean.join(translations)
            assignments_clean.loc[:, "assigned_by"] = "anarci"
            assignments_clean.loc[:, "region_definiton"] = definition
            assignments = assign_index_position(assignments_clean, molecule="aa")
            assignments = assign_index_position(assignments_clean, molecule="nt")
            logger.info(f"Writing {definition} V segment to {imgt_db_engine}")
            assignments.to_sql(f"v_segment_{definition}", con=imgt_db_engine, if_exists="replace")


def translate_j_segment(X):
    nt_seq = X[0]
    reading_frame = X[1].strip()
    if reading_frame:
        aa_seq = Seq.Seq(nt_seq[int(reading_frame) - 1 :]).translate().__str__()
    else:
        longest = 0
        reading_frame = 0
        for potential_reading_frame in [1, 2, 3]:
            aa_seq = Seq.Seq(nt_seq[int(potential_reading_frame) - 1 :]).translate().__str__()
            if "*" in aa_seq:

                stop_index = aa_seq.index("*")
            else:
                stop_index = len(aa_seq)
            if stop_index > longest:
                reading_frame = potential_reading_frame
                longest = stop_index

            # if nt_seq == "GACTATAAAGGAGGTGACACCCTTTTGACCGTGAAA":
            # print("here")
        aa_seq = Seq.Seq(nt_seq[int(reading_frame) - 1 :]).translate().__str__()
    return pd.Series({"aa_sequence_no_gaps": aa_seq, "reading_frame": reading_frame})


def split_fw4(X):
    common = X[0]
    aa_seq = X[1]
    gene = X[2]
    reading_frame = X[3]
    nt = X[4]
    short_name = gene.split("*")[0]
    gene_type = gene[:4]
    return_series = pd.Series(
        {
            "cdr3_nt": "",
            "frw4_nt": "",
            "cdr3_aa": "",
            "frw4_aa": "",
            "end_cdr3_nt_index": "",
        }
    )

    if common in MOTIF_LOOKUP.keys():
        if gene_type in MOTIF_LOOKUP[common].keys():
            motif = MOTIF_LOOKUP[common][gene_type]
        else:
            logger.warning(
                f"{gene_type}-{short_name} for {common} not defined, but other genes are for this species..."
            )
            return return_series
    else:
        return return_series
    expression = re.findall(motif, aa_seq)
    ignore = MOTIF_LOOKUP[common]["ignore"]
    if not expression:
        if short_name in ignore or aa_seq in ignore:
            logger.info(f"Settign file says to skip {common}-{short_name}-{aa_seq}")
            return pd.Series(
                {
                    "cdr3_nt": "",
                    "frw4_nt": "",
                    "cdr3_aa": "",
                    "frw4_aa": "",
                    "end_cdr3_nt_index": "",
                }
            )
        else:
            # if common == "sharks":
            #     print(f"{aa_seq}")
            #     return return_series
            raise Exception(
                f"""
                    Cant find a matching motif for {common}-{short_name} with regex {motif} in sequence {aa_seq}
                    Should this sequence be ignored? Add it to ignored under settigns.motifs file
                    """
            )

    end_motif = expression[0]
    end_cdr3_index = aa_seq.index(end_motif)
    end_index_codon = end_cdr3_index * 3 + int(reading_frame) - 1
    return pd.Series(
        {
            "cdr3_nt": nt[:end_index_codon],
            "frw4_nt": nt[end_index_codon:],
            "cdr3_aa": aa_seq[:end_cdr3_index],
            "frw4_aa": aa_seq[end_cdr3_index:],
            "end_cdr3_nt_index": end_index_codon,
        }
    )


def generate_j_gene_data(engine):
    """[summary]

    Args:
        blast_path ([type]): [description]
        aux_path ([type]): [description]

    Returns:
        [type]: [description]
    """

    imgt_all_df = pd.read_sql("imgt_all", con=engine, index_col="index")
    j_region = imgt_all_df.loc[imgt_all_df["label"] == "J-REGION"]
    j_region.loc[:, "reading_frame"] = j_region["fasta_header"].str.split("|").str.get(7)
    # show(j_region, settings={"block": True})
    chain_dfs_append = []
    for _, species_df in j_region.groupby(["common"]):
        species_df = species_df.drop(["reading_frame", "aa_sequence_no_gaps", "aa_sequence_gaps"], axis=1).join(
            species_df[["nt_sequence_no_gaps", "reading_frame"]].apply(translate_j_segment, axis=1)
        )
        for _, chain_df in species_df.groupby(species_df["gene"].str[0:4]):
            split_fw4_segments = chain_df.loc[
                chain_df.index,
                [
                    "common",
                    "aa_sequence_no_gaps",
                    "gene",
                    "reading_frame",
                    "nt_sequence_no_gaps",
                ],
            ].apply(split_fw4, axis=1)
            chain_dfs_append.append(chain_df.join(split_fw4_segments, how="outer"))
            # show(j_region, settings={"block": True})
            # j_region.update(split_fw4_segments)

    j_segment_df = pd.concat(chain_dfs_append)
    # spread.df_to_sheet(j_segment_df, sheet="j_segment_imgt")
    j_segment_df.to_sql("j_segment_imgt", con=engine, if_exists="replace")
