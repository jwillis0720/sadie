from Levenshtein import distance
import pandas as pd


def compute_lev(x, y):
    try:
        return distance(x, y)
    except TypeError:
        # print("problem with retu ",x)
        return len(y)


def compute_lev_norm(x, y):
    try:
        return y / len(x)
    except TypeError:
        return 1


def compute_lev_table_from_airr(airr, target_antibody):
    lev_df = pd.DataFrame()
    # FW1
    lev_df["fwr1_aa_lev"] = airr["fwr1_aa"].apply(lambda x: compute_lev(x, target_antibody["fwr1_aa"].get_aa()))

    # FW1 Norm
    lev_df["fwr1_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["fwr1_aa"], lev_df["fwr1_aa_lev"]))

    # CDR1
    lev_df["cdr1_aa_lev"] = airr["cdr1_aa"].apply(lambda x: compute_lev(x, target_antibody["cdr1_aa"].get_aa()))

    # CDR1 Norm
    lev_df["cdr1_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["cdr1_aa"], lev_df["cdr1_aa_lev"]))

    # FW2
    lev_df["fwr2_aa_lev"] = airr["fwr2_aa"].apply(lambda x: compute_lev(x, target_antibody["fwr2_aa"].get_aa()))

    # FW2 Norm
    lev_df["fwr2_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["fwr2_aa"], lev_df["fwr2_aa_lev"]))

    # CDR2
    lev_df["cdr2_aa_lev"] = airr["cdr2_aa"].apply(lambda x: compute_lev(x, target_antibody["cdr2_aa"].get_aa()))

    # CDR2 Norm
    lev_df["cdr2_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["cdr2_aa"], lev_df["cdr2_aa_lev"]))

    # FW3
    lev_df["fwr3_aa_lev"] = airr["fwr3_aa"].apply(lambda x: compute_lev(x, target_antibody["fwr3_aa"].get_aa()))

    # FW3 Norm
    lev_df["fwr3_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["fwr3_aa"], lev_df["fwr3_aa_lev"]))

    # CDR3
    lev_df["cdr3_aa_lev"] = airr["cdr3_aa"].apply(lambda x: compute_lev(x, target_antibody["cdr3_aa"].get_aa()))

    # CDR3 Norm
    lev_df["cdr3_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["cdr3_aa"], lev_df["cdr3_aa_lev"]))

    # FW4
    lev_df["fwr4_aa_lev"] = airr["fwr4_aa"].apply(lambda x: compute_lev(x, target_antibody["fwr4_aa"].get_aa()))

    # CDR4 Norm
    lev_df["fwr4_aa_lev_norm"] = list(map(lambda x, y: y / len(x), airr["fwr4_aa"], lev_df["fwr4_aa_lev"]))

    lev_df["vh_fwr_lev_sum"] = lev_df["fwr1_aa_lev"] + lev_df["fwr2_aa_lev"] + lev_df["fwr3_aa_lev"]
    lev_df["vh_fwr_lev_norm"] = lev_df["fwr1_aa_lev_norm"] + lev_df["fwr2_aa_lev_norm"] + lev_df["fwr3_aa_lev_norm"]
    lev_df["vh_cdr_lev_sum"] = lev_df["cdr1_aa_lev"] + lev_df["cdr2_aa_lev"]
    lev_df["vh_cdr_lev_norm"] = lev_df["cdr1_aa_lev_norm"] + lev_df["cdr2_aa_lev_norm"]
    lev_df["vh_cdr_lev_sum_w_cdr3"] = lev_df["cdr1_aa_lev"] + lev_df["cdr2_aa_lev"] + lev_df["cdr3_aa_lev"]
    lev_df["vh_cdr_lev_norm_w_cdr3"] = (
        lev_df["cdr1_aa_lev_norm"] + lev_df["cdr2_aa_lev_norm"] + lev_df["cdr3_aa_lev_norm"]
    )
    lev_df["vh_sum"] = lev_df["vh_fwr_lev_sum"] + lev_df["vh_cdr_lev_sum"]
    lev_df["vh_norm"] = lev_df["vh_fwr_lev_norm"] + lev_df["vh_cdr_lev_norm"]

    # And set index
    lev_df.index = airr.index
    return lev_df
