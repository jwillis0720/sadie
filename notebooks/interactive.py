# %%
import sadie
import pandas as pd

pd.set_option("display.max_rows", 45)
pd.set_option("display.max_columns", 500)
pd.set_option("display.width", 100000)
airr_api = sadie.airr.Airr("human")
results = airr_api.run_file("../tests/integration/airr/fixtures/catnap_nt_heavy.fasta.gz")
#%%
# %%
normal = airr_api.run_single(
    "normal",
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
)
# %%
deletion = airr_api.run_single(
    "deletion",
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
)
deletion

# %%
airr_api.igblast.gap_open = 2
insertion = airr_api.run_single(
    "insertion",
    "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTGAAGCCTGGAGAGGGATCCCTTAGACTCTCATGTTCAGCCTCTGGTTTCGACTTCGATAACGCCTGGATGACTTGGGTCCGCCAGCCTCCAGGGAAGGGCCTCGAATGGGTTGGTCGTATTACGGGTCCAGGTGAAGGTTGGTCAGTGGACTATGCTGCACCCGTGGAAGGCAGATTTACCATCTCGAGACTCAATTCAATAAATTTCTTATATTTGGAGATGAACAATTTAAGAATGGAAGACTCAGGCCTTTACTTCTGTGCCCGCACGGGAAAATATTATGATTTTTGGAGTGGCTATCCGCCGGGAGAAGAATACTTCCAAGACTGGGGCCGGGGCACCCTGGTCACCGTCTCCTCA",
)
insertion

# %%
correctioon = pd.concat(
    [normal.table, deletion.table, insertion.table]
)  # .to_csv("~/DropboxPersonal/Dropbox/look_at_me.csv")
# %%


def correct_sequence_alignments(x):
    sequence_alignment = x
    print(sequence_alignment)
    for x in range(0, len(sequence_alignment), 3):
        codon = sequence_alignment[x : x + 3]


results.table.head(1)["sequence_alignment"].apply(correct_sequence_alignments)


# %%
# %%

# %%
# %%
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from gspread_pandas import Spread

# gspread = Spread("SadieDebug", create_spread=True)
matrix = matlist.blosum62
suspect_results = results[results["sequence_alignment_aa"].str.len() != results["germline_alignment_aa"].str.len()]


def get_alignments(x, y):
    alignments = pairwise2.align.globalds(
        x, y, matrix, -12, -4, penalize_extend_when_opening=True, penalize_end_gaps=False
    )
    corrected_x = alignments[0][0]
    corrected_y = alignments[0][1]
    return pd.Series({"sequence_alignment_aa": corrected_x, "germline_alignment_aa": corrected_y})


d = suspect_results.apply(lambda x: get_alignments(x["sequence_alignment_aa"], x["germline_alignment_aa"]), axis=1)
d
# df = suspect_results.join(d)
# for col in df.dtypes[
#     df.dtypes.apply(lambda x: x in ["float64", "float16", "int16"] or isinstance(x, pd.Int16Dtype))
# ].index:
#     # print(col,df.dtypes[col])
#     df[col] = df[col].astype("float")

# gspread.df_to_sheet(df, sheet="corrected")
# d.iloc[0]["y"]


# %%
suspect_results[
    [
        "sequence_id",
        "sequence",
        "sequence_alignment",
        "germline_alignment",
        "sequence_alignment_aa",
        "germline_alignment_aa",
    ]
].sample(5).to_markdown()
# %%

# %%
