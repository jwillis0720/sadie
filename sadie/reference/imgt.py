# Std libs
import glob
import logging
import os
from io import StringIO

# third party
import pandas as pd
import requests
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Data import CodonTable
from bs4 import BeautifulSoup

# from gspread_pandas import Spread, Client

# module level
from .settings import (
    BLAST_CONVENTION,
    IMGT_GB_LOOKUP,
    REVERSE_IMGT_LOOKUP,
    IMGT_DEF_nt,
    RENAME_DICT,
)

base = "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
logger = logging.getLogger(__name__)
codon_table = CodonTable.standard_dna_table.forward_table
codon_table["TAG"] = "*"
codon_table["TAA"] = "*"
codon_table["TGA"] = "*"


def flatten(to_flatten):
    """Flatten list of list"""
    return [item for sublist in to_flatten for item in sublist]


def get_imgt_file_structure(path):
    """Get IMGT file structure as a dataframe"""
    logging.basicConfig(level=10)
    imgt_file_structure = []
    for (dirpath, _, filenames) in os.walk(path):
        if filenames:
            if not filenames:
                continue
            for file in filenames:
                if file.split(".")[-1] != "fasta":
                    continue
                imgt_file_structure.append(
                    {
                        "path": os.path.join(dirpath, file),
                        "file": file,
                        "rootname": file.split(".")[0],
                        "gene_segment": file.split(".")[0][-1],
                        "receptor": dirpath.split("/")[-1],
                        "latin": dirpath.split("/")[-2],
                        "species": REVERSE_IMGT_LOOKUP[dirpath.split("/")[-2]],
                        "blast_conventirion_receptor": BLAST_CONVENTION[dirpath.split("/")[-1]],
                    }
                )
    # This has the exact IMGT file path from the website but local
    # I don't want to keep taking it from their since they will ban my IP
    imgt_file_structure_df = pd.DataFrame(imgt_file_structure)
    return imgt_file_structure_df


def get_imgt_file(species, receptor, segment):
    """Get final IMGT endpoint"""
    imgt_endpoint = "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/{}/{}/{}".format(
        species, receptor, segment
    )
    return imgt_endpoint


def get_imgt_files(species, receptor):
    """List all files within a species receptor structure"""
    files_endpoint = base + "/" + species + "/" + receptor
    response = requests.get(files_endpoint)
    bs = BeautifulSoup(response.text, features="lxml")
    tb = bs.find_all("table")[0]
    imgt_files = list(pd.read_html(str(tb))[0].fillna("").iloc[2:-1]["Name"])
    return imgt_files


def get_imgt_species():
    """get all species from imgt.org/V-QUEST/IMGT_V-Quest_reference_directory/"""
    species_endpoint = base
    response = requests.get(species_endpoint)
    if response.status_code == 403:
        logger.critical("Forbbiden access to {}".format(species_endpoint))
        raise Exception("Forbidden Access, they banned your IP!!!")
    bs = BeautifulSoup(response.text, features="lxml")
    tb = bs.find_all("table")[0]
    available_receptors = list(pd.read_html(str(tb))[0].fillna("").iloc[2:-1]["Name"].str[:-1])
    return available_receptors


def get_receptor_for_species(species):
    """get the receptors given the speceis from imgt.org/V-QUEST/IMGT_V-Quest_reference_directory/{species}/"""
    receptor_endpoint = base + "/" + species
    logger.debug(receptor_endpoint)
    response = requests.get(receptor_endpoint)
    bs = BeautifulSoup(response.text, features="lxml")
    tb = bs.find_all("table")[0]
    available_species = list(pd.read_html(str(tb))[0].fillna("").iloc[2:-1]["Name"].str[:-1])
    return available_species


def write_out_fasta(file_path, endpoint):
    """Write out IMGT by endpoint to file"""
    logger.debug("Writing {} to {}".format(file_path, endpoint))

    # directory of file path
    basepath = os.path.dirname(file_path)

    # If exists warn or create
    if not os.path.exists(basepath):
        logger.info("Creating path {}".format(basepath))
        os.makedirs(basepath)
    if os.path.exists(file_path):
        logger.warning("{} exists!! overwriting".format(file_path))

    response = requests.get(endpoint)
    string_buffer = StringIO(response.text)

    # Parse the fasta file
    imgt_fasta_of_seqs = list(SeqIO.parse(string_buffer, "fasta"))
    logger.debug("Writing fasta to {}".format(file_path))
    # SeqIO write
    SeqIO.write(imgt_fasta_of_seqs, open(file_path, "w"), "fasta")


def write_imgt_local(directory):
    """
    Traverses IMGT file Structure and writes it to a directory

    direcotry = file path
    """
    for species in get_imgt_species():
        for receptor in get_receptor_for_species(species):
            for file in get_imgt_files(species, receptor):
                endpoint = get_imgt_file(species, receptor, file)
                FILE_PATH = os.path.join(directory, species, receptor, file)
                write_out_fasta(FILE_PATH, endpoint)


def multiindex_pivot(df, columns=None, values=None):
    # https://github.com/pandas-dev/pandas/issues/23955
    names = list(df.index.names)
    df = df.reset_index()
    list_index = df[names].values
    tuples_index = [tuple(i) for i in list_index]  # hashable
    df = df.assign(tuples_index=tuples_index)
    df = df.pivot(index="tuples_index", columns=columns, values=values)
    tuples_index = df.index  # reduced
    index = pd.MultiIndex.from_tuples(tuples_index, names=names)
    df.index = index
    return df


def show_pairwise(seq1, seq2):
    s1 = Seq(seq1)
    s2 = Seq(seq2)
    alignments = pairwise2.align.globalxs(s1, s2, -10, -10)
    return pairwise2.format_alignment(*alignments[0])


def apply_index_per_row(X, molecule="nt"):
    """Apply method to get ranges and indexes

    Parameters
    ----------
    X : [type]
        [description]
    molecule : str, optional
        [description], by default "nt"

    Returns
    -------
        pd.Series
    """
    v_gene = X[f"v_gene_{molecule}"]
    fw1 = X[f"fwr1_{molecule}"]
    cdr1 = X[f"cdr1_{molecule}"]
    fw2 = X[f"fwr2_{molecule}"]
    cdr2 = X[f"cdr2_{molecule}"]
    fw3 = X[f"fwr3_{molecule}"]

    fw1_start = 0
    fw1_end = len(fw1)
    cdr1_start = len(fw1)
    cdr1_end = cdr1_start + len(cdr1)
    fw2_start = cdr1_end
    fw2_end = fw2_start + len(fw2)
    cdr2_start = fw2_end
    cdr2_end = cdr2_start + len(cdr2)
    fw3_start = cdr2_end
    fw3_end = fw3_start + len(fw3)
    assert fw1 == v_gene[fw1_start:fw1_end]
    assert fw2 == v_gene[fw2_start:fw2_end]
    assert fw3 == v_gene[fw3_start:fw3_end]
    assert cdr1 == v_gene[cdr1_start:cdr1_end]
    assert cdr2 == v_gene[cdr2_start:cdr2_end]
    return pd.Series(
        {
            f"fwr1_{molecule}_range_start": fw1_start,
            f"fwr1_{molecule}_range_end": fw1_end,
            f"cdr1_{molecule}_range_start": cdr1_start,
            f"cdr1_{molecule}_range_end": cdr1_end,
            f"fwr2_{molecule}_range_start": fw2_start,
            f"fwr2_{molecule}_range_end": fw2_end,
            f"cdr2_{molecule}_range_start": cdr2_start,
            f"cdr2_{molecule}_range_end": cdr2_end,
            f"fwr3_{molecule}_range_start": fw3_start,
            f"fwr3_{molecule}_range_end": fw3_end,
            f"fwr1_{molecule}_index_start": fw1_start + 1,
            f"fwr1_{molecule}_index_end": fw1_end,
            f"cdr1_{molecule}_index_start": cdr1_start + 1,
            f"cdr1_{molecule}_index_end": cdr1_end,
            f"fwr2_{molecule}_index_start": fw2_start + 1,
            f"fwr2_{molecule}_index_end": fw2_end,
            f"cdr2_{molecule}_index_start": cdr2_start + 1,
            f"cdr2_{molecule}_index_end": cdr2_end,
            f"fwr3_{molecule}_index_start": fw3_start + 1,
            f"fwr3_{molecule}_index_end": fw3_end,
        }
    )


def assign_index_position(df, molecule="nt"):
    """Take in dataframe and apply nucleotide and aa index positionss

    Parameters
    ----------
    df : [type]
        [description]
    molecule : str, optional
        [description], by default "nt"

    Returns
    -------
        pd.DataFrame
    """
    df = df.fillna("")
    df_clean = df[
        (
            df[f"v_gene_{molecule}"]
            == df[f"fwr1_{molecule}"]
            + df[f"cdr1_{molecule}"]
            + df[f"fwr2_{molecule}"]
            + df[f"cdr2_{molecule}"]
            + df[f"fwr3_{molecule}"]
            + df[f"cdr3_{molecule}"]
        )
    ]
    if len(df_clean) < len(df):
        print(
            "incorrect trans",
            list(df[~df.index.isin(df_clean.index)][["gene", "common"]].to_numpy()[:, 0:2]),
        )
    segment_indexes = df_clean.apply(lambda x: apply_index_per_row(x, molecule), axis=1)
    return df.join(segment_indexes, how="left")


def parse_imgt_gb(gb_path):
    """Parse the IMGT GB fasta files and make a database"""
    # these are feature strings that are preloaded into features straight from IMGT

    # Where we will load IMGT GB Data
    gb_dataframe = []

    # First thing is to go through all the fasta sequence. Each entry corresponds to a feature of a certain gene
    fasta_files = glob.glob(gb_path + "/*/*.fasta")
    logger.debug("Parsing %s", fasta_files)

    for gb_file in glob.glob(gb_path + "/*/*.fasta"):

        # we have to keep track of a the seen list becaue IMGT has double entries for mice and some other species
        seen = []
        for sequence in SeqIO.parse(gb_file, "fasta"):
            parse_desc = sequence.description.split("|")
            """
            parse desc Description looks like this. We can grab everythign we need.
            ex.
            ['IMGT000050', 'IGKV(II)-8*01', 'Felis catus_Abyssinian', 'P', 'DONOR-SPLICE',
                                 '70502..70504', '3 nt', '1', ' ', ' ', ' ', ' ', '3+0=3', ' ', ...]
            """
            imgt_designation = parse_desc[0]
            gene = parse_desc[1]
            gene_segment = gene[3]
            latin = parse_desc[2].replace(" ", "_")
            functional = parse_desc[3]
            label = parse_desc[4]
            partial = parse_desc[-3].strip()

            # only want these labels, sometimes there are splice sites and intron sites which we don't need
            if label not in [
                "FR1-IMGT",
                "CDR1-IMGT",
                "FR2-IMGT",
                "CDR2-IMGT",
                "FR3-IMGT",
                "CDR3-IMGT",
            ]:
                continue

            # check we have not seen this gene,species,feature..
            if (imgt_designation, gene, latin, label) in seen:
                logger.warning("seen %s", (gene, latin, label))
                continue
            else:
                seen.append((imgt_designation, gene, latin, label))

            # Grab common name from latin name eg Homo_Sapiens = human
            common = IMGT_GB_LOOKUP[latin]
            gb_dataframe.append(
                {
                    "imgt_designation": imgt_designation,
                    "gene": gene,
                    "latin": latin,
                    "common": common,
                    "functional": functional,
                    "gene_segment": gene_segment,
                    "label": RENAME_DICT[label],
                    "nt_sequence_gaps": "",
                    "nt_sequence_no_gaps": str(sequence.seq).upper(),
                    "partial": partial,
                    "file": os.path.relpath(gb_file),
                    "fasta_header": sequence.description,
                    "source": "imgt_db",
                }
            )

    # set the index at gene, latin, common name and functional
    gb_dataframe = pd.DataFrame(gb_dataframe).set_index(["imgt_designation", "gene", "latin", "common", "functional"])

    # Pivot around label so FW1,FW2 etc will be columns insetead of rows
    # Multi index pivot since pandas freaks if you try and pivot with more than one column as index
    # gb_dataframe_pivot = multiindex_pivot(gb_dataframe, columns="label", values="nt_sequence").rename(
    #     RENAME_DICT, axis=1
    # )[PREFERRED_COLUMN_ORDER]

    # Tell everyone we parsed this by the IMGT_DB rather than by the Vgene fasta which we will do in the next part
    # gb_dataframe = gb_dataframe_pivot.reset_index()
    return gb_dataframe


def get_imgt_aa_gap(imgt_gapped_nt):
    raw_string = imgt_gapped_nt
    add_to_string = ""
    take_off_cdr3_overhange = raw_string[0 : len(raw_string) // 3 * 3]
    # add_to_nt = raw_string[len(raw_string)//3 * 3:]
    for codon_num in range(0, len(take_off_cdr3_overhange), 3):
        codon = raw_string[codon_num : codon_num + 3].upper()
        if codon == "...":
            add_to_string += "."
            continue

        try:
            add_to_string += codon_table[codon]
        except KeyError:
            raise CodonTable.TranslationError(f"{codon}, can't be translated")
    return add_to_string


def parse_imgt_v_quest_fasta(imgt_fasta):
    # Since IMGT v quest fastas is structured into latin named director
    imgt_v_quest_df = []
    for latin_species_path in glob.glob(imgt_fasta + "/*"):
        latin_species = latin_species_path.split("/")[-1]
        # grab the common name from the latin species
        common_name = REVERSE_IMGT_LOOKUP[latin_species]

        # There will be heavy,lambda,kappa and TCR V genes
        V_files_to_combine = glob.glob(latin_species_path + "/*/*V.fasta")
        D_files_to_combine = glob.glob(latin_species_path + "/*/*D.fasta")
        J_files_to_combine = glob.glob(latin_species_path + "/*/*J.fasta")
        # combine all those into one SEQIO object
        # combined_fasta = flatten([list(SeqIO.parse(file, "fasta")) for file in files_to_combine])

        # Sub dataframe will just act as lookup to see if we have defined this gene before
        # sub_dataframe = gb_dataframe[gb_dataframe["common"] == common_name].set_index("gene")
        """
        parse desc Description looks like this. We can grab everythign we need.
        ex.
        ['IMGT000050', 'IGKV(II)-8*01', 'Felis catus_Abyssinian', 'P', 'DONOR-SPLICE',
                                '70502..70504', '3 nt', '1', ' ', ' ', ' ', ' ', '3+0=3', ' ', ...]
         """

        # Keep a seen lookup since sometimes there are more than one entry because...well becuase IMGT is poorly constructued
        for file in V_files_to_combine + D_files_to_combine + J_files_to_combine:
            for sequence in SeqIO.parse(file, "fasta"):

                # this has the same parsed description above except without a feature label. We have to define that ourselvs
                parse_desc = sequence.description.split("|")
                imgt_designation = parse_desc[0]
                gene = parse_desc[1]
                gene_segment = gene[3]
                latin = parse_desc[2].replace(" ", "_")
                functional = parse_desc[3]
                partial = parse_desc[-3].strip()
                raw_string = str(sequence.seq)

                # IGHV or TRDV etc
                if gene_segment == "V":
                    for feature in IMGT_DEF_nt:
                        # Go through IMGT numbering
                        imgt_nt_start = IMGT_DEF_nt[feature]["start"]
                        imgt_nt_end = IMGT_DEF_nt[feature]["end"]
                        if not imgt_nt_end:
                            # from the raw string, grab just the feature part
                            imgt_feature_string = raw_string[imgt_nt_start:]
                        else:
                            imgt_feature_string = raw_string[imgt_nt_start - 1 : imgt_nt_end]

                        try:
                            imgt_feature_string_aa = get_imgt_aa_gap(imgt_feature_string)
                        except CodonTable.TranslationError as e:
                            logger.debug(f"{e}")
                            imgt_feature_string_aa = ""

                        imgt_v_quest_df.append(
                            {
                                "imgt_designation": imgt_designation,
                                "gene": gene,
                                "gene_segment": gene_segment,
                                "latin": latin,
                                "common": common_name,
                                "functional": functional,
                                "label": RENAME_DICT[feature],
                                "nt_sequence_no_gaps": imgt_feature_string.replace(".", "").upper(),
                                "nt_sequence_gaps": imgt_feature_string.upper(),
                                "aa_sequence_no_gaps": imgt_feature_string_aa.replace(".", "").upper(),
                                "aa_sequence_gaps": imgt_feature_string_aa.upper(),
                                "partial": partial,
                                "file": os.path.relpath(file),
                                "fasta_header": sequence.description,
                            }
                        )
                else:
                    feature = sequence.description.split("|")[4]
                    imgt_v_quest_df.append(
                        {
                            "imgt_designation": imgt_designation,
                            "gene": gene,
                            "gene_segment": gene_segment,
                            "latin": latin,
                            "common": common_name,
                            "functional": functional,
                            "label": feature,
                            "nt_sequence_no_gaps": raw_string.replace(".", "").upper(),
                            "nt_sequence_gaps": raw_string.upper(),
                            "aa_sequence_no_gaps": "",
                            "aa_sequence_gaps": "",
                            "partial": partial,
                            "file": os.path.relpath(file),
                            "fasta_header": sequence.description,
                        }
                    )

    # Turn this into a dataframe
    imgt_v_quest_df = pd.DataFrame(imgt_v_quest_df).set_index(
        ["imgt_designation", "gene", "latin", "common", "functional"]
    )

    # Pivot on labels to get cdrs,fws into columns
    # imgt_v_quest_df_pivot = multiindex_pivot(imgt_v_quest_df, columns="label", values="nt_sequence").rename(
    #     RENAME_DICT, axis=1
    # )[PREFERRED_COLUMN_ORDER]

    # Add a little note that tells us that these came from Vquest fasta
    # imgt_v_quest_df_pivot["source"] = "vquest_fasta"

    # Reset index
    # imgt_v_quest_df_pivot = imgt_v_quest_df_pivot.reset_index()
    imgt_v_quest_df.loc[:, "source"] = "vquest_fasta"
    return imgt_v_quest_df


def make_imgt_db(imgt_fasta, gb_path, engine):
    """Make IMGT DB file structure.

    Arguments:
       imgt_fasta {str} -- the path to Vquest fasta files
       gb_path {str} -- the path to IMGT_DB/
       engine {sqlalchmey.engine} -- the sqlite engine to write data to
    """

    gb_dataframe = parse_imgt_gb(gb_path)
    v_quest_dataframe = parse_imgt_v_quest_fasta(imgt_fasta)

    combined_df = pd.concat([v_quest_dataframe, gb_dataframe]).reset_index()
    total_entries = len(combined_df)
    logger.info(f"total information {total_entries}")
    per_source_index = ["imgt_designation", "gene", "common", "label", "source"]
    per_source_count = combined_df.groupby(per_source_index).size().sort_values()
    dups = per_source_count[per_source_count > 1].index
    logger.debug(
        f"{list(dups)} are duplicate. This can happen because IMGT puts them in two different IMGT fasta files. e.g, TRAV4-4/DV10*01 will be in TRAV and TRDV"
    )
    combined_df_no_dups = combined_df.drop_duplicates(combined_df.columns.drop("file"))

    # Make sure we removed the multiple fasta file entries coming from same source
    per_source_index = ["imgt_designation", "gene", "common", "label", "source"]
    per_source_count = combined_df_no_dups.groupby(per_source_index).size().sort_values()
    dups = per_source_count[per_source_count > 1].index
    if not dups.empty:
        raise Exception(f"{dups} have multiple entries in the same source. Check your fasta files")

    # Now we need to figure out the things with that have the same gene, label and imgt accession number but disagree with each other at the sequence level
    per_source_index = ["imgt_designation", "gene", "common", "label"]
    per_source_gb = combined_df_no_dups.groupby(per_source_index)
    per_source_count = per_source_gb.size().sort_values()
    multiple_entries_per_gene = per_source_count[per_source_count > 1].index
    logger.debug(f"{multiple_entries_per_gene} has been defined twice, once in vquest and once in imgt db")

    # place holder for things we will keep
    unique_indexes = []

    # keep track of things we discarded
    discarded_indexes = []

    # use groupby statement
    for group, group_df in per_source_gb:
        # if there is only 1 entry, then we keep the index
        if len(group_df) == 1:
            unique_indexes.append(group_df.index[0])
        # If there is two entries, then they come from vquest and imgt db
        elif len(group_df) == 2:
            # if the sequences agree with each other lets take vquest one
            if len(group_df["nt_sequence_no_gaps"].unique()) == 1:
                index_to_keep = group_df[group_df["source"] == "vquest_fasta"].index[0]
                index_to_discard = group_df[group_df["source"] == "imgt_db"].index[0]
                unique_indexes.append(index_to_keep)
                discarded_indexes.append(index_to_discard)
            else:  # if the sequences don't agree, lets yield to imgt_db
                index_to_keep = group_df[group_df["source"] == "imgt_db"].index[0]
                unique_indexes.append(index_to_keep)
                index_to_discard = group_df[group_df["source"] == "vquest_fasta"].index[0]
                discarded_indexes.append(index_to_discard)
        else:
            raise Exception(f"{group} has three entries listed, {group_df}, that should never happen")

    unique_combined_df = combined_df_no_dups.loc[unique_indexes]
    discarded_entries_df = combined_df_no_dups.loc[discarded_indexes]
    if (len(unique_combined_df) + len(discarded_entries_df)) != len(combined_df_no_dups):
        raise Exception(f"{len(unique_combined_df)} != {len(discarded_entries_df) + len(unique_combined_df)}")

    # our final source of ambiguity is if genes per species had two different accession numbers which is a...you guesed it...an imgt error:
    per_source_index = ["gene", "common", "label", "latin"]
    per_source_gb = unique_combined_df.groupby(per_source_index)
    per_source_count = per_source_gb.size().sort_values()
    genes_dont_match = per_source_count[per_source_count > 1].index

    logger.debug(f"droping {genes_dont_match} since they have differing accesion numbers in imgt")
    # add the duplicated entries
    discarded_entries_df = discarded_entries_df.append(
        unique_combined_df[unique_combined_df.duplicated(per_source_index)]
    )
    # drop the other ones
    unique_combined_df = unique_combined_df.drop_duplicates(per_source_index)
    if not (unique_combined_df.groupby(per_source_index).size().sort_values() == 1).all():
        raise Exception(
            f"{unique_combined_df.groupby(per_source_index).size()[unique_combined_df.groupby(per_source_index).size() != 1]} are duplicated"
        )
    # strip parenthensis and bracket
    unique_combined_df["functional"] = (
        unique_combined_df["functional"].str.strip("(").str.strip(")").str.strip("[").str.strip("]")
    )

    logger.info(f"Unique={len(unique_combined_df)} Total={len(combined_df)}")

    # No parse the othe rpath

    unique_combined_df.to_sql("imgt_unique", con=engine, if_exists="replace")
    combined_df.to_sql("imgt_all", con=engine, if_exists="replace")


def add_custom_sequences(custom_sequence_path, engine):
    custom_sequence_df = []
    for common_name_path in glob.glob(custom_sequence_path + "/*"):
        genes_files = glob.glob(common_name_path + "/*/*.fasta")
        common_name = os.path.basename(common_name_path)
        for fasta_file in genes_files:
            for sequence in SeqIO.parse(fasta_file, "fasta"):
                gene = sequence.description
                imgt_designation = "custom"
                gene_segment = gene[3]
                latin = ""
                functional = ""
                partial = ""
                raw_string = str(sequence.seq)

                # IGHV or TRDV etc
                if gene_segment == "V":
                    for feature in IMGT_DEF_nt:
                        # Go through IMGT numbering
                        imgt_nt_start = IMGT_DEF_nt[feature]["start"]
                        imgt_nt_end = IMGT_DEF_nt[feature]["end"]
                        if not imgt_nt_end:
                            # from the raw string, grab just the feature part
                            imgt_feature_string = raw_string[imgt_nt_start:]
                        else:
                            imgt_feature_string = raw_string[imgt_nt_start - 1 : imgt_nt_end]

                        custom_sequence_df.append(
                            {
                                "imgt_designation": imgt_designation,
                                "gene": gene,
                                "gene_segment": gene_segment,
                                "latin": latin,
                                "common": common_name,
                                "functional": functional,
                                "label": RENAME_DICT[feature],
                                "nt_sequence_no_gaps": imgt_feature_string.replace(".", "").upper(),
                                "nt_sequence_gaps": imgt_feature_string.upper(),
                                "partial": partial,
                                "file": os.path.relpath(fasta_file),
                                "fasta_header": sequence.description,
                            }
                        )
            else:
                if gene_segment.upper() == "J":
                    feature = "J-REGION"
                elif gene_segment == "D":
                    feature = "D-REGION"
                else:
                    raise Exception(f"Can't detct {sequence.description}")
                custom_sequence_df.append(
                    {
                        "imgt_designation": imgt_designation,
                        "gene": gene,
                        "gene_segment": gene_segment,
                        "latin": latin,
                        "common": common_name,
                        "functional": functional,
                        "label": feature,
                        "nt_sequence_no_gaps": raw_string.replace(".", "").upper(),
                        "nt_sequence_gaps": raw_string.upper(),
                        "partial": partial,
                        "file": os.path.relpath(fasta_file),
                        "fasta_header": sequence.description,
                    }
                )

    # Turn this into a dataframe
    custom_sequence_df = (
        pd.DataFrame(custom_sequence_df)
        .set_index(["imgt_designation", "gene", "latin", "common", "functional"])
        .reset_index()
    )
    custom_sequence_df.loc[:, "source"] = "custom"
    custom_sequence_df.to_sql("custom", con=engine, if_exists="replace")
