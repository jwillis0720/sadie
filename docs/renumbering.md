#Renumbering Antibody Sequences

One of the hardest parts of working with antibody sequences is that have different definitions of numbering. There is also different definitions of where the framework and CDR regions are depending on the scheme. SADIE provides a simple interface to renumber antibody sequences to a common numbering scheme. We borrow heavily from the Antigen receptor Numbering and Receptor Classification ([ANARCI](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/anarci/))

---

## Single Sequence Annotation

```Python
{!docs_src/renumbering/simple_usage.py!}

```

|     | Id       | sequence                                                                                                                                          | domain_no | hmm_species | chain_type |     e-value |  score | seqstart_index | seqend_index | identity_species | v_gene       | v_identity | j_gene    | j_identity | Chain | Numbering                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | Insertion                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | Numbered_Sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | scheme  | region_definition | allowed_species | allowed_chains | fwr1_aa_gaps              | fwr1_aa_no_gaps           | cdr1_aa_gaps | cdr1_aa_no_gaps | fwr2_aa_gaps      | fwr2_aa_no_gaps   | cdr2_aa_gaps | cdr2_aa_no_gaps | fwr3_aa_gaps                           | fwr3_aa_no_gaps                        | cdr3_aa_gaps                            | cdr3_aa_no_gaps                         | fwr4_aa_gaps | fwr4_aa_no_gaps | leader | follow |
| --: | :------- | :------------------------------------------------------------------------------------------------------------------------------------------------ | --------: | :---------- | :--------- | ----------: | -----: | -------------: | -----------: | :--------------- | :----------- | ---------: | :-------- | ---------: | :---- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------ | :---------------- | :-------------- | :------------- | :------------------------ | :------------------------ | :----------- | :-------------- | :---------------- | :---------------- | :----------- | :-------------- | :------------------------------------- | :------------------------------------- | :-------------------------------------- | :-------------------------------------- | :----------- | :-------------- | :----- | :----- |
|   0 | VRC26.27 | QKQLVESGGGVVQPGRSLTLSCAASQFPFSHYGMHWVRQAPGKGLEWVASITNDGTKKYHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFCVRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDIWGHGTMVIVS |         0 | human       | H          | 1.65353e-43 | 134.25 |              0 |          144 | human            | IGHV3-30\*03 |        0.8 | IGHJ3\*02 |       0.64 | H     | [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 82, 82, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112] | ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'A', 'B', 'C', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', '', '', '', '', '', '', '', '', '', '', '', ''] | ['Q', 'K', 'Q', 'L', 'V', 'E', 'S', 'G', 'G', 'G', 'V', 'V', 'Q', 'P', 'G', 'R', 'S', 'L', 'T', 'L', 'S', 'C', 'A', 'A', 'S', 'Q', 'F', 'P', 'F', 'S', 'H', 'Y', 'G', 'M', 'H', 'W', 'V', 'R', 'Q', 'A', 'P', 'G', 'K', 'G', 'L', 'E', 'W', 'V', 'A', 'S', 'I', 'T', 'N', 'D', 'G', 'T', 'K', 'K', 'Y', 'H', 'G', 'E', 'S', 'V', 'W', 'D', 'R', 'F', 'R', 'I', 'S', 'R', 'D', 'N', 'S', 'K', 'N', 'T', 'L', 'F', 'L', 'Q', 'M', 'N', 'S', 'L', 'R', 'A', 'E', 'D', 'T', 'A', 'L', 'Y', 'F', 'C', 'V', 'R', 'D', 'Q', 'R', 'E', 'D', 'E', 'C', 'E', 'E', 'W', 'W', 'S', 'D', 'Y', 'Y', 'D', 'F', 'G', 'K', 'E', 'L', 'P', 'C', 'R', 'K', 'F', 'R', 'G', 'L', 'G', 'L', 'A', 'G', 'I', 'F', 'D', 'I', 'W', 'G', 'H', 'G', 'T', 'M', 'V', 'I', 'V', 'S'] | chothia | imgt              | human           | H,K,L          | QKQLVESGGGVVQPGRSLTLSCAAS | QKQLVESGGGVVQPGRSLTLSCAAS | QFPFSHYG     | QFPFSHYG        | MHWVRQAPGKGLEWVAS | MHWVRQAPGKGLEWVAS | ITNDGTKK     | ITNDGTKK        | YHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFC | YHGESVWDRFRISRDNSKNTLFLQMNSLRAEDTALYFC | VRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDI | VRDQREDECEEWWSDYYDFGKELPCRKFRGLGLAGIFDI | WGHGTMVIVS   | WGHGTMVIVS      |        |        |

The output will contain `<sadie.renumbering.result.NumberingResults'>` object. This object contains the following fields:

| Field             | Description                                      |
| :---------------- | :----------------------------------------------- |
| Id                | The sequence ID                                  |
| sequence          | sequence                                         |
| domain_no         | not used                                         |
| hmm_species       | the top species found in the HMM                 |
| chain_type        | the chain type, e.g 'H' or 'L'                   |
| e-value           | The e-value of the alignment                     |
| score             | The score for the alignment                      |
| seqstart_index    | where in the sequence does the alignment start   |
| seqend_index      | where in the sequence does the alignment end     |
| identity_species  | what species does the sequence aligns to best    |
| v_gene            | The top V gene                                   |
| v_identity        | V gene identity                                  |
| j_gene            | The top J gene in alignment                      |
| j_identity        | J gene identity                                  |
| Chain             | not used                                         |
| Numbering         | The numbering of the sequence stored as an array |
| Insertion         | The insertions if any stored as an array         |
| Numbered_Sequence | The matched sequence stored as an array          |
| scheme            | scheme, e.g. "kabat"                             |
| region_definition | CDR/FW definition                                |
| allowed_species   | allowed_species                                  |
| allowed_chains    | allowed_chains                                   |
| fwr1_aa_gaps      | fwr1_aa_gaps                                     |
| fwr1_aa_no_gaps   | fwr1_aa_no_gaps                                  |
| cdr1_aa_gaps      | cdr1_aa_gaps                                     |
| cdr1_aa_no_gaps   | cdr1_aa_no_gaps                                  |
| fwr2_aa_gaps      | fwr2_aa_gaps                                     |
| fwr2_aa_no_gaps   | fwr2_aa_no_gaps                                  |
| cdr2_aa_gaps      | cdr2_aa_gaps                                     |
| cdr2_aa_no_gaps   | cdr2_aa_no_gaps                                  |
| fwr3_aa_gaps      | fwr3_aa_gaps                                     |
| fwr3_aa_no_gaps   | fwr3_aa_no_gaps                                  |
| cdr3_aa_gaps      | cdr3_aa_gaps                                     |
| cdr3_aa_no_gaps   | cdr3_aa_no_gaps                                  |
| fwr4_aa_gaps      | fwr4_aa_gaps                                     |
| fwr4_aa_no_gaps   | fwr4_aa_no_gaps                                  |
| leader            | what sequences come before the alignment         |
| follow            | what sequences come after the alignment          |

The `NumberingResults` is a pandas dataframe instance so it can be used like one. It also contains an alignment table that looks like the following.

```Python
{!docs_src/renumbering/alignment_table.py!}
```

!!! warning multiprocessing

    Multiprocessing must be wrapped in a function at the current time if you set run_multi=True. It will also work inside a Jupyter notebook cell.

which stores a handy alignment table of the sequence.

|     | Id       | chain_type | scheme  | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  | 13  | 14  | 15  | 16  | 17  | 18  | 19  | 20  | 21  | 22  | 23  | 24  | 25  | 26  | 27  | 28  | 29  | 30  | 31  | 32  | 33  | 34  | 35  | 36  | 37  | 38  | 39  | 40  | 41  | 42  | 43  | 44  | 45  | 46  | 47  | 48  | 49  | 50  | 51  | 52  | 52A | 53  | 54  | 55  | 56  | 57  | 58  | 59  | 60  | 61  | 62  | 63  | 64  | 65  | 66  | 67  | 68  | 69  | 70  | 71  | 72  | 73  | 74  | 75  | 76  | 77  | 78  | 79  | 80  | 81  | 82  | 82A | 82B | 82C | 83  | 84  | 85  | 86  | 87  | 88  | 89  | 90  | 91  | 92  | 93  | 94  | 95  | 96  | 97  | 98  | 99  | 100 | 100A | 100B | 100C | 100D | 100E | 100F | 100G | 100H | 100I | 100J | 100K | 100L | 100M | 100N | 100O | 100P | 100Q | 100R | 100S | 100T | 100U | 100V | 100W | 100X | 100Y | 100Z | 100a | 100b | 100c | 101 | 102 | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 110 | 111 | 112 |
| --: | :------- | :--------- | :------ | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
|   0 | VRC26.27 | H          | chothia | Q   | K   | Q   | L   | V   | E   | S   | G   | G   | G   | V   | V   | Q   | P   | G   | R   | S   | L   | T   | L   | S   | C   | A   | A   | S   | Q   | F   | P   | F   | S   | H   | Y   | G   | M   | H   | W   | V   | R   | Q   | A   | P   | G   | K   | G   | L   | E   | W   | V   | A   | S   | I   | T   | N   | D   | G   | T   | K   | K   | Y   | H   | G   | E   | S   | V   | W   | D   | R   | F   | R   | I   | S   | R   | D   | N   | S   | K   | N   | T   | L   | F   | L   | Q   | M   | N   | S   | L   | R   | A   | E   | D   | T   | A   | L   | Y   | F   | C   | V   | R   | D   | Q   | R   | E   | D   | E   | C    | E    | E    | W    | W    | S    | D    | Y    | Y    | D    | F    | G    | K    | E    | L    | P    | C    | R    | K    | F    | R    | G    | L    | G    | L    | A    | G    | I    | F    | D   | I   | W   | G   | H   | G   | T   | M   | V   | I   | V   | S   |

## Multiple Sequence Numbering


You can also renumber a fasta file.

=== ":material-console-line: Command Line Usage"

    <div class="termy">

    ```console
    $ sadie renumbering -q catnap_aa_heavy_sample.fasta
    ```
    </div>

    The output will be `catnap_aa_heavy_sample_numbering_segment.csv` which will be the table from the Numbering Results and `catnap_aa_heavy_sample_numbering_alignment.csv` which will be the alignment table.

=== " :material-api: Python"

    ```Python
    {!>docs_src/renumbering/multi_use.py!}
    ```

## Schemes

These are the current numbering schemes we have implemented.

| Scheme  | Description                                                                                              |
| :------ | :------------------------------------------------------------------------------------------------------- |
| chothia | [Chothia numbering scheme](https://www.chemogenomix.com/chothia-antibody-numbering)                      |
| kabat   | [Kabat numbering scheme](https://www.chemogenomix.com/kabat-antibody-numbering)                          |
| imgt    | [IMGT numbering scheme](https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html) |

## Region definitions


Given a numbering scheme we can define CDRS and frameworks with the following definitions.

['imgt', 'kabat', 'chothia', 'abm', 'contact', 'scdr']

| Region Definition | Description |
| :---------------- | :---------- |
| imgt              | IMGT        |
| kabat             | Kabat       |
| chothia           | Chothia     |
| abm               | ABM         |
| contact           | Contact     |
| scdr              | SCDR        |

The following is a description of each definition taken from [the Martin group](http://www.bioinf.org.uk/abs/info.html)

- The Kabat definition is based on sequence variability and is the most commonly used
- The Chothia definition is based on the location of the structural loop regions - see more detail at the bottom of this section
- The AbM definition is a compromise between the two used by Oxford Molecular's AbM antibody modelling software
- The contact definition has been recently introduced by us and is based on an analysis of the available complex crystal structures. This definition is likely to be the most useful for people wishing to perform mutagenesis to modify the affinity of an antibody since these are residues which take part in interactions with antigen. Lists of CDR residues making contact in each antibody with summary data for each CDR
- SCDR is the longest CDR definition for each region. It's used in industry

For a great review of the numbering schemes and region definitions, see this [paper](https://www.frontiersin.org/articles/10.3389/fimmu.2018.02278/full)
