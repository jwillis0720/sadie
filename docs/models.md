# Antibody Objects

This module provides convenience classes to store the annotations of an antibody or TCR. First, let's start at a very low level. These are pythonic objects that model the data we expect in an AIRR compliant data format. They are divided by [AIRR 1.3 Rearrangement category](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

- Input Sequence
- Primary Annotations
- Alignment Annotations
- Alignment Positions
- RegionSequences
- RegionPositions
- Junction Lengths

All of these are combined as a `Receptor Chain` Object.

# Required Models

These are the required models when constructing a `Receptor Chain` object.

## Input Sequence

The input sequence object contains the input sequence to the V(D)J assignment process.

```Python
{!docs_src/object/tutorial001.py!}
```

## Primary Annotations

```Python
{!docs_src/object/tutorial002.py!}
```

## Alignment Annotations

```Python
{!docs_src/object/tutorial003.py!}
```

# Optional Models

These are optional but recommended models when constructing a `Receptor Chain` object.

## Alignment Positions

```Python
{!docs_src/object/tutorial004.py!}
```

## Region Sequences

```Python
{!docs_src/object/tutorial005.py!}
```

## Region Positions

```Python
{!docs_src/object/tutorial006.py!}
```

## Junction Lengths

```Python
{!docs_src/object/tutorial007.py!}
```

# Receptor Chain

The `ReceptorChain` object is a container for all of the above models. It is the primary object that is used to store the annotations of an antibody or TCR. This is what it looks like all combined.

```Python
{!docs_src/object/tutorial008.py!}
```

outputs:

```
input_sequence
--------------
sequence_id : VRC01_heavy
sequence : CAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCTGGATATGAATTTATTGATTGTACGCTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGGCTGAAGCCTCGGGGGGGGGCCGTCAACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGGGGCCGGGGCACCCCGGTCATCGTCTCATCAGGGGGG
raw_sequence : CAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCTGGATATGAATTTATTGATTGTACGCTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGGCTGAAGCCTCGGGGGGGGGCCGTCAACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGGGGCCGGGGCACCCCGGTCATCGTCTCATCAGGGGGG
sequence_aa : None
category : category='input'


primary_annotations
-------------------
rev_comp : False
productive : True
vj_in_frame : True
stop_codon : False
complete_vdj : True
locus : IGH
v_call : IGHV1-2*02
d_call : ['IGHD3-16*01', 'IGHD3-16*02']
d2_call : None
j_call : IGHJ1*01
c_call : IGHG1*01
v_call_top : IGHV1-2*02
v_call_top_gene : None
v_call_top_allele : None
d_call_top : IGHD3-16*01
d_call_gene : None
d_call_allele : None
j_call_top : IGHJ1*01
j_call_top_gene : None
j_call_top_allele : None
c_call_allele : None
reference_name : None
category : category='primary_annotations'


alignment_annotations
---------------------
sequence_alignment : CAGGTGCAGCTGGTGCAGTCTGGGGGTCAGATGAAGAAGCCTGGCGAGTCGATGAGAATTTCTTGTCGGGCTTCTGGATATGAATTTATTGATTGTACGCTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGGCTGAAGCCTCGGGGGGGGGCCGTCAACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGGGGCCGGGGCACCCCGGTCATCGTCTCATCAG
sequence_alignment_aa : QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS
germline_alignment : CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGNNNNNNNNNNNNTGATTACGTTTGGGACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG
germline_alignment_aa : QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAXXXXXDYVWDFQHWGQGTLVTVSS
v_score : 168.2
v_identity : 0.6825
v_support : 6.796e-44
v_cigar : 6S293M76S3N
d_score : 17.8
d_identity : 0.85
d_support : 0.5755
d_cigar : 311S6N14M50S17N
d2_score : None
d2_identity : None
d2_support : None
d2_cigar : None
j_score : 52.6
j_identity : 0.86
j_support : 5.727e-11
j_cigar : 325S7N45M5S
junction : TGTACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACACTGG
junction_aa : CTRGKNCDYNWDFEHW
np1 : GGGAAAAAACTG
np1_aa : None
np2 : None
np2_aa : None
np3 : None
np3_aa : None
c_score : 100.0
c_identity : 1.0
c_support : 1e-44
c_cigar : 6S293M76S3N
category : category='alignment_annotations'


alignment_positions
-------------------
v_sequence_start : 7
v_sequence_end : 299
v_germline_start : 1
v_germline_end : 293
v_alignment_start : 1
v_alignment_end : 293
d_sequence_start : 312
d_sequence_end : 325
d_germline_start : 7
d_germline_end : 20
d_alignment_start : 306
d_alignment_end : 319
d2_sequence_start : None
d2_sequence_end : None
d2_germline_start : None
d2_germline_end : None
d2_alignment_start : None
d2_alignment_end : None
j_sequence_start : 326
j_sequence_end : 370
j_germline_start : 8
j_germline_end : 52
j_alignment_start : 320
j_alignment_end : 364
category : category='alignment_positions'


region_sequences
----------------
fwr1 : None
fwr1_aa : QVQLVQSGGQMKKPGESMRISCRAS
cdr1 : GGATATGAATTTATTGATTGTACG
cdr1_aa : GYEFIDCT
fwr2 : CTAAATTGGATTCGTCTGGCCCCCGGAAAAAGGCCTGAGTGGATGGGATGG
fwr2_aa : LNWIRLAPGKRPEWMGW
cdr2 : CTGAAGCCTCGGGGGGGGGCCGTC
cdr2_aa : LKPRGGAV
fwr3 : AACTACGCACGTCCACTTCAGGGCAGAGTGACCATGACTCGAGACGTTTATTCCGACACAGCCTTTTTGGAGCTGCGCTCGTTGACAGTAGACGACACGGCCGTCTACTTTTGT
fwr3_aa : NYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFC
cdr3 : ACTAGGGGAAAAAACTGTGATTACAATTGGGACTTCGAACAC
cdr3_aa : TRGKNCDYNWDFEH
fwr4 : TGGGGCCGGGGCACCCCGGTCATCGTCTCATCA
fwr4_aa : WGRGTPVIVSS
category : category='region_sequence_annotations'


region_positions
----------------
fwr1_start : 7
fwr1_end : 81
cdr1_start : 82
cdr1_end : 105
fwr2_start : 106
fwr2_end : 156
cdr2_start : 157
cdr2_end : 180
fwr3_start : 181
fwr3_end : 294
cdr3_start : 295
cdr3_end : 336
fwr4_start : 337
fwr4_end : 369
category : category='region_positions'


junction_lengths
----------------
junction_length : 48
junction_aa_length : None
np1_length : None
np2_length : None
np3_length : None
n1_length : None
n2_length : None
n3_length : None
p3v_length : None
p5d_length : None
p3d_length : None
p5d2_length : None
p3d2_length : None
p5j_length : None
category : category='junction_lengths'
```

# Convenience Functions

But what if you don't want to fill in all that information? Where do you even get all that meta-info.

```Python
{!docs_src/object/tutorial009.py!}
```

!!! info

    You can use all these objects for pythonic syntatic sugar but you will probably just want to use the `ReceptorChain` object for most of your use cases.
