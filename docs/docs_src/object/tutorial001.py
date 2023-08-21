from pprint import pprint

from Bio import SeqIO

from sadie.receptor.rearrangment import InputSequence

vrc01_heavy_sequecne = SeqIO.read("vrc01_heavy.fasta", "fasta")

# make an input sequence model
input_sequence_model = InputSequence(
    sequence_id=vrc01_heavy_sequecne.name,
    sequence=vrc01_heavy_sequecne.seq,
    raw_sequence=vrc01_heavy_sequecne.seq,
)

# Print out dictionary to see
pprint(input_sequence_model.__dict__)
