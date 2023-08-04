from sadie.receptor.rearrangment import ReceptorChain
from Bio import SeqIO

# get a sequence
vrc01_heavy_sequecne = SeqIO.read("vrc01_heavy.fasta", "fasta")

# make a receptor chain model from a single sequence...much easier!
receptor_chain = ReceptorChain.from_single("vrc01_heavy", vrc01_heavy_sequecne.seq)
