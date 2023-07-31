import tempfile

from sadie.reference import Reference, References

# create empty reference object
ref_class = Reference()
with tempfile.TemporaryDirectory() as tmpdirectory:
    # Add genes one at a time
    ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "source": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "source": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHJ6*01", "source": "imgt"})

    # call make_airr database on a path
    references = References()
    references.add_reference("human", ref_class)
    references.make_airr_database(tmpdirectory)
