from sadie.reference import Reference
import tempfile

# create empty reference object
ref_class = Reference()
with tempfile.TemporaryDirectory() as tmpdirectory:
    # Add genes one at a time
    ref_class.add_gene({"species": "human", "gene": "IGHV1-69*01", "database": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHD3-3*01", "database": "imgt"})
    ref_class.add_gene({"species": "human", "gene": "IGHJ6*01", "database": "imgt"})

    # call make_airr database on a path
    ref_class.make_airr_database(tmpdirectory)
    # use tempdirectory using AIRR module
    # ...
