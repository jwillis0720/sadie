from sadie.reference import Reference
from sadie.reference.yaml import YamlRef
import tempfile

# enter no file to use reference.yml
yml_ref = YamlRef()

# create empty reference object
ref_class = Reference()

# Iterate through YamlRef
for entry in yml_ref:
    # these are dictionary entries
    db = entry["database"]
    species = entry["species"]
    sub_species = entry["sub_species"]

    # gene is a list of genes
    genes = entry["gene"]

    # only keep custom
    if db == "custom" and species == "cat":  # only want cat and custom
        only_vs = list(filter(lambda x: x[3] == "V", genes))  # only get v genes, lookup third letter for this
        for gene in only_vs[:5]:  # only getting first 5
            ref_class.add_gene({"gene": gene, "species": species, "sub_species": sub_species, "database": db})

# now add J6 to make a cat/dog. Can't mix and match database types!
ref_class.add_gene({"gene": "IGHJ6*01", "species": "dog", "sub_species": "dog", "database": "custom"})

with tempfile.TemporaryDirectory() as named_temp:
    ref_class.make_airr_database(named_temp)
