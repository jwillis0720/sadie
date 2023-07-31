import tempfile

from sadie.reference import Reference, References
from sadie.reference.yaml import YamlRef

# enter no file to use reference.yml
yml_ref = YamlRef()

# create empty reference object
ref_class = Reference()

# references class
references = References()

# Iterate through YamlRef
for name in yml_ref:
    # these are dictionary entries
    for database in yml_ref[name]:
        species = yml_ref[name][database]
        for single in species:
            single_species = single
            # gene is a list of genes
            genes = yml_ref[name][database][single_species]
            if database == "custom" and single_species == "macaque":  # only want cat and custom
                only_vs = list(filter(lambda x: x[3] == "V", genes))  # only get v genes, lookup third letter for this
                for gene in only_vs[:5]:  # only getting first 5
                    ref_class.add_gene({"gene": gene, "species": single_species, "source": database})

references.add_reference("small_macaque", ref_class)
