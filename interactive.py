# %%
from yaml import load
from pprint import pprint
from yaml import Loader
from sadie.reference import loaded_database

# %%
pprint(load(open("reference.yml"), Loader=Loader))


# %%
with open("custom.yml", "w") as f:
    for species in ["human", "macaque", "cat", "dog", "mouse", "rat", "alpaca"]:
        f.write("\t\t" + species + ":\n")
        f.write(f"\t\t\t{species}:\n")
        for x in list(
            map(
                lambda x: "- " + x["gene"],
                sorted(
                    list(
                        filter(
                            lambda x: x["common"] == species and x["receptor"] == "Ig",
                            loaded_database,
                        )
                    ),
                    key=lambda x: x["gene"],
                ),
            )
        ):
            f.write("\t\t\t" + x + "\n")
# %%
set(map(lambda x: x["common"], loaded_database))
# %%
