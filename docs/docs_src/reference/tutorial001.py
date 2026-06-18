import gzip
import json
from importlib.resources import files

# The live G3 API is retired; gene records ship inside the package and are read offline.
collection = files("sadie.reference") / "data" / "imgt-g3.json.gz"
with gzip.open(collection, "rt") as handle:
    genes = json.load(handle)

human_v = [g for g in genes if g["common"] == "human" and g["gene_segment"] == "V"][:5]
print(json.dumps(human_v, indent=4))
json.dump(human_v, open("human_v.json", "w"), indent=4)
