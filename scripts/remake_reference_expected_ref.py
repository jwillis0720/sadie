# Use this script if you ever have changes with major structure of the reference database that are expected in the integration tests

# import os
import json
from glob import glob
from pathlib import Path

# import pprint


# BLAST DIR
PATH_TO_DATA = "./src/sadie/airr/data/germlines/Ig/blastdb/**/*"
OUTPUT_JSON = "./tests/data/fixtures/reference/blast_dir.json"
blast_db_struct = [
    "/".join(list(Path(i).parts[5:])) for i in glob(PATH_TO_DATA, recursive=True) if not Path(i).is_dir()
]
json.dump(blast_db_struct, open(OUTPUT_JSON, "w"), indent=2)

# Aux DIR
PATH_TO_DATA = "./src/sadie/airr/data/germlines/aux_db/**/*.aux"
OUTPUT_JSON = "./tests/data/fixtures/reference/aux.json"
aux_db_struct = ["/".join(list(Path(i).parts[5:])) for i in glob(PATH_TO_DATA, recursive=True) if not Path(i).is_dir()]
json.dump(aux_db_struct, open(OUTPUT_JSON, "w"), indent=2)

# internal DIR
PATH_TO_DATA = "./src/sadie/airr/data/germlines/Ig/internal_data/**/*"
OUTPUT_JSON = "./tests/data/fixtures/reference/internal.json"
aux_db_struct = ["/".join(list(Path(i).parts[5:])) for i in glob(PATH_TO_DATA, recursive=True) if not Path(i).is_dir()]
json.dump(aux_db_struct, open(OUTPUT_JSON, "w"), indent=2)
