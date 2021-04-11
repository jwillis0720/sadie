__version__ = "0.2.14"
import os
import json
import gzip

database_path = os.path.join(os.path.dirname(__file__), "data/ig_database.json.gz")
loaded_database = json.load(gzip.open(database_path, "rt"))
__all__ = ["database_path", "loaded_database"]
