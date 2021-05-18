__version__ = "0.3.1"
from pathlib import Path
import json
import gzip


# handle database IO
database_root_path = Path(__file__).parent.joinpath("data")
custom_database_path = database_root_path.joinpath("custom-g3.json.gz")
imgt_database_path = database_root_path.joinpath("imgt-g3.json.gz")
loaded_database = {
    "custom": json.load(gzip.open(custom_database_path, "rt")),
    "imgt": json.load(gzip.open(imgt_database_path, "rt")),
}


def get_loaded_database() -> dict:
    return loaded_database


__all__ = ["database_root_path", "custom_database_path", "imgt_database_path", loaded_database, get_loaded_database]
