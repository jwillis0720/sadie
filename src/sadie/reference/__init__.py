__version__ = "0.3.6"
import requests
import logging

logger = logging.getLogger("reference")
# handle database IO
# database_root_path = Path(__file__).parent.joinpath("data")
# custom_database_path = database_root_path.joinpath("custom-g3.json.gz")
# imgt_database_path = database_root_path.joinpath("imgt-g3.json.gz")


def get_database(source: str) -> list:
    if source not in ["custom", "imgt"]:
        raise ValueError("Invalid database source, needs to be 'custom' or 'imt'")
    response = requests.get(f"https://g3.jordanrwillis.com/api/v1/genes?source={source}&limit=-1")
    logger.info(f"{source} database response: {response.status_code}")
    return response.json()


def get_loaded_database() -> dict:
    """Get G3 database , wrapped in this function so we don't call on response when module is loaded

    Returns
    -------
    dict
        Dictionary of G3 database. custom or imgt
    """
    return {"custom": get_database("custom"), "imgt": get_database("imgt")}


__all__ = ["get_database", "get_loaded_database"]
