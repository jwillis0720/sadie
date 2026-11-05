import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional
from collections.abc import Mapping

# Module-level logger
logger = logging.getLogger(__name__)

RECEPTORS = ["IG", "TR"]
IMGT_DEF_nt = {
    "FW1": {"start": 1, "end": 78},
    "CDR1": {"start": 79, "end": 114},
    "FW2": {"start": 115, "end": 165},
    "CDR2": {"start": 166, "end": 195},
    "FW3": {"start": 196, "end": 312},
    "CDR3": {"start": 312, "end": ""},
    "V-REGION": {"start": 0, "end": ""},
}

IMGT_GB_LOOKUP = {
    "Canis_lupus_familiaris_boxer": "dog",
    "Felis_catus_Abyssinian": "cat",
    "Canis_lupus_familiaris_Canis_lupus_familiaris_boxer": "dog",
    "Rattus_norvegicus_BN;_Sprague-Dawley": "rat",
    "Rattus_norvegicus_BN/SsNHsdMCW": "rat",
    "Rattus_norvegicus": "rat",
    "Mus_musculus_C57BL/6": "mouse",
    "Mus_musculus_BALB/c": "mouse",
    "Mus_musculus_C57BL/6J": "mouse",
    "Mus_musculus_MRL/lpr": "mouse",
    "Mus_musculus": "mouse",
    "Mus_musculus_A/J": "mouse",
    "Mus_musculus_C57BL/10": "mouse",
    "Mus_musculus_129/Sv": "mouse",
    "Mus_musculus_NZB": "mouse",
    "Mus_musculus_I/St": "mouse",
    "Mus_musculus_NFS": "mouse",
    "Mus_musculus_BALB.K": "mouse",
    "Mus_musculus_C3H": "mouse",
    "Mus_musculus_NZB/BINJ": "mouse",
    "Mus_musculus_CE/J": "mouse",
    "Mus_musculus_PERU": "mouse",
    "Mus_musculus_AKR": "mouse",
    "Mus_musculus_domesticus": "mouse",
    "Mus_musculus_O20/A": "mouse",
    "Mus_musculus_castaneus": "mouse",
    "Mus_musculus_molossinus_MOLF/Ei": "mouse",
    "Mus_musculus_musculus": "mouse",
    "Mus_musculus_castaneus_CAST/Ei": "mouse",
    "Mus_musculus_C58": "mouse",
    "Mus_musculus_SK": "mouse",
    "Mus_musculus_PERA": "mouse",
    "Mus_musculus_MRL": "mouse",
    "Mus_musculus_MBK": "mouse",
    "Mus_musculus_PWK": "mouse",
    "Mus_musculus_MAI": "mouse",
    "Homo_sapiens": "human",
    "Vicugna_pacos": "alpaca",
}


IMGT_LOOKUP = {
    "human": "Homo_sapiens",
    "cow": "Bos_taurus",
    "camel": "Camelus_dromedarius",
    "dog": "Canis_lupus_familiaris",
    "cat": "Felis_catus",
    "junglefowl": "Gallus_gallus",
    "night_monkey": "Aotus_nancymaae",
    "goat": "Capra_hircus",
    "sharks": "Chondrichthyes",
    "zebrafish": "Danio_rerio",
    "horse": "Equus_caballus",
    "cod": "Gadus_morhua",
    "catfish": "Ictalurus_punctatus",
    "crabmacaque": "Macaca_fascicularis",
    "macaque": "Macaca_mulatta",
    "mouse": "Mus_musculus",
    "ferret": "Mustela_putorius_furo",
    "nhp": "Nonhuman_primates",
    "trout": "Oncorhynchus_mykiss",
    "platypus": "Ornithorhynchus_anatinus",
    "rabbit": "Oryctolagus_cuniculus",
    "sheep": "Ovis_aries",
    "rat": "Rattus_norvegicus",
    "salmon": "Salmo_salar",
    "boar": "Sus_scrofa",
    "teleosts": "Teleostei",
    "dolphin": "Tursiops_truncatus",
    "alpaca": "Vicugna_pacos",
}

REVERSE_IMGT_LOOKUP = {v: k for k, v in IMGT_LOOKUP.items()}

BLAST_CONVENTION = {"IG": "Ig", "TR": "TCR"}


SEGMENTS_INTERNAL_DATA = {
    "TR": ["TRAV", "TRBV", "TRDV", "TRGV"],
    "IG": ["IGHV", "IGKV", "IGLV"],
}
SEGMENTS = {
    "TR": ["TRAV", "TRBD", "TRBJ", "TRBV", "TRDD", "TRDJ", "TRDV", "TRGV", "TRGJ"],
    "IG": ["IGHD", "IGHJ", "IGHV", "IGKJ", "IGKV", "IGLJ", "IGLV"],
}


# JSON motif registry loader and cache
_MOTIF_REGISTRY_CACHE: Optional[Dict[str, Any]] = None
_MOTIF_LOOKUP_CACHE: Optional[Dict[str, Any]] = None
_MOTIF_PROVENANCE_CACHE: Optional[Dict[str, Any]] = None


def _get_motif_json_path() -> Path:
    """Get the path to the motif JSON file."""
    return Path(__file__).parent / "data" / "j_gene_motif.json"


def _load_motif_registry() -> Dict[str, Any]:
    """
    Load the complete motif registry from JSON file.

    Returns:
        Complete motif registry with provenance metadata

    Raises:
        FileNotFoundError: If JSON file is not found
        json.JSONDecodeError: If JSON file is malformed
    """
    global _MOTIF_REGISTRY_CACHE

    if _MOTIF_REGISTRY_CACHE is not None:
        return _MOTIF_REGISTRY_CACHE

    json_path = _get_motif_json_path()

    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            registry = json.load(f)

        # Validate loaded data structure
        if not isinstance(registry, dict):
            raise ValueError(f"Motif registry must be a dictionary, got {type(registry)}")

        if len(registry) == 0:
            raise ValueError("Motif registry is empty")

        # Basic structure validation
        for species, species_data in registry.items():
            if not isinstance(species_data, dict):
                raise ValueError(f"Species data for '{species}' must be a dictionary, got {type(species_data)}")

            if '_provenance' not in species_data:
                raise ValueError(f"Species '{species}' missing required _provenance metadata")

            # Check for at least one pattern
            pattern_count = sum(1 for k, v in species_data.items()
                              if k not in ('_provenance', 'ignore') and isinstance(v, str))
            if pattern_count == 0:
                logger.warning(f"Species '{species}' has no motif patterns")

        _MOTIF_REGISTRY_CACHE = registry
        logger.info(f"Loaded motif registry with {len(registry)} species from {json_path}")

        # Log validation counts
        try:
            imgt_validated = sum(1 for species_data in registry.values()
                                if species_data.get('_provenance', {}).get('imgt_validated', False))
            legacy_count = len(registry) - imgt_validated
            logger.info(f"Motif registry: {imgt_validated} IMGT-validated, {legacy_count} legacy species")
        except Exception as e:
            logger.warning(f"Could not compute validation counts: {e}")

        return registry

    except FileNotFoundError:
        logger.error(f"Motif registry file not found: {json_path}")
        raise FileNotFoundError(
            f"Motif registry file not found at {json_path}. "
            f"Ensure the file exists or check the module installation."
        )
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in motif registry: {e}")
        raise json.JSONDecodeError(
            f"Malformed JSON in motif registry {json_path}: {e}",
            e.doc, e.pos
        )
    except (ValueError, TypeError) as e:
        logger.error(f"Invalid data structure in motif registry: {e}")
        raise ValueError(
            f"Motif registry data is malformed at {json_path}: {e}. "
            f"Check that the file contains valid motif registry data."
        )
    except PermissionError as e:
        logger.error(f"Permission denied reading motif registry: {json_path}")
        raise PermissionError(
            f"Cannot read motif registry file {json_path}: {e}. "
            f"Check file permissions."
        )
    except Exception as e:
        logger.error(f"Unexpected error loading motif registry: {e}")
        raise RuntimeError(
            f"Failed to load motif registry from {json_path}: {e}. "
            f"This may indicate a corrupted file or installation issue."
        )


def _filter_provenance_metadata(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Filter out provenance metadata to maintain backward compatibility.

    Args:
        data: Complete motif registry data

    Returns:
        Filtered data without _provenance fields
    """
    filtered = {}

    for species, species_data in data.items():
        # Create copy without _provenance
        filtered_species = {
            key: value for key, value in species_data.items()
            if key != "_provenance"
        }
        filtered[species] = filtered_species

    return filtered


def _extract_provenance_metadata(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Extract provenance metadata from the complete registry.

    Args:
        data: Complete motif registry data

    Returns:
        Dictionary mapping species to provenance metadata
    """
    provenance = {}

    for species, species_data in data.items():
        if "_provenance" in species_data:
            provenance[species] = species_data["_provenance"]

    return provenance


def get_motif_lookup() -> Dict[str, Dict[str, Any]]:
    """
    Get the backward-compatible MOTIF_LOOKUP dictionary.

    Returns:
        Motif lookup dictionary without provenance metadata
    """
    global _MOTIF_LOOKUP_CACHE

    if _MOTIF_LOOKUP_CACHE is not None:
        return _MOTIF_LOOKUP_CACHE

    registry = _load_motif_registry()
    lookup = _filter_provenance_metadata(registry)
    _MOTIF_LOOKUP_CACHE = lookup

    logger.debug(f"Generated MOTIF_LOOKUP for {len(lookup)} species (provenance filtered)")
    return lookup


def get_motif_provenance() -> Dict[str, Dict[str, Any]]:
    """
    Get the provenance metadata for all species.

    Returns:
        Dictionary mapping species to provenance metadata
    """
    global _MOTIF_PROVENANCE_CACHE

    if _MOTIF_PROVENANCE_CACHE is not None:
        return _MOTIF_PROVENANCE_CACHE

    registry = _load_motif_registry()
    provenance = _extract_provenance_metadata(registry)
    _MOTIF_PROVENANCE_CACHE = provenance

    logger.debug(f"Generated MOTIF_PROVENANCE for {len(provenance)} species")
    return provenance


# Lazy-loaded backward-compatible MOTIF_LOOKUP
# This preserves the original API while loading from JSON
class _MotifLookupProxy(Mapping):
    """Proxy object that loads MOTIF_LOOKUP on first access."""

    def __getitem__(self, key):
        try:
            return get_motif_lookup()[key]
        except KeyError:
            available_species = list(get_motif_lookup().keys())[:10]  # Show first 10
            raise KeyError(
                f"Species '{key}' not found in MOTIF_LOOKUP. "
                f"Available species: {available_species}... "
                f"(total: {len(get_motif_lookup())}). "
                f"Check species name spelling (lowercase with underscores)."
            )
        except Exception as e:
            logger.error(f"Error accessing MOTIF_LOOKUP['{key}']: {e}")
            raise RuntimeError(
                f"Failed to access motif data for species '{key}': {e}. "
                f"This may indicate a data loading issue."
            )

    def __iter__(self):
        try:
            return iter(get_motif_lookup())
        except Exception as e:
            logger.error(f"Error iterating MOTIF_LOOKUP: {e}")
            raise RuntimeError(f"Failed to iterate motif data: {e}")

    def __len__(self):
        try:
            return len(get_motif_lookup())
        except Exception as e:
            logger.error(f"Error getting MOTIF_LOOKUP length: {e}")
            raise RuntimeError(f"Failed to get motif data length: {e}")

    def keys(self):
        try:
            return get_motif_lookup().keys()
        except Exception as e:
            logger.error(f"Error getting MOTIF_LOOKUP keys: {e}")
            raise RuntimeError(f"Failed to get motif data keys: {e}")

    def values(self):
        try:
            return get_motif_lookup().values()
        except Exception as e:
            logger.error(f"Error getting MOTIF_LOOKUP values: {e}")
            raise RuntimeError(f"Failed to get motif data values: {e}")

    def items(self):
        try:
            return get_motif_lookup().items()
        except Exception as e:
            logger.error(f"Error getting MOTIF_LOOKUP items: {e}")
            raise RuntimeError(f"Failed to get motif data items: {e}")

    def get(self, key, default=None):
        try:
            return get_motif_lookup().get(key, default)
        except Exception as e:
            logger.error(f"Error getting MOTIF_LOOKUP key '{key}': {e}")
            # Return default rather than raise, to maintain get() semantics
            logger.warning(f"Returning default value due to error: {e}")
            return default

    def __contains__(self, key):
        try:
            return key in get_motif_lookup()
        except Exception as e:
            logger.error(f"Error checking MOTIF_LOOKUP contains '{key}': {e}")
            # Return False rather than raise, to maintain contains semantics
            return False


# Backward-compatible MOTIF_LOOKUP that loads from JSON
MOTIF_LOOKUP = _MotifLookupProxy()


J_SEGMENTS = {"IG": ["IGHJ", "IGKJ", "IGKJ"], "TR": ["TRAJ", "TRBJ", "TRGJ", "TRDJ"]}
RENAME_DICT = {
    "CDR1-IMGT": "cdr1_nt",
    "FR1-IMGT": "fwr1_nt",
    "FR2-IMGT": "fwr2_nt",
    "FR3-IMGT": "fwr3_nt",
    "CDR2-IMGT": "cdr2_nt",
    "CDR3-IMGT": "cdr3_nt",
    "FW1": "fwr1_nt",
    "FW2": "fwr2_nt",
    "FW3": "fwr3_nt",
    "CDR1": "cdr1_nt",
    "CDR2": "cdr2_nt",
    "CDR3": "cdr3_nt",
    "V-REGION": "v_gene_nt",
    "D-REGION": "d_gene_nt",
    "J-REGION": "j_gene_nt",
}
RENAME_DICT_TRANSLATE = {
    "cdr1_nt": "cdr1_aa",
    "fwr1_nt": "fwr1_aa",
    "fwr2_nt": "fwr2_aa",
    "fwr3_nt": "fwr3_aa",
    "cdr2_nt": "cdr2_aa",
    "cdr3_nt": "cdr3_aa",
}
