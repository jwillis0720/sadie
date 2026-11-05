from sadie.reference.cache import DatabaseCache, compute_cache_key
from sadie.reference.generate import generate_reference_yaml, get_g3_curated_species
from sadie.reference.reference import G3Error, Reference, References
from sadie.reference.settings import (
    BLAST_CONVENTION,
    IMGT_GB_LOOKUP,
    IMGT_LOOKUP,
    J_SEGMENTS,
    MOTIF_LOOKUP,
    get_motif_provenance,
    RECEPTORS,
    RENAME_DICT,
    RENAME_DICT_TRANSLATE,
    REVERSE_IMGT_LOOKUP,
    SEGMENTS,
    SEGMENTS_INTERNAL_DATA,
    IMGT_DEF_nt,
)

# Create MOTIF_PROVENANCE alias for the public API
MOTIF_PROVENANCE = get_motif_provenance()
from sadie.reference.yaml import YamlRef

__all__ = [
    "DatabaseCache",
    "compute_cache_key",
    "References",
    "Reference",
    "RECEPTORS",
    "IMGT_DEF_nt",
    "IMGT_GB_LOOKUP",
    "IMGT_LOOKUP",
    "G3Error",
    "BLAST_CONVENTION",
    "J_SEGMENTS",
    "MOTIF_LOOKUP",
    "MOTIF_PROVENANCE",
    "RENAME_DICT",
    "RENAME_DICT_TRANSLATE",
    "REVERSE_IMGT_LOOKUP",
    "SEGMENTS",
    "SEGMENTS_INTERNAL_DATA",
    "YamlRef",
    "generate_reference_yaml",
    "get_g3_curated_species",
]
