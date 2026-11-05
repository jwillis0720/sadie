"""
Reference YAML Auto-Generation
===============================

Generates a reference YAML config by combining:
1. Curated IMGT alleles from reference.g3.yml (for species it covers)
2. ALL IMGT alleles for species NOT in reference.g3.yml
3. ALL OGRDB alleles for every species with OGRDB data
4. ALL VDJbase alleles for every species with VDJbase data

The output follows the existing YAML format:
    name → provider → species → [allele_list]

Usage
-----
    from sadie.reference.generate import generate_reference_yaml
    config = generate_reference_yaml()

    # Or write directly to file:
    generate_reference_yaml(output_path=Path("reference.yml"))
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import yaml

from sadie.germlines.manager import GermlineManager

logger = logging.getLogger(__name__)

# Default path for the curated g3 YAML baseline
_DEFAULT_G3_PATH = Path(__file__).parent / "data" / "reference.g3.yml"

# YAML dict type alias
YamlConfig = Dict[str, Dict[str, Dict[str, List[str]]]]


def get_g3_curated_species(g3_yaml_path: Optional[Path] = None) -> YamlConfig:
    """Load the reference.g3.yml baseline and return its structure.

    Parameters
    ----------
    g3_yaml_path : Path, optional
        Path to the reference.g3.yml file. Defaults to repo root.

    Returns
    -------
    dict
        The parsed YAML data: name → provider → species → [allele_list]

    Raises
    ------
    FileNotFoundError
        If the g3 YAML file does not exist.
    """
    if g3_yaml_path is None:
        g3_yaml_path = _DEFAULT_G3_PATH

    if not g3_yaml_path.exists():
        raise FileNotFoundError(f"Curated g3 YAML not found at {g3_yaml_path}")

    with open(g3_yaml_path) as f:
        data: YamlConfig = yaml.safe_load(f)

    logger.info(f"Loaded g3 curated baseline with {len(data)} reference names from {g3_yaml_path}")
    return data


def _get_all_alleles_for_provider(provider_name: str, species: str) -> List[str]:
    """Get all allele names for a species from a specific provider.

    Parameters
    ----------
    provider_name : str
        Provider name (e.g., "imgt", "ogrdb", "vdjbase")
    species : str
        Species name (e.g., "human", "mouse")

    Returns
    -------
    list of str
        Sorted list of allele names
    """
    manager = GermlineManager(providers=[provider_name])
    alleles: Set[str] = set()

    for segment in ["V", "D", "J", "C"]:
        for chain in ["H", "K", "L"]:
            genes = manager.get_genes(species, segment, chain, functional_only=False)
            alleles.update(g.name for g in genes)

    return sorted(alleles)


def _collect_g3_species(g3_data: YamlConfig) -> Set[str]:
    """Collect all species that appear as sub-keys in the g3 curated data.

    This includes species used across all reference names and providers.

    Parameters
    ----------
    g3_data : dict
        The g3 YAML data.

    Returns
    -------
    set of str
        All species names used in g3 data.
    """
    species_set: Set[str] = set()
    for name in g3_data:
        for source in g3_data[name]:
            for species in g3_data[name][source]:
                species_set.add(species)
    return species_set


def generate_reference_yaml(
    g3_yaml_path: Optional[Path] = None,
    output_path: Optional[Path] = None,
) -> YamlConfig:
    """Generate a reference YAML config from germline sources.

    Generation logic:
    1. Load reference.g3.yml as the baseline for curated species.
       All reference names and their IMGT allele lists are preserved as-is.
    2. For species available in IMGT but NOT covered as reference names in
       reference.g3.yml, create new reference entries with ALL IMGT alleles.
    3. Add ALL OGRDB alleles for every species with OGRDB data.
    4. Add ALL VDJbase alleles for every species with VDJbase data.

    Parameters
    ----------
    g3_yaml_path : Path, optional
        Path to the reference.g3.yml baseline. Defaults to repo root.
    output_path : Path, optional
        If provided, write the YAML to this file path.

    Returns
    -------
    dict
        The generated YAML config: name → provider → species → [allele_list]
    """
    # Step 1: Load g3 curated baseline
    g3_data = get_g3_curated_species(g3_yaml_path)

    # Start with the g3 curated data as the foundation
    result: YamlConfig = {}
    for name in g3_data:
        result[name] = {}
        for source in g3_data[name]:
            result[name][source] = {}
            for species in g3_data[name][source]:
                result[name][source][species] = list(g3_data[name][source][species])

    # Collect reference names and species already covered by g3
    g3_ref_names: Set[str] = set(g3_data.keys())
    g3_species_used: Set[str] = _collect_g3_species(g3_data)

    logger.info(f"G3 curated reference names: {sorted(g3_ref_names)}")
    logger.info(f"G3 species used: {sorted(g3_species_used)}")

    # Step 2: Add ALL IMGT alleles for uncurated species
    imgt_manager = GermlineManager(providers=["imgt"])
    all_imgt_species = set(imgt_manager.get_available_species())

    # Species that have IMGT data but are NOT reference names in g3
    # AND are not already used as sub-species in g3 references
    uncurated_species = all_imgt_species - g3_ref_names - g3_species_used

    for species in sorted(uncurated_species):
        alleles = _get_all_alleles_for_provider("imgt", species)
        if alleles:
            # Create a new reference entry using species name as reference name
            if species not in result:
                result[species] = {}
            if "imgt" not in result[species]:
                result[species]["imgt"] = {}
            result[species]["imgt"][species] = alleles
            logger.info(f"Added {len(alleles)} IMGT alleles for uncurated species '{species}'")

    # Step 3: Add ALL OGRDB alleles
    ogrdb_manager = GermlineManager(providers=["ogrdb"])
    ogrdb_species = ogrdb_manager.get_available_species()

    for species in sorted(ogrdb_species):
        alleles = _get_all_alleles_for_provider("ogrdb", species)
        if not alleles:
            continue

        # Find the reference name for this species
        ref_name = _find_ref_name_for_species(result, species)
        if ref_name is None:
            # No existing reference for this species - create one
            ref_name = species
            if ref_name not in result:
                result[ref_name] = {}

        if "ogrdb" not in result[ref_name]:
            result[ref_name]["ogrdb"] = {}
        result[ref_name]["ogrdb"][species] = alleles
        logger.info(f"Added {len(alleles)} OGRDB alleles for species '{species}' under reference '{ref_name}'")

    # Step 4: Add ALL VDJbase alleles
    vdjbase_manager = GermlineManager(providers=["vdjbase"])
    vdjbase_species = vdjbase_manager.get_available_species()

    for species in sorted(vdjbase_species):
        alleles = _get_all_alleles_for_provider("vdjbase", species)
        if not alleles:
            continue

        # Find the reference name for this species
        ref_name = _find_ref_name_for_species(result, species)
        if ref_name is None:
            ref_name = species
            if ref_name not in result:
                result[ref_name] = {}

        if "vdjbase" not in result[ref_name]:
            result[ref_name]["vdjbase"] = {}
        result[ref_name]["vdjbase"][species] = alleles
        logger.info(f"Added {len(alleles)} VDJbase alleles for species '{species}' under reference '{ref_name}'")

    # Write to file if output path provided
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            yaml.dump(result, f, default_flow_style=False, sort_keys=False)
        logger.info(f"Wrote generated reference YAML to {output_path}")

    _log_summary(result)
    return result


def _find_ref_name_for_species(config: YamlConfig, species: str) -> Optional[str]:
    """Find the reference name that already has data for a given species.

    Searches through the existing config to find a reference name that
    contains the species in any provider section. If the species itself
    is a reference name, return it.

    Parameters
    ----------
    config : dict
        Current reference config.
    species : str
        Species to find.

    Returns
    -------
    str or None
        Reference name if found, None otherwise.
    """
    # First, check if species name is directly a reference name
    if species in config:
        return species

    # Search through all reference names for the species
    for ref_name in config:
        for source in config[ref_name]:
            if species in config[ref_name][source]:
                return ref_name

    return None


def _log_summary(config: YamlConfig) -> None:
    """Log a summary of the generated reference config.

    Parameters
    ----------
    config : dict
        The generated YAML config.
    """
    total_alleles = 0
    for name in config:
        for source in config[name]:
            for species in config[name][source]:
                count = len(config[name][source][species])
                total_alleles += count
                logger.debug(f"  {name}/{source}/{species}: {count} alleles")

    logger.info(f"Generated reference YAML: {len(config)} references, {total_alleles} total alleles")
