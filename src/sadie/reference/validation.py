"""
Validation utilities for MOTIF_LOOKUP data and schema validation.

This module provides utilities for validating JSON structure, motif patterns,
and data integrity for the MOTIF_LOOKUP registry.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import logging

try:
    import jsonschema
    from jsonschema import Draft7Validator, ValidationError
    JSONSCHEMA_AVAILABLE = True
except ImportError:
    JSONSCHEMA_AVAILABLE = False

logger = logging.getLogger(__name__)


class MotifValidationError(Exception):
    """Custom exception for motif validation errors."""
    pass


class JSONValidationError(Exception):
    """Custom exception for JSON schema validation errors."""
    pass


def validate_motif_pattern(pattern: str) -> bool:
    """
    Validate that a motif pattern is a valid regex.

    Args:
        pattern: Regex pattern string to validate

    Returns:
        True if pattern is valid, False otherwise

    Raises:
        MotifValidationError: If pattern is invalid
    """
    try:
        re.compile(pattern)
        return True
    except re.error as e:
        raise MotifValidationError(f"Invalid regex pattern '{pattern}': {e}")


def validate_species_name(species: str) -> bool:
    """
    Validate species name format (lowercase with underscores and numbers).

    Args:
        species: Species name to validate

    Returns:
        True if species name is valid

    Raises:
        MotifValidationError: If species name format is invalid
    """
    pattern = r'^[a-z0-9_]+$'
    if not re.match(pattern, species):
        raise MotifValidationError(
            f"Invalid species name '{species}'. Must match pattern: {pattern}"
        )
    return True


def validate_locus_name(locus: str) -> bool:
    """
    Validate locus name format.

    Args:
        locus: Locus name to validate (e.g., IGHJ, IGKJ, IGLJ)

    Returns:
        True if locus name is valid

    Raises:
        MotifValidationError: If locus name format is invalid
    """
    pattern = r'^(IGH|IGK|IGL|TRA|TRB|TRD|TRG)[JV]$'
    if not re.match(pattern, locus):
        raise MotifValidationError(
            f"Invalid locus name '{locus}'. Must match pattern: {pattern}"
        )
    return True


def validate_motif_data_structure(data: Dict[str, Any]) -> bool:
    """
    Validate the basic structure of motif data.

    Args:
        data: Motif data dictionary to validate

    Returns:
        True if structure is valid

    Raises:
        MotifValidationError: If data structure is invalid
    """
    if not isinstance(data, dict):
        raise MotifValidationError("Motif data must be a dictionary")

    for species, species_data in data.items():
        # Validate species name
        validate_species_name(species)

        if not isinstance(species_data, dict):
            raise MotifValidationError(f"Species data for '{species}' must be a dictionary")

        for locus, pattern in species_data.items():
            if locus == "_provenance":
                # Skip provenance validation here, will be handled separately
                continue
            elif locus == "ignore":
                # Validate ignore list
                if not isinstance(pattern, list):
                    raise MotifValidationError(
                        f"Ignore list for '{species}' must be a list, got {type(pattern)}"
                    )
            else:
                # Validate locus and pattern
                validate_locus_name(locus)
                if not isinstance(pattern, str):
                    raise MotifValidationError(
                        f"Pattern for '{species}.{locus}' must be a string, got {type(pattern)}"
                    )
                validate_motif_pattern(pattern)

    return True


def validate_provenance_metadata(provenance: Dict[str, Any]) -> bool:
    """
    Validate provenance metadata structure.

    Args:
        provenance: Provenance metadata dictionary

    Returns:
        True if provenance metadata is valid

    Raises:
        MotifValidationError: If provenance metadata is invalid
    """
    required_fields = {"source", "imgt_validated", "last_reviewed"}
    optional_fields = {"notes", "validation_count", "validation_rate"}

    if not isinstance(provenance, dict):
        raise MotifValidationError("Provenance metadata must be a dictionary")

    # Check required fields
    missing_fields = required_fields - set(provenance.keys())
    if missing_fields:
        raise MotifValidationError(f"Missing required provenance fields: {missing_fields}")

    # Check field types
    if not isinstance(provenance["source"], str) or not provenance["source"]:
        raise MotifValidationError("Provenance 'source' must be a non-empty string")

    if not isinstance(provenance["imgt_validated"], bool):
        raise MotifValidationError("Provenance 'imgt_validated' must be a boolean")

    # Validate date format (YYYY-MM-DD)
    date_pattern = r'^\d{4}-\d{2}-\d{2}$'
    if not re.match(date_pattern, provenance["last_reviewed"]):
        raise MotifValidationError(
            f"Provenance 'last_reviewed' must be in YYYY-MM-DD format, "
            f"got '{provenance['last_reviewed']}'"
        )

    # Validate optional fields if present
    if "validation_count" in provenance:
        if not isinstance(provenance["validation_count"], int) or provenance["validation_count"] < 0:
            raise MotifValidationError("Provenance 'validation_count' must be a non-negative integer")

    if "validation_rate" in provenance:
        rate = provenance["validation_rate"]
        if not isinstance(rate, (int, float)) or not (0 <= rate <= 1):
            raise MotifValidationError("Provenance 'validation_rate' must be a number between 0 and 1")

    return True


def load_json_schema(schema_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load JSON schema from file.

    Args:
        schema_path: Path to JSON schema file

    Returns:
        Loaded schema dictionary

    Raises:
        JSONValidationError: If schema cannot be loaded
    """
    try:
        with open(schema_path, 'r', encoding='utf-8') as f:
            schema = json.load(f)
        return schema
    except (FileNotFoundError, PermissionError) as e:
        raise JSONValidationError(f"Cannot load schema from '{schema_path}': {e}")
    except json.JSONDecodeError as e:
        raise JSONValidationError(f"Invalid JSON in schema '{schema_path}': {e}")


def validate_against_schema(data: Dict[str, Any], schema: Dict[str, Any]) -> bool:
    """
    Validate data against JSON schema.

    Args:
        data: Data to validate
        schema: JSON schema to validate against

    Returns:
        True if data is valid

    Raises:
        JSONValidationError: If validation fails or jsonschema not available
    """
    if not JSONSCHEMA_AVAILABLE:
        raise JSONValidationError(
            "jsonschema package not available. Install with: pip install jsonschema"
        )

    try:
        validator = Draft7Validator(schema)
        errors = list(validator.iter_errors(data))

        if errors:
            error_messages = []
            for error in errors:
                path = ".".join(str(p) for p in error.path)
                error_messages.append(f"Path '{path}': {error.message}")

            raise JSONValidationError(
                f"Schema validation failed:\n" + "\n".join(error_messages)
            )

        return True
    except ValidationError as e:
        raise JSONValidationError(f"Schema validation error: {e}")


def comprehensive_motif_validation(
    data: Dict[str, Any],
    schema_path: Optional[Union[str, Path]] = None
) -> bool:
    """
    Perform comprehensive validation of motif data.

    Args:
        data: Motif data to validate
        schema_path: Optional path to JSON schema for validation

    Returns:
        True if all validations pass

    Raises:
        MotifValidationError, JSONValidationError: If any validation fails
    """
    logger.info("Starting comprehensive motif data validation")

    # Validate basic data structure
    validate_motif_data_structure(data)
    logger.info("✓ Data structure validation passed")

    # Validate provenance metadata for each species
    provenance_count = 0
    for species, species_data in data.items():
        if "_provenance" in species_data:
            validate_provenance_metadata(species_data["_provenance"])
            provenance_count += 1

    logger.info(f"✓ Provenance metadata validation passed for {provenance_count} species")

    # Schema validation if schema provided
    if schema_path:
        schema = load_json_schema(schema_path)
        validate_against_schema(data, schema)
        logger.info(f"✓ JSON schema validation passed against {schema_path}")

    logger.info("✓ Comprehensive validation completed successfully")
    return True


def quick_validation_summary(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate a quick validation summary without raising exceptions.

    Args:
        data: Motif data to analyze

    Returns:
        Summary dictionary with validation statistics
    """
    summary = {
        "total_species": len(data),
        "species_with_provenance": 0,
        "imgt_validated_species": 0,
        "legacy_species": 0,
        "total_patterns": 0,
        "validation_errors": []
    }

    for species, species_data in data.items():
        try:
            validate_species_name(species)

            if "_provenance" in species_data:
                summary["species_with_provenance"] += 1
                provenance = species_data["_provenance"]

                try:
                    validate_provenance_metadata(provenance)
                    if provenance.get("imgt_validated", False):
                        summary["imgt_validated_species"] += 1
                    else:
                        summary["legacy_species"] += 1
                except MotifValidationError as e:
                    summary["validation_errors"].append(f"{species}._provenance: {e}")

            # Count patterns
            for key, value in species_data.items():
                if key not in ("_provenance", "ignore") and isinstance(value, str):
                    summary["total_patterns"] += 1
                    try:
                        validate_motif_pattern(value)
                    except MotifValidationError as e:
                        summary["validation_errors"].append(f"{species}.{key}: {e}")

        except MotifValidationError as e:
            summary["validation_errors"].append(f"{species}: {e}")

    summary["validation_passed"] = len(summary["validation_errors"]) == 0
    return summary