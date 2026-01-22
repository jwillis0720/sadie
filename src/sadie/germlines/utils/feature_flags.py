"""
Feature flag utilities for germlines module migration.

This module provides feature flag controls for the gradual migration from G3 API
to the local germlines module, following Constitution Principle V (Integration Compatibility).
"""

import os
import logging
from functools import lru_cache

logger = logging.getLogger(__name__)

SADIE_USE_GERMLINES_MODULE = "SADIE_USE_GERMLINES_MODULE"


@lru_cache(maxsize=1)
def use_germlines_module() -> bool:
    """
    Determine whether to use the germlines module or fall back to G3 API.

    This feature flag controls the migration strategy from G3 API to the local
    germlines module. During the validation period, users can toggle between
    implementations for testing and rollback capability.

    Environment Variable:
        SADIE_USE_GERMLINES_MODULE: Set to "true" (default) or "false"

    Returns:
        bool: True to use germlines module, False to use G3 API

    Behavior:
        - Default: "true" (use germlines module)
        - "true" -> Use germlines module exclusively (no G3 API calls)
        - "false" -> Use G3 API exclusively (for rollback during validation period only)

    Notes:
        - No automatic fallback implemented
        - G3 fallback is only for validation period
        - After validation period, G3 dependencies will be removed
    """
    env_value = os.environ.get(SADIE_USE_GERMLINES_MODULE, "true").lower()
    use_germlines = env_value in ("true", "1", "yes", "on")
    if not use_germlines:
        logger.warning(
            "G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true. "
            "G3 will be removed after 2026-06-01."
        )
    logger.debug(
        f"Feature flag SADIE_USE_GERMLINES_MODULE={env_value} -> "
        f"use_germlines={use_germlines}"
    )
    return use_germlines


def clear_feature_flag_cache():
    """Clear the cached feature flag value to allow re-evaluation."""
    use_germlines_module.cache_clear()
