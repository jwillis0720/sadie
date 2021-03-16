from .chain import (
    AntibodyChainAA,
    AntibodyChainNT,
    HeavyChainAA,
    HeavyChainNT,
    KappaChainAA,
    KappaChainNT,
    LambdaChainAA,
    LambdaChainNT,
)

from .antibody import AntibodyAA, AntibodyNT
from .genetable import VGene, JGene

__all__ = [
    "AntibodyAA",
    "AntibodyChainAA",
    "AntibodyChainNT",
    "HeavyChainAA",
    "HeavyChainNT",
    "KappaChainAA",
    "KappaChainNT",
    "LambdaChainAA",
    "LambdaChainNT",
    "AntibodyAA",
    "AntibodyNT",
    "VGene",
    "JGene",
]
