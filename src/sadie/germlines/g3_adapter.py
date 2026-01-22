"""
G3 API Adapter
==============

Transforms GermlineGene objects from the germlines module into G3 API response format.

This adapter allows existing code that expects G3 API responses to work seamlessly
with the local germlines database.

Design:
- Input: GermlineGene objects from germlines module
- Output: Dict matching G3 API JSON structure
- Handles region extraction and formatting
"""

import logging
from typing import Dict, List, Any
from sadie.germlines.models import GermlineGene

logger = logging.getLogger(__name__)


class GermlineToG3Adapter:
    """
    Adapter to transform GermlineGene objects to G3 API format.

    The G3 API returns germline data as JSON with nested IMGT fields.
    This adapter converts GermlineGene objects to match that structure.

    Examples
    --------
    >>> from sadie.germlines import get_gene_by_name
    >>> adapter = GermlineToG3Adapter()
    >>> gene = get_gene_by_name("human", "IGHV1-69*01")
    >>> g3_dict = adapter.to_g3_format(gene)
    >>> g3_dict["gene"]
    'IGHV1-69*01'
    >>> g3_dict["imgt"]["imgt_functional"]
    'F'
    """

    def to_g3_format(self, gene: GermlineGene) -> Dict[str, Any]:
        """
        Convert GermlineGene to G3 API response format.

        Parameters
        ----------
        gene : GermlineGene
            Germline gene object from germlines module

        Returns
        -------
        Dict[str, Any]
            Dictionary matching G3 API JSON structure

        Examples
        --------
        >>> adapter = GermlineToG3Adapter()
        >>> gene = GermlineGene(
        ...     name="IGHV1-69*01",
        ...     species="human",
        ...     segment="V",
        ...     chain="H",
        ...     sequence="CAGGTG...",
        ...     functionality="F",
        ...     source="imgt"
        ... )
        >>> g3_dict = adapter.to_g3_format(gene)
        """
        # Extract chain from gene name (e.g., "IGHV1-69*01" -> "H")
        chain = self._extract_chain_from_name(gene.name, gene.chain)

        # Build base structure
        g3_dict = {
            "source": gene.source,
            "common": gene.species,
            "latin": self._get_latin_name(gene.species),
            "gene": gene.name,
            "label": f"{gene.segment}-REGION",
            "gene_segment": gene.segment,
            "receptor": "IG",  # Immunoglobulin (vs TCR for T-cell receptors)
            "sequence": gene.sequence,
            "species": gene.species,  # Added for compatibility
        }

        # Build IMGT nested structure
        imgt_dict = {
            "sequence": gene.sequence,
            "sequence_gapped": gene.sequence_gapped or "",
            "sequence_gapped_aa": gene.sequence_aa_gapped or "",
            "imgt_functional": gene.functionality,
            "contrived_functional": gene.functionality,
        }

        # Add regions if available
        if gene.regions and gene.region_positions:
            self._add_regions_to_imgt(imgt_dict, gene)

        g3_dict["imgt"] = imgt_dict

        return g3_dict

    def to_g3_format_batch(self, genes: List[GermlineGene]) -> List[Dict[str, Any]]:
        """
        Convert multiple GermlineGene objects to G3 format.

        Parameters
        ----------
        genes : List[GermlineGene]
            List of germline genes

        Returns
        -------
        List[Dict[str, Any]]
            List of G3-formatted dictionaries
        """
        return [self.to_g3_format(gene) for gene in genes]

    def _extract_chain_from_name(self, gene_name: str, default_chain: str) -> str:
        """
        Extract chain type from gene name.

        Parameters
        ----------
        gene_name : str
            Gene name (e.g., "IGHV1-69*01")
        default_chain : str
            Default chain if extraction fails

        Returns
        -------
        str
            Chain type ("H", "K", or "L")
        """
        # Gene names follow pattern: IG{chain}{segment}...
        # e.g., IGHV1-69*01 -> H, IGKV1-39*01 -> K, IGLV1-40*01 -> L
        if len(gene_name) >= 3 and gene_name.startswith("IG"):
            chain_char = gene_name[2].upper()
            if chain_char in ["H", "K", "L"]:
                return chain_char

        return default_chain

    def _get_latin_name(self, common_name: str) -> str:
        """
        Get Latin (scientific) name for species.

        Parameters
        ----------
        common_name : str
            Common name (e.g., "human")

        Returns
        -------
        str
            Latin name with underscores (e.g., "Homo_sapiens")
        """
        latin_names = {
            "human": "Homo_sapiens",
            "mouse": "Mus_musculus",
            "rat": "Rattus_norvegicus",
            "rabbit": "Oryctolagus_cuniculus",
            "dog": "Canis_lupus_familiaris",
            "cat": "Felis_catus",
            "macaque": "Macaca_mulatta",
            "alpaca": "Vicugna_pacos",
        }
        return latin_names.get(common_name.lower(), common_name.capitalize())

    def _add_regions_to_imgt(self, imgt_dict: Dict[str, Any], gene: GermlineGene) -> None:
        """
        Add CDR/FWR regions to IMGT dict.

        Modifies imgt_dict in place to add region sequences and positions.

        Parameters
        ----------
        imgt_dict : Dict[str, Any]
            IMGT dictionary to modify
        gene : GermlineGene
            Gene with regions and region_positions
        """
        if not gene.regions or not gene.region_positions:
            return

        # V genes have: FWR1, CDR1, FWR2, CDR2, FWR3, CDR3
        # J genes have: FWR4
        # D genes typically don't have regions

        for region_name in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]:
            # Get nucleotide sequence
            if region_name in gene.regions:
                imgt_dict[region_name] = gene.regions[region_name]

            # Get amino acid sequence (if available)
            aa_key = f"{region_name}_aa"
            if region_name in gene.regions and gene.sequence_aa:
                # Extract AA from full AA sequence using positions
                if region_name in gene.region_positions:
                    start, end = gene.region_positions[region_name]
                    # Convert nucleotide positions to amino acid positions
                    aa_start = start // 3
                    aa_end = end // 3
                    if gene.sequence_aa_gapped:
                        # Use gapped sequence if available
                        imgt_dict[aa_key] = gene.sequence_aa_gapped[aa_start:aa_end]
                    elif gene.sequence_aa:
                        imgt_dict[aa_key] = gene.sequence_aa[aa_start:aa_end]

            # Get positions
            if region_name in gene.region_positions:
                start, end = gene.region_positions[region_name]
                imgt_dict[f"{region_name}_start"] = start
                imgt_dict[f"{region_name}_end"] = end


def create_g3_adapter() -> GermlineToG3Adapter:
    """
    Factory function to create G3 adapter instance.

    Returns
    -------
    GermlineToG3Adapter
        Adapter instance
    """
    return GermlineToG3Adapter()
