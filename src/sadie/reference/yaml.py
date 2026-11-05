"""yaml object for reference data"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Set, Type

import pandas as pd
from yaml import load

logger = logging.getLogger("Reference")

try:
    from yaml import CLoader, Loader

    cload: Type[CLoader] | Type[Loader] = CLoader
except ImportError:
    from yaml import Loader

    cload = Loader


class YamlRef:
    def __init__(self, filepath: None | Path | str = None):
        if not filepath:
            filepath = Path(__file__).parent.joinpath("data/reference.yml")
        self.ref_path = filepath
        self.yaml = load(open(self.ref_path), Loader=Loader)
        self.yaml_df = self._normalize_and_verify_yaml()

    def get_names(self) -> Set[str]:
        """Return a list of all names whcih can be annotated

        Example
        -------
        yaml.get_names
        >>> ['human','mouse','macque']

        Returns
        -------
        set
            unique types
        """
        return set(self.yaml_df["name"].to_list())

    def get_genes(self, name: str, source: str, species: str) -> List[str]:
        """Get the genes associated with a name, source, and species

        Parameters
        ----------
        name : str
            ex. 'human'
        source: str
            ex. 'imgt'
        species : str
            ex.'human'

        Returns
        -------
        list
            list of genes

        Examples
        --------
        # get all annotated class from a human imgt genes
        object.get_genes('human','imgt','human')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        _a: List[str] = self.yaml.get(name).get(source).get(species)
        return _a

    def get_gene_segment(self, name: str, source: str, species: str, gene_segment: str) -> List[str]:
        """Get the genes associated with these keys

        Parameters
        ----------
        name : str
            ex. 'human'
        source: str
            ex. 'imgt'
        species : str
            ex.'human'
        gene_segment: str
            ex. V

        Returns
        -------
        list
            list of genes of the gene segment

        Examples
        --------
        object.get_gene_segment('human','imgt','human','V')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        return list(filter(lambda x: x[3] == gene_segment, self.get_genes(name, source, species)))

    def get_yaml_as_dataframe(self) -> pd.DataFrame:
        """Return yaml as a normalized dataframe"""
        return self.yaml_df

    def __repr__(self) -> Any:
        return self.yaml.__repr__()

    def __iter__(self) -> Any:
        """Iter method will step through the yaml file"""
        return self.yaml.__iter__()

    def __getitem__(self, key: str) -> Any:
        return self.yaml[key]

    def __len__(self) -> int:
        return len(self.yaml_df)

    def _normalize_and_verify_yaml(self) -> pd.DataFrame:
        dataframe_loader: List[Dict[str, str]] = []
        data = self.yaml
        for name in data:
            for source in data.get(name):
                for species in data.get(name).get(source):
                    dataframe_loader.append(
                        {
                            "name": name,
                            "source": source,
                            "species": species,
                            "genes": data.get(name).get(source).get(species),
                        }
                    )

        _df = pd.DataFrame(dataframe_loader).explode("genes").reset_index(drop=True)
        lookup: List[str] = ["name", "species", "genes"]
        duplicated = _df.set_index(lookup).loc[_df.groupby(lookup).size() > 1]  # type: ignore
        if not duplicated.empty:
            if len(duplicated["source"].unique()) == 1:
                # Same gene duplicated within a single source - this is an error
                raise ValueError(f"{duplicated}\nappears twice")
            else:
                # Cross-provider duplicates: keep first-come-first-serve, log and deduplicate
                dup_genes = duplicated.reset_index()["genes"].unique()
                for gene in dup_genes:
                    logger.debug(f"Skipping duplicate allele '{gene}' from later provider (first-come-first-serve)")
                # Drop duplicate (name, species, genes) rows, keeping the first occurrence
                _df = _df.drop_duplicates(subset=["name", "species", "genes"], keep="first").reset_index(drop=True)
                # Also update self.yaml to reflect deduplication so from_yaml sees deduped gene lists
                self._deduplicate_yaml_data(data)
        return _df

    def _deduplicate_yaml_data(self, data: Dict[str, Any]) -> None:
        """Remove cross-provider duplicate alleles from the raw YAML data in-place.

        For each reference name, tracks seen allele names per (name, species) and removes
        duplicates from later providers, preserving first-come-first-serve ordering.

        Parameters
        ----------
        data : dict
            The raw YAML data dictionary (name → source → species → [genes]).
        """
        for name in data:
            # Track seen genes per (name, species) key
            seen: Dict[str, Set[str]] = {}
            for source in data.get(name, {}):
                for species in data.get(name, {}).get(source, {}):
                    if species not in seen:
                        seen[species] = set()
                    gene_list: List[str] = data[name][source][species]
                    deduped: List[str] = []
                    for gene in gene_list:
                        if gene not in seen[species]:
                            seen[species].add(gene)
                            deduped.append(gene)
                        else:
                            logger.debug(
                                f"Duplicate allele '{gene}' in '{source}' for '{species}' "
                                f"under '{name}' skipped (already seen from earlier provider)"
                            )
                    data[name][source][species] = deduped
