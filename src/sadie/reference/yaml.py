"""yaml object for reference data"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Set, Type
import pandas as pd
from yaml import load

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
        duplicated = _df.set_index(lookup).loc[_df.groupby(lookup).size() > 1]
        if not duplicated.empty:
            if len(duplicated["source"].unique()) == 1:
                raise ValueError(f"{duplicated}\nappears twice")
            else:
                raise ValueError(f"{duplicated}\nappears twice from two difference sources")
        return _df
