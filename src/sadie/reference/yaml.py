from typing import Union
from yaml import load
from pathlib import Path

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


class YamlRef:
    def __init__(self, filepath: Union[None, Path, str] = None):
        if not filepath:
            self.ref_path = Path(__file__).parent.joinpath("data/reference.yml")
        else:
            self.ref_path = filepath

        self.yaml = load(open(self.ref_path), Loader=Loader)

    def get_database_types(self) -> set:
        """Return database types in current yaml

        Example
        -------
        yaml.get_reference_types()
        >>> {'imgt','custom'}

        Returns
        -------
        set
            unique types
        """
        return set(self.yaml.keys())

    def get_species_keys(self, reference_type: str) -> list:
        """Return functional types from reference types and functional type

        Example
        -------
        yaml.get_species_keys('imgt','functional')
        >>> {'alpaca','cat','dog'...}

        Returns
        -------
        set
            unique types
        """
        return list(self.yaml.get(reference_type).keys())

    def get_sub_species(self, reference_type: str, species: str) -> list:
        """
        get sub species keys. For instance for humanized mouse, a sub species will be mouse and human

        Parameters
        ----------
        reference_type : str
            ex. imgt
        functional : str
            ex. all
        species : str
            cat

        Returns
        -------
        int
            number of sub species
        """
        return list(self.yaml.get(reference_type).get(species).keys())

    def get_genes(self, reference_type: str, species: str, sub_species: str) -> list:
        """Get the genes associated with these keys

        Parameters
        ----------
        reference_type : str
            ex. imgt
        species : str
            ex. human
        sub_species : str
            ex. human

        Returns
        -------
        list
            list of genes

        Examples
        --------
        object.get_genes('imgt'','human','human')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        return self.yaml.get(reference_type).get(species).get(sub_species)

    def get_gene_segment(self, reference_type: str, species: str, sub_species: str, gene_segment: str) -> list:
        """Get the genes associated with these keys

        Parameters
        ----------
        reference_type : str
            ex. imgt
        species : str
            ex. human
        sub_species : str
            ex. human
        gene_segment: str
            ex. V
        Returns
        -------
        list
            list of genes of the gene segment

        Examples
        --------
        object.get_gene_segment('imgt','human','human','V')
        >>> ['IGHV1-2*01','IGHV1-69*01'....]
        """
        return list(filter(lambda x: x[3] == gene_segment, self.get_genes(reference_type, species, sub_species)))

    def __repr__(self):
        return self.yaml.__repr__()

    def __iter__(self):
        """Iter method will step through the yaml file"""
        for database in self.yaml:
            for species in self.yaml[database]:
                for sub_species in self.yaml[database][species]:
                    yield {
                        "database": database,
                        "species": species,
                        "sub_species": sub_species,
                        "gene": self.yaml[database][species][sub_species],
                    }
