from yaml import load
from pathlib import Path

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


class YamlRef:
    def __init__(self, filepath=None):

        if not filepath:
            self.ref_path = Path(__file__).parent.joinpath("data/reference.yml")
        else:
            self.ref_path = filepath

        self.yaml = load(open(self.ref_path), Loader=Loader)

    def get_reference_types(self) -> set:
        """Return reference type. ex imgt or custom

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

    def get_genes(self, reference_type: str, functional: str, species: str, sub_species: str) -> list:
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
        object.get_genes('imgt','functional','human','human')
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


if __name__ == "__main__":
    y = YamlRef()
    print(y)
