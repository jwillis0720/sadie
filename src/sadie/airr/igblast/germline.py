import os
import warnings
from pathlib import Path

# package/module level
from sadie.reference.yaml import YamlRef
from sadie.airr.igblast.igblast import ensure_prefix_to


class GermlineData:
    """
    The germline data paths are extremely cumbersome to workwith. This class will abstract away their paths to make it easier to fold into IgBLAST

    Examples
    --------
    >>> gd = GermlineData('human')
    >>> gd.base_dir
    /Users/jwillis/repos/sadie/airr/data/germlines
    >>> gd.v_gene_dir
    /Users/jwillis/repos/sadie/airr/data/germlines/blastdb/Ig/human/human_V'
    >>> gd.aux_path
    /Users/jwillis/repos/sadie/airr/data/germlines/aux_data/human_gl.aux
    """

    def __init__(
        self,
        species: str,
        database: str = "imgt",
        receptor: str = "Ig",
    ):
        """

        Parameters
        ----------
        species : str
            The species of interest, e.g. human
        receptor : str, optional
            the receptor type, by default "Ig"
        """
        self.species = species
        self.base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data/germlines"))
        self.blast_dir = os.path.join(self.base_dir, f"{database}/{receptor}/blastdb/{species}_")
        self.v_gene_dir = self.blast_dir + "V"
        self.d_gene_dir = self.blast_dir + "D"
        self.j_gene_dir = self.blast_dir + "J"
        self.aux_path = os.path.join(self.base_dir, f"{database}/aux_db/{species}_gl.aux")
        self.igdata = os.path.join(self.base_dir, f"{database}/{receptor}/")

    @property
    def base_dir(self) -> Path:
        """The base dir

        Returns
        -------
        Path
            The base directory path that contains all the germline data
        """
        return self._base_dir

    @base_dir.setter
    def base_dir(self, directory: str):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"Base directory, {directory} not found")
        self._base_dir = directory

    @property
    def blast_dir(self) -> Path:
        return self._blast_dir

    @blast_dir.setter
    def blast_dir(self, directory: str):
        # Must be a parent since this is not a valid path yet
        _path = Path(directory).parent
        if not _path.exists():
            raise FileNotFoundError(f"Blast directory, {directory} not found")
        self._blast_dir = directory

    @property
    def v_gene_dir(self) -> Path:
        """The V gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           human_V does not exists but it's the prefix to human_V.nod and other files used by blast
        """
        return self._v_gene_dir

    @v_gene_dir.setter
    def v_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            raise FileNotFoundError(f"V gene directory glob, {directory} not found")
        self._v_gene_dir = _path

    @property
    def d_gene_dir(self) -> Path:
        """The D gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex: human_D does not exists but it's the prefix to human_D.nod and other files used by blast
        """
        return self._d_gene_dir

    @d_gene_dir.setter
    def d_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            warnings.warn(f"D gene directory not found for {self.species}", UserWarning)
        self._d_gene_dir = _path

    @property
    def j_gene_dir(self) -> Path:
        """The J gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex: human_J does not exists but it's the prefix to human_j.nod and other files used by blast
        """
        return self._j_gene_dir

    @j_gene_dir.setter
    def j_gene_dir(self, directory: str):
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._j_gene_dir = _path

    @property
    def aux_path(self) -> Path:
        """The auxillary data path used to reconstruct CDR3 regions.

        Returns
        -------
        Path
           the fully qualified path to the species auxilary data
           ex:/Users/jwillis/repos/sadie/airr/data/germlines/aux_data/human_gl.aux
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, directory: str):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._aux_path = _path

    @property
    def igdata(self) -> Path:
        return self._igdata

    @igdata.setter
    def igdata(self, directory: Path):
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"IGDATA, {directory} not found")
        self._igdata = _path

    @staticmethod
    def get_available_datasets() -> list:
        """A static non-instantiated method to get a list of avaialble species with the builtin data

        Returns
        -------
        list
           available datasets (common_name, custom|imgt, functional|all)
        """
        y = YamlRef()
        db_types = []
        for database_type in y.yaml:
            for common in y.yaml[database_type]:
                if (common, database_type) not in db_types:
                    db_types.append((common, database_type))
        return db_types
