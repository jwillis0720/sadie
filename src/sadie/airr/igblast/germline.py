from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Set

from sadie.airr.igblast.igblast import ensure_prefix_to
from sadie.reference import YamlRef

logger = logging.getLogger(__name__)


def _use_germlines_module() -> bool:
    import os
    env_value = os.environ.get("SADIE_USE_GERMLINES_MODULE", "true").lower()
    use_germlines = env_value in ("true", "1", "yes")
    if not use_germlines:
        logger.warning(
            "G3 API is deprecated. Set SADIE_USE_GERMLINES_MODULE=true. "
            "G3 will be removed after 2026-06-01."
        )
    return use_germlines


def _get_germlines_igblast_dir() -> Path:
    from sadie.germlines import get_germlines_base_dir
    return get_germlines_base_dir() / "igblast"


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
        name: str,
        receptor: str = "Ig",
        database_dir: Optional[str | Path] = None,
        scheme: str = "imgt",
    ):
        """

        Parameters
        ----------
        species : str
            The species of interest, e.g. human
        receptor : str, optional
            the receptor type, by default "Ig"
        """
        self.name = name

        # Determine base directory based on feature flag
        if database_dir:
            self.base_dir = Path(database_dir).absolute()
        elif _use_germlines_module():
            # Use germlines module paths (new default)
            germlines_igblast = _get_germlines_igblast_dir()
            internal_data_species = germlines_igblast / "Ig" / "internal_data" / name

            # Check if this species has databases in germlines module
            if internal_data_species.exists():
                self.base_dir = germlines_igblast
                # Germlines module uses IgBLAST-compatible directory structure:
                # igblast/Ig/internal_data/{species}/ contains both .ndm.imgt and BLAST databases
                self.blast_dir = internal_data_species / f"{name}_"
                self.v_gene_dir = Path(self.blast_dir.__str__() + "V")
                self.d_gene_dir = Path(self.blast_dir.__str__() + "D")
                self.j_gene_dir = Path(self.blast_dir.__str__() + "J")
                self.c_gene_dir = Path(self.blast_dir.__str__() + "C")
                self.aux_path = germlines_igblast / "aux_db" / f"{name}_gl.aux"
                # IGDATA points to the directory containing internal_data/
                self.igdata = germlines_igblast / "Ig"
            else:
                # NFR-002: No silent fallback to G3 when germlines is selected
                raise ValueError(
                    f"Species '{name}' not found in germlines module at "
                    f"{internal_data_species}. "
                    f"Build germlines databases with: update_databases('{name}'). "
                    f"To use legacy G3 paths, set SADIE_USE_GERMLINES_MODULE=false."
                )
        else:
            # Legacy G3 paths (deprecated, for backwards compatibility)
            self._use_legacy_paths(name, receptor, scheme)

    def _use_legacy_paths(self, name: str, receptor: str, scheme: str) -> None:
        """Set paths to legacy G3 directory structure.

        Parameters
        ----------
        name : str
            Species name
        receptor : str
            Receptor type (Ig or TCR)
        scheme : str
            Numbering scheme
        """
        self.base_dir = Path(__file__).absolute().parent / "../data/germlines/"
        self.blast_dir = Path(str(self.base_dir) + f"/{receptor}/blastdb/{name}/{name}_")
        self.v_gene_dir = Path(self.blast_dir.__str__() + "V")
        self.d_gene_dir = Path(self.blast_dir.__str__() + "D")
        self.j_gene_dir = Path(self.blast_dir.__str__() + "J")
        self.c_gene_dir = Path(self.blast_dir.__str__() + "C")
        self.aux_path = self.base_dir / f"aux_db/{scheme}/{name}_gl.aux"
        # the literal 'internal_data/{name}` must be discovered by IgBLAST
        self.igdata = self.base_dir / f"{receptor}/"

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
    def base_dir(self, directory: str | Path) -> None:
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"Base directory, {directory} not found")
        self._base_dir = _path

    @property
    def blast_dir(self) -> Path:
        return self._blast_dir

    @blast_dir.setter
    def blast_dir(self, directory: str | Path) -> None:
        # Must be a parent since this is not a valid path yet
        if not Path(directory).parent.exists():
            raise FileNotFoundError(f"Blast directory, {directory} not found")
        self._blast_dir = Path(directory)

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
    def v_gene_dir(self, directory: str | Path) -> None:
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
    def d_gene_dir(self, directory: str | Path) -> None:
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            warnings.warn(f"D gene directory not found for {self.name}", UserWarning)
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
    def j_gene_dir(self, directory: str | Path) -> None:
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._j_gene_dir = _path

    @property
    def c_gene_dir(self) -> Path:
        """The C gene directory prefix for the species of interest

        Returns
        -------
        str
           this is not a qualified path but a glob path.
           ex: human_C does not exists but it's the prefix to human_C.nod and other files used by blast
        """
        return self._c_gene_dir

    @c_gene_dir.setter
    def c_gene_dir(self, directory: str | Path) -> None:
        _path = Path(directory)
        if not ensure_prefix_to(_path):
            warnings.warn(f"C gene directory not found for {self.name}", UserWarning)
        self._c_gene_dir = _path

    @property
    def aux_path(self) -> Path:
        """The auxillary data path used to reconstruct CDR3 regions.

        Returns
        -------
        Path
           the fully qualified path to the species auxilary data
           ex:/Users/jwillis/repos/sadie/airr/data/germlines/aux_data/{scheme}/human_gl.aux
        """
        return self._aux_path

    @aux_path.setter
    def aux_path(self, directory: str | Path) -> None:
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"J gene directory glob, {directory} not found")
        self._aux_path = _path

    @property
    def igdata(self) -> Path:
        return self._igdata

    @igdata.setter
    def igdata(self, directory: Path) -> None:
        _path = Path(directory)
        if not _path.exists():
            raise FileNotFoundError(f"IGDATA, {directory} not found")
        self._igdata = _path

    @staticmethod
    def get_available_datasets() -> Set[str]:
        """A static non-instantiated method to get a list of available species with the builtin data

        Returns
        -------
        Set[str]
            Set of available species names
        """
        datasets: Set[str] = set()
        
        # Add germlines module species if feature flag is enabled
        if _use_germlines_module():
            germlines_internal_data = _get_germlines_igblast_dir() / "Ig" / "internal_data"
            if germlines_internal_data.exists():
                for species_dir in germlines_internal_data.iterdir():
                    if species_dir.is_dir() and not species_dir.name.startswith('.'):
                        datasets.add(species_dir.name)
        
        # Also add legacy G3 species for backwards compatibility
        y = YamlRef()
        datasets.update(y.get_names())
        
        return datasets
