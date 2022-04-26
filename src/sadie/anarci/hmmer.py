from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Union, List

# from Bio.SeqRecord import SeqRecord
import pyhmmer

# Sequence = Union[List[Path, SeqRecord, str], Path, SeqRecord, str]


class HMMER:
    """
    Extension of Pyhmmer with built-in alignment models that are specie specific.
    """
    
    def __init__(self):
        self.hmm_folder = Path(__file__).parent / "data/anarci/HMMs"
        self.hmm_paths = self.hmm_folder.glob('*.hmm')
        self.species_to_paths = self.get_species_to_paths()
        self.available_species = sorted(self.species_to_paths.keys())
        
    def get_species_to_paths(self) -> dict:
        """
        Return a dictionary with the available species and their paths.

        Returns
        -------
        dict
            Dictionary with the available species and their paths.
        """
        species_to_paths = {}
        for path in self.hmm_paths:
            species = path.stem.split('_')[0]
            try:
                species_to_paths[species].append(path)
            except KeyError:
                species_to_paths[species] = [path]
        return species_to_paths
                
    # @lru_cache(maxsize=None)
    def get_hmm_models(self, species: Union[List[str], str]) -> List[pyhmmer.plan7.HMMFile]:
        """
        Return a HMMER model for a given specie.

        Parameters
        ----------
        get_species_model: Union[list, str]
            Available species: ALL, alpaca, cat, cow, dog, human, mouse, pig, rabbit, rat, rhesus

        Returns
        -------
        pyhmmer.plan7.HMMERModel
            HMM model for a specific species 
        """
        hmms = []
        if isinstance(species, str):
            if species.lower() in ['all', '*']:
                species = self.available_species
            else:
                species = [species]
        for _species in species:
            try:
                hmm_paths = self.species_to_paths[_species]
            except KeyError:
                raise KeyError(f"{_species} is not a valid species from {self.available_species}")
            for hmm_path in hmm_paths:
                with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                    hmm = next(hmm_file)
                    hmms.append(hmm)
        return hmms
    
    def __transform_seq(self, seq_obj: Path) -> pyhmmer.easel.SequenceFile:
        """
        Transform a sequence into a Easel sequence file.

        Parameters
        ----------
        sequence: Union[Path, list, str]
            Sequence to be transformed.

        Returns
        -------
        pyhmmer.easel.SequenceFile
            Easel sequence file.
        """
        # TODO: add support for lists of SeqRecord, Paths, and stings with random given IDs
        with pyhmmer.easel.SequenceFile(seq_obj, digital=True) as seq_file:
            sequences = list(seq_file)
        return sequences

    def hmmersearch(self, sequences: str, species: Union[list, str]) -> List[str, int]:
        """
        Perform a HMMER search for a given sequence.

        Parameters
        ----------
        sequences : str
            Sequence to be searched.
        species : Union[list, str]
            Available species: ALL, alpaca, cat, cow, dog, human, mouse, pig, rabbit, rat, rhesus

        Returns
        -------
        List[str, int]
            List of HMMER hit attributes for ANARCI numbering.
        """
        sequences = self.__transform_seq(sequences)
        hmms = self.get_hmm_models(species)
        results = []
        for top_hits in pyhmmer.hmmsearch(hmms, sequences):
            for hit in top_hits:
                cog = hit.best_domain.alignment.hmm_name.decode()
                specie, chain_type = cog.split('_')
                results.append((hit.name.decode(), specie, chain_type, hit.score))
        return results