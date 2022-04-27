# from functools import lru_cache TODO: see if this is worth it
from collections import namedtuple
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
        self.hmm_paths = self.hmm_folder.glob("*.hmm")
        self.species_to_paths = self.get_species_to_paths()
        self.available_species = sorted(self.species_to_paths.keys())
        self.Result = namedtuple("Result", ['query', 'id', 'description', 'evalue', 'bitscore', 'bias', 'query_start', 'query_end', 'species', 'chain_type'])

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
            species = path.stem.split("_")[0]
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
        # if isinstance(species, str):
        #     if species.lower() in ["all", "*"]:
        #         species = self.available_species
        #     else:
        #         species = [species]
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

    def hmmsearch(
        self, sequences: str, 
        species: List[Union[str, int]] = 'ALL', 
        bit_score_threshold: int = 80,  # TODO: doesn't work for now
        limit: int = 5,  # todo: doesnt work for now
        vector_state_format: bool = False,  # TODO: refactor legacy
    ) -> List[Union[str, int]]:
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
            for i, hit in enumerate(top_hits):
                results.append(self.Result(
                    hit.name.decode(),
                    hit.best_domain.alignment.hmm_name.decode(),
                    hit.description,
                    hit.best_domain.c_evalue,
                    round(hit.best_domain.score, 1),
                    round(hit.best_domain.bias, 1),
                    hit.best_domain.env_from - 1,
                    hit.best_domain.env_to - 1,
                    hit.best_domain.alignment.hmm_name.decode().split('_')[0],
                    hit.best_domain.alignment.hmm_name.decode().split('_')[1],
                ))
        # TODO: this is a bit of a hack, but it works for now
        best_results = {}
        keep_query = set()
        for result in results:
            if result.query in best_results:
                previous_bitscore = best_results[result.query].bitscore
                best_results[result.query].bitscore
                if result.bitscore > previous_bitscore:
                    best_results[result.query] = result
                    keep_query.add(result.query)
                elif result.bitscore == previous_bitscore:
                    if best_results[result.query].id != result.id:
                        keep_query.remove(result.query)
            else:
                best_results[result.query] = result
                keep_query.add(result.query)
        filtered_results = [best_results[k] for k in sorted(best_results) if k in keep_query]
            
        return filtered_results

    def get_vector_state(self, alignment: pyhmmer.plan7.Alignment, seq: str):
        # dont let target_from be larger than 5
        # pad 'm' to begining and end of target so we can format it properly
        # padding can be assummed acceptable since threshold is set
        vector_state = []
        
        hmm_aln = alignment.hmm_sequence
        seq_aln = alignment.target_sequence
        hmm_step = alignment.hmm_from  # starts from 1; the real world index; idk anarci did this
        seq_step = alignment.target_from - 1  # starts at 0 for pythonic index
        
        # pad front
        if seq_step > 0:
            pad_hmm_step = hmm_step - seq_step
            for i in range(seq_step):
                vector_state.append(((pad_hmm_step + i, 'm'), i))
        
        for i in range(len(hmm_aln)):
            if hmm_aln[i] == '.':
                vector_state.append(((hmm_step, 'i'), seq_step))
                seq_step += 1
            elif seq_aln[i] == '-':
                vector_state.append(((hmm_step, 'd'), None))
                hmm_step += 1
            else:
                vector_state.append(((hmm_step, 'm'), seq_step))
                hmm_step += 1
                seq_step += 1
                
        # pad end
        hmm_step = vector_state[-1][0][0]
        seq_step = vector_state[-1][1]
        if seq_step < len(seq) - 1:
            for seq_step in range(seq_step, len(seq)):
                vector_state.append(((hmm_step, 'm'), seq_step))
                hmm_step += 1
                
        return vector_state
