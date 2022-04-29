# from functools import lru_cache TODO: see if this is worth it
from collections import defaultdict
from operator import itemgetter
from pathlib import Path
from typing import Union, List, Tuple

# from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyhmmer

# Sequence = Union[List[Path, SeqRecord, str], Path, SeqRecord, str]


class HMMER:
    """
    Extension of Pyhmmer with built-in alignment models that are specie specific.
    """

    def __init__(self):
        # pathing 
        self.hmm_folder = Path(__file__).parent.parent / "data/anarci/HMMs"
        self.hmm_paths = self.hmm_folder.glob("*.hmm")
        # species specific
        self.species_to_paths = self.get_species_to_paths()
        self.available_species = sorted(self.species_to_paths.keys())
        # place holders for hmmer
        self.alphabet = pyhmmer.easel.Alphabet.amino()

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
    
    def __digitize_seq(self, name: str, seq: str) -> pyhmmer.easel.DigitalSequence:
        """
        Digitize a sequence for hmmer. 

        Parameters
        ----------
        name : str
            id or name of the sequence
        seq : str
            amino acid sequence to be digitized

        Returns
        -------
        pyhmmer.easel.DigitalSequence
            Digitized sequence.
        """
        if isinstance(name, str):
            name = name.encode()
        if isinstance(seq, Seq):
            seq = str(seq)
        return pyhmmer.easel.TextSequence(name=name, sequence=seq).digitize(self.alphabet)
        
    def __transform_seq(
        self, 
        seq_objs: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str]
    ) -> List[pyhmmer.easel.DigitalSequence]:
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
        sequences = []
        
        if not isinstance(seq_objs, (list, set, tuple)):
            seq_objs = [seq_objs]
        
        for sudo_name, seq_obj in enumerate(seq_objs):
            
            # If the sequence is a path, open it and digitize
            if isinstance(seq_obj, Path):
                with pyhmmer.easel.SequenceFile(seq_obj, digital=True) as seq_file:
                    sequences.extend(list(seq_file))
                    continue
            if isinstance(seq_obj, str):
                if len(seq_obj) < 4096:
                    if Path(seq_obj).is_file():
                        with pyhmmer.easel.SequenceFile(seq_obj, digital=True) as seq_file:
                            sequences.extend(list(seq_file))
                            continue
                        
            # If sequence is a string, digitize it directly and add it to the list
            if isinstance(seq_obj, str):
                sequences.append(self.__digitize_seq(name=str(sudo_name), seq=seq_obj))
                continue
            if isinstance(seq_obj, SeqRecord):
                sequences.append(self.__digitize_seq(name=seq_obj.id, seq=seq_obj.seq))
                continue
            
            raise ValueError(f'seq_obj {seq_obj} is not a valid sequence or path')      
                
        return sequences

    def hmmsearch(
        self,
        sequences: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str],
        species: List[Union[str, int]] = [
            "alpaca",
            "cat",
            "cow",
            "dog",
            "human",
            "mouse",
            "pig",
            "rabbit",
            "rat",
            "rhesus",
        ],
        bit_score_threshold: int = 80,  # TODO: doesn't work for now
        limit: int = 1,
        for_anarci: bool = False,
    ) -> List[Union[str, int]]:
        """
        Perform a HMMER search for a given sequence.

        Notes
        -----
        1. Seperate HMMs is not a true chimeric search; the hmm need to be created together for that.
        2. There is a HMM loading overheard so ALL.hmm is much faster than grabbing specific species.

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
        if not sequences:
            return None
        hmms = self.get_hmm_models(species)
        results = defaultdict(list)
        for top_hits in pyhmmer.hmmsearch(hmms, sequences):
            for sequence_index, hit in enumerate(top_hits):
                domain = hit.best_domain
                ali = hit.best_domain.alignment
                results[hit.name.decode()].append(
                    {
                        "query": hit.name.decode(),
                        "hmm_seq": ali.hmm_sequence,
                        "hmm_start": ali.hmm_from,  # hmm seq starts with 1 and not 0
                        "hmm_end": ali.hmm_to,
                        # "seq": sequences[sequence_index].textize().sequence,
                        "id": ali.hmm_name.decode(),
                        "description": hit.description or "",
                        "evalue": float("{:.2e}".format(domain.c_evalue)),
                        # "bitscore": round(hit.score, 1),
                        "bitscore": round(domain.score, 1),
                        "bias": round(hit.bias, 1),
                        "query_seq": ali.target_sequence,
                        "query_start": ali.target_from - 1,  # target seq starts on pythonic index
                        "query_end": ali.target_to,
                        "species": ali.hmm_name.decode().split("_")[0],
                        "chain_type": ali.hmm_name.decode().split("_")[1],
                    }
                )
        best_results = []
        for query_id, result in results.items():
            best_results.append(sorted(results[query_id], key=lambda x: x["bitscore"], reverse=True)[:limit])

        # TODO: This will be deprecated in the future and will be replaced direct formatting
        if for_anarci:
            anarci_keys = [
                "id",
                "description",
                "evalue",
                "bitscore",
                "bias",
                "query_start",
                "query_end",
                "species",
                "chain_type",
            ]
            best_results = [result[0] for result in best_results]  # Anarci only expects best result per query
            return [
                (
                    [anarci_keys, itemgetter(*anarci_keys)(result)],
                    [self.get_vector_state(**result)],
                    [result],
                )
                for result in best_results
            ]

        return best_results

    def get_vector_state(
        self, 
        hmm_seq: str, 
        query_seq: str, 
        hmm_start: int, 
        hmm_end: int, 
        query_start: int, 
        query_end: int,
        **kwargs,
    ) -> List[Tuple[Tuple[int, str], int]]:
        # assert target_start < 5, "Query start should not be more than 5"
        assert len(hmm_seq) == len(query_seq), "The 2 seqs should be alignments of eachother"

        vector_state = []

        hmm_step = hmm_start
        query_step = query_start

        for i in range(len(hmm_seq)):
            if hmm_seq[i] == ".":
                vector_state.append(((hmm_step, "i"), query_step))
                query_step += 1
            elif query_seq[i] == "-":
                vector_state.append(((hmm_step, "d"), None))
                hmm_step += 1
            else:
                vector_state.append(((hmm_step, "m"), query_step))
                hmm_step += 1
                query_step += 1

        return vector_state
