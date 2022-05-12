# from functools import lru_cache TODO: see if this is worth it for get_hmm_models
from operator import itemgetter
from pathlib import Path
from typing import Union, Optional, Dict, List, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyhmmer


class HMMER:
    """
    Extension of Pyhmmer with built-in alignment models that are species specific.
    """

    def __init__(self):
        # pathing
        self.hmm_folder = Path(__file__).parent / "data/anarci/HMMs"
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
            species, *chain = path.stem.split("_")
            try:
                species_to_paths[species].append(path)
            except KeyError:
                species_to_paths[species] = [path]

        return species_to_paths

    def get_hmm_models(
        self,
        species: Optional[Union[List[str], str]] = None,
    ) -> List[pyhmmer.plan7.HMMFile]:
        """
        Return a HMMER model for a given specie.

        Parameters
        ----------
        species: Optional[Union[List[str], str]]
            Available species: ALL, alpaca, cat, cow, dog, human, mouse, pig, rabbit, rat, rhesus

        Returns
        -------
        pyhmmer.plan7.HMMERModel
            HMM model for specific species
        """
        hmms = []

        if not species:
            species = self.available_species
        if isinstance(species, str):
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
        self, seq_objs: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str]
    ) -> List[pyhmmer.easel.DigitalSequence]:
        """
        Transform sequences or Fasta files into a list of Easel Digital objects for hmmer.

        Parameters
        ----------
        sequence: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str]
            Sequence to be transformed.

        Example
        -------
        >>> file = Path("test.fasta")  # only 1 sequence in the file
        >>> seq = Seq('DIVMTQSPLSLPVTPGEPASISCRSSQSLLYS')
        >>> seqrecord = SeqRecord('DIVMTQSPLSLPVTPGEPASISCRSSQSLLYS', id="unique name", description="long description")
        >>> loose_str = 'DIVMTQSPLSLPVTPGEPASISCRSSQSLLYS'

        >>> __transform_seq(file)  # read fasta file with amino acids
        [<pyhmmer.easel.DigitalSequence at 0x1396614c0>]

        >>> __transform_seq(seq)  # read Biopython Seq object; auto-assigns id via a counter starting at 0
        [<pyhmmer.easel.DigitalSequence at 0x139661540>]

        >>> __transform_seq(seqrecord)  # read Biopython SeqRecord object
        [<pyhmmer.easel.DigitalSequence at 0x139de98c0>]

        >>> __transform_seq(loose_str)  # read string of amino acids; auto-assigns id via a counter starting at 0
        [<pyhmmer.easel.DigitalSequence at 0x139deb800>]

        >>> __transform_seq([file, seq, seqrecord, loose_str])  # read all of the above together as 1 list output
        [
            <pyhmmer.easel.DigitalSequence at 0x1396614c0>,
            <pyhmmer.easel.DigitalSequence at 0x139661540>,
            <pyhmmer.easel.DigitalSequence at 0x139de98c0>,
            <pyhmmer.easel.DigitalSequence at 0x139deb800>
        ]

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
            if isinstance(seq_obj, (Seq, str)):
                sequences.append(self.__digitize_seq(name=str(sudo_name), seq=seq_obj))
                continue
            if isinstance(seq_obj, SeqRecord):
                sequences.append(self.__digitize_seq(name=seq_obj.id, seq=seq_obj.seq))
                continue

            raise ValueError(f"seq_obj {seq_obj} is not a valid sequence or path")

        if not sequences:
            raise ValueError(f"No valid sequences were found in {seq_objs}")

        return sequences

    def hmmsearch(
        self,
        sequences: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str],
        species: Optional[List[Union[str, int]]] = None,
        bit_score_threshold: int = 80,
        limit: int = 1,
        for_anarci: bool = False,
        for_j_region: bool = False,
    ) -> List[List[Dict[str, Union[str, int]]]]:
        """
        Perform a HMMER search for a given sequence.

        Notes
        -----
        1. Seperate HMMs is not a true chimeric search; the hmm need to be created together for that.
        2. There is a HMM loading overheard so ALL.hmm is much faster than grabbing specific species.

        Parameters
        ----------
        sequences : Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str]
            Sequences and/or Fasta files to be searched.
        species : Union[list, str]
            Available species: alpaca, cat, cow, dog, human, mouse, pig, rabbit, rat, rhesus
        bit_score_threshold : int
            Domain bit score threshold, default is 80.
        limit : int
            Number of domain hits to be returned, default is 1.
        for_anarci : bool
            If True, return the ANARCI expected state_vector object.
        for_j_region : bool
            If True, return the J-region expected state_vector object.
            Anarci reruns hmmsearch with the J-region to correct possible mistakes.

        Examples
        >>> seq = 'DVQLVESGGDLAKPGGSLRLTCVASGLSVTSNSMSWVRQAPGKGLRWVSTIWSKGGTYYADSVKGRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYDADYLHWYFDFWGQGALVTVSF'
        >>> HMMER().hmmsearch(seq, species=['human', 'mouse']', bit_score_threshold=80, limit=1)
        [[{
            'query': '0',
            'hmm_seq': 'qvqLvesGalelvkpgeslklsCaasGftlsllsssyalsWvrqapgkgLewvglisssaesgsteYaeslklgrvtisrdtskntlylqlsslraeDtavYyCarklll....llllfdvWGqGtlvtvs',
            'hmm_start': 1,
            'hmm_end': 127,
            'id': 'human_H',
            'description': '',
            'evalue': 8.02e-53,
            'bitscore': 164.4,
            'bias': 1.6,
            'query_seq': 'DVQLVESGG-DLAKPGGSLRLTCVASGLSV----TSNSMSWVRQAPGKGLRWVSTIWSK---GGTYYADSVK-GRFTVSRDSAKNTLYLQMDSLATEDTATYYCASIYHYdadyLHWYFDFWGQGALVTVS',
            'query_start': 0,
            'query_end': 122,
            'species': 'human',
            'chain_type': 'H'
        }]]

        Returns
        -------
        List[List[Dict[str, Union[str, int]]]]
            List of HMMER hits for ANARCI numbering.
        """
        # Convert sequences to Easel sequences
        sequences = self.__transform_seq(sequences)

        # Load Models by species
        hmms = self.get_hmm_models(species)

        # Maintain order of sequences since pyhmmer is async
        results = {seq.name.decode(): [] for seq in sequences}

        for top_hits in pyhmmer.hmmsearch(hmms, sequences):
            for hit in top_hits:

                domain = hit.best_domain
                ali = hit.best_domain.alignment

                if domain.score < bit_score_threshold:
                    continue

                results[hit.name.decode()].append(
                    {
                        "query": hit.name.decode(),
                        "hmm_seq": ali.hmm_sequence,
                        "hmm_start": ali.hmm_from,  # hmm seq starts with 1 and not 0
                        "hmm_end": ali.hmm_to,
                        "id": ali.hmm_name.decode(),
                        "description": hit.description or "",
                        "evalue": float("{:.2e}".format(domain.c_evalue)),
                        # "bitscore": round(hit.score, 1),  # TODO: anarci doesnt use hit score, but domain score; maybe use this as an option later?
                        "bitscore": round(domain.score, 1),
                        "bias": round(hit.bias, 1),
                        "query_seq": ali.target_sequence,
                        "query_start": ali.target_from - 1,  # target seq starts on pythonic index
                        "query_end": ali.target_to,
                        "species": ali.hmm_name.decode().split("_")[0],
                        "chain_type": ali.hmm_name.decode().split("_")[1],
                    }
                )

        # Sort by bitscore and limit results
        best_results = []
        for query_id in results.keys():
            best_results.append(sorted(results[query_id], key=lambda x: x["bitscore"], reverse=True)[:limit])

        # TODO: This will be deprecated in the future and will be replaced direct formatting
        if for_anarci:
            default_out = ([["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"]], [], [])
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
            # we want to keep every empty results as well for Anarci
            best_results = [
                result[0] if result else None for result in best_results
            ]  # Anarci only expects best result per query
            if not best_results:
                return [default_out]
            return [
                (
                    [anarci_keys, itemgetter(*anarci_keys)(result)],
                    [self.get_vector_state(**result)],
                    [result],
                )
                # if result['query_start'] < 4 and for_j_region is False else default_out
                if result else default_out
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
        """
        Get the ANARCI state_vector object.

        Parameters
        ----------
        hmm_seq : str
            HMM sequence.
        query_seq : str
            Query sequence.
        hmm_start : int
            HMM start alignmnet index; starts at 1
        hmm_end : int
            HMM end position.
        query_start : int
            Query start position; starts at 0 for pythonic index
        query_end : int
            Query end position.

        Notes
        -----
        'i' = insertion in HMM seq (aka a '.')
        'd' = deletion in Query seq (aka a '-')
        'm' = match or accepted missmatch where both are showing an amino acid in their seqs

        Example
        -------
        >>> hmm_seq = 'divltqsPsslsvsvgdrvtisCrasqsilesddgssylaWyqqkpgkapklliyaalllllllsslasGvPlsrfsGsGllsGtdftltissleaedvavyyCqqaklllllllltfGqGtkveik'
        >>> query_seq = 'DIVMTQSPLSLPVTPGEPASISCRSSQSLLYS-IGYNYLDWYLQKSGQSPQLLIYLG-------SNRASGVP-DRFSGSG--SGTDFTLKISRVEAEDVGFYYCMQAL----QTPYTFGQGTKLEIK'
        >>> hmm_start = 1
        >>> hmm_end = 127
        >>> query_start = 0
        >>> query_end = 122
        >>> get_vector_state(hmm_seq, query_seq, hmm_start, hmm_end, query_start, query_end)
        [
            ((1, 'm'), 0),
            ((2, 'm'), 1),
            ((3, 'm'), 2),
            ((4, 'm'), 3),
            ((5, 'm'), 4),
            ((6, 'm'), 5),
            ((7, 'm'), 6),
            ((8, 'm'), 7),
            ((9, 'm'), 8),
            ((10, 'm'), 9),
            ((11, 'm'), 10),
            ((12, 'm'), 11),
            ((13, 'm'), 12),
            ((14, 'm'), 13),
            ((15, 'm'), 14),
            ((16, 'm'), 15),
            ((17, 'm'), 16),
            ((18, 'm'), 17),
            ((19, 'm'), 18),
            ((20, 'm'), 19),
            ((21, 'm'), 20),
            ((22, 'm'), 21),
            ((23, 'm'), 22),
            ((24, 'm'), 23),
            ((25, 'm'), 24),
            ((26, 'm'), 25),
            ((27, 'm'), 26),
            ((28, 'm'), 27),
            ((29, 'm'), 28),
            ((30, 'm'), 29),
            ((31, 'm'), 30),
            ((32, 'm'), 31),
            ((33, 'd'), None),
            ((34, 'm'), 32),
            ((35, 'm'), 33),
            ((36, 'm'), 34),
            ((37, 'm'), 35),
            ((38, 'm'), 36),
            ((39, 'm'), 37),
            ((40, 'm'), 38),
            ((41, 'm'), 39),
            ((42, 'm'), 40),
            ((43, 'm'), 41),
            ((44, 'm'), 42),
            ((45, 'm'), 43),
            ((46, 'm'), 44),
            ((47, 'm'), 45),
            ((48, 'm'), 46),
            ((49, 'm'), 47),
            ((50, 'm'), 48),
            ((51, 'm'), 49),
            ((52, 'm'), 50),
            ((53, 'm'), 51),
            ((54, 'm'), 52),
            ((55, 'm'), 53),
            ((56, 'm'), 54),
            ((57, 'm'), 55),
            ((58, 'd'), None),
            ((59, 'd'), None),
            ((60, 'd'), None),
            ((61, 'd'), None),
            ((62, 'd'), None),
            ((63, 'd'), None),
            ((64, 'd'), None),
            ((65, 'm'), 56),
            ((66, 'm'), 57),
            ((67, 'm'), 58),
            ((68, 'm'), 59),
            ((69, 'm'), 60),
            ((70, 'm'), 61),
            ((71, 'm'), 62),
            ((72, 'm'), 63),
            ((73, 'd'), None),
            ((74, 'm'), 64),
            ((75, 'm'), 65),
            ((76, 'm'), 66),
            ((77, 'm'), 67),
            ((78, 'm'), 68),
            ((79, 'm'), 69),
            ((80, 'm'), 70),
            ((81, 'd'), None),
            ((82, 'd'), None),
            ((83, 'm'), 71),
            ((84, 'm'), 72),
            ((85, 'm'), 73),
            ((86, 'm'), 74),
            ((87, 'm'), 75),
            ((88, 'm'), 76),
            ((89, 'm'), 77),
            ((90, 'm'), 78),
            ((91, 'm'), 79),
            ((92, 'm'), 80),
            ((93, 'm'), 81),
            ((94, 'm'), 82),
            ((95, 'm'), 83),
            ((96, 'm'), 84),
            ((97, 'm'), 85),
            ((98, 'm'), 86),
            ((99, 'm'), 87),
            ((100, 'm'), 88),
            ((101, 'm'), 89),
            ((102, 'm'), 90),
            ((103, 'm'), 91),
            ((104, 'm'), 92),
            ((105, 'm'), 93),
            ((106, 'm'), 94),
            ((107, 'm'), 95),
            ((108, 'm'), 96),
            ((109, 'd'), None),
            ((110, 'd'), None),
            ((111, 'd'), None),
            ((112, 'd'), None),
            ((113, 'm'), 97),
            ((114, 'm'), 98),
            ((115, 'm'), 99),
            ((116, 'm'), 100),
            ((117, 'm'), 101),
            ((118, 'm'), 102),
            ((119, 'm'), 103),
            ((120, 'm'), 104),
            ((121, 'm'), 105),
            ((122, 'm'), 106),
            ((123, 'm'), 107),
            ((124, 'm'), 108),
            ((125, 'm'), 109),
            ((126, 'm'), 110),
            ((127, 'm'), 111)
        ]

        Returns
        -------
        List[Tuple[Tuple[int, str], int]]
            List of ANARCI state_vector objects.
        """
        assert len(hmm_seq) == len(query_seq), "The 2 seqs should be alignments of eachother"

        vector_state = []

        hmm_step = hmm_start  # real world index starting at 1 to match HMMER output
        query_step = query_start  # pythonic index starting at 0

        for i in range(len(hmm_seq)):
            # HMM seq insertion
            if hmm_seq[i] == ".":
                vector_state.append(((hmm_step, "i"), query_step))
                query_step += 1
            # Query seq deletion
            elif query_seq[i] == "-":
                vector_state.append(((hmm_step, "d"), None))
                hmm_step += 1
            # match or accepted missmatch
            else:
                vector_state.append(((hmm_step, "m"), query_step))
                hmm_step += 1
                query_step += 1

        return vector_state