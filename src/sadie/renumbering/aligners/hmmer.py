# from functools import lru_cache TODO: see if this is worth it for get_hmm_models
from __future__ import annotations

from operator import itemgetter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pyhmmer


def _ensure_str(value: str | bytes) -> str:
    """Convert bytes to str, handling both pyhmmer <0.12 (bytes) and >=0.12 (str)."""
    return value if isinstance(value, str) else value.decode()


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sadie.renumbering.clients import G3
from sadie.renumbering.numbering_translator import NumberingTranslator
from sadie.typing import Chain, Source, Species

# Species aliases for HMM lookup: maps backward-compat names to canonical HMM filenames
_HMM_SPECIES_ALIASES: Dict[str, str] = {
    "rhesus": "macaque",
}


def _use_local_hmm_builder() -> bool:
    """Check if germlines module should be used for HMM building."""
    try:
        from sadie.germlines.renumbering_integration import use_local_hmm_builder

        return use_local_hmm_builder()
    except ImportError:
        return False


class HMMER:
    """
    Extension of Pyhmmer to accept local HMMs (from Numbering) and external HMMs (from G3 or germlines).
    """

    g3 = G3()
    numbering = (
        NumberingTranslator()
    )  # TODO: merge this with G3 created HMMs and record which ones are legacy and are not built live via G3
    _local_hmm_builder = None  # Lazy initialization

    def __init__(
        self,
        species: Optional[Union[List[Species], Species]] = None,
        chains: Optional[Union[List[Chain], Chain]] = None,
        source: Source = "imgt",
        use_numbering_hmms: bool = False,
        hmm_dir: Optional[Path] = None,
        validate_pairs: bool = True,
    ) -> None:
        """Initialize HMMER with optional custom HMM directory.

        Parameters
        ----------
        species : Optional[Union[List[Species], Species]], optional
            Species to load HMMs for, by default None (loads all)
        chains : Optional[Union[List[Chain], Chain]], optional
            Chains to load HMMs for, by default None (loads all)
        source : Source, optional
            Source database (e.g., "imgt"), by default "imgt"
        use_numbering_hmms : bool, optional
            If True, force use of legacy Numbering HMMs, by default False
        hmm_dir : Optional[Path], optional
            Path to custom HMM directory. When provided, HMMs are loaded from
            `{hmm_dir}/{species}_{chain}.hmm`. Falls back to default HMM sources
            if custom HMM not found. By default None.
        validate_pairs : bool, optional
            If True (default), validate that every explicitly requested
            species/chain pair has a loaded HMM. Set to False when the caller
            has its own validation (e.g., Renumbering). By default True.
        """
        # Store custom HMM directory path
        self._hmm_dir = Path(hmm_dir) if hmm_dir else None

        # Force Numbering local HMMs to be used -- mostely for primiary testing
        self.hmms = self.get_hmm_models(
            species=species,
            chains=chains,
            source=source,
            use_numbering_hmms=use_numbering_hmms,
            validate_pairs=validate_pairs,
        )
        # place holders for hmmer
        self.alphabet = pyhmmer.easel.Alphabet.amino()

    def get_hmm_models(
        self,
        species: Optional[Union[List[Species], Species]] = None,
        chains: Optional[Union[List[Chain], Chain]] = None,
        source: Source = "imgt",
        use_numbering_hmms: bool = False,
        validate_pairs: bool = True,
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
        # Normalize species to list
        if species is None:
            species_list: List[str] = list(Species.species.values())
        elif isinstance(species, str):
            species_list = [str(species)]
        else:
            species_list = [str(single_species) for single_species in species]
        # Normalize chains to list
        if chains is None:
            chain_list: List[str] = list(Chain.chains)
        elif isinstance(chains, str):
            chain_list = [str(chains)]
        else:
            chain_list = [str(single_chain) for single_chain in chains]

        # Check if we should use local HMM builder from germlines
        use_local = _use_local_hmm_builder() and not use_numbering_hmms

        # Track which (species, chain) pairs successfully loaded an HMM.
        # Used by callers (e.g. Renumbering) to validate requested pairs.
        self._loaded_pairs: set[tuple[str, str]] = set()

        for single_species in species_list:
            # Resolve species aliases (e.g., "rhesus" → "macaque") so that
            # LocalHMMBuilder and custom-dir lookups find the correct HMM files.
            canonical_species = _HMM_SPECIES_ALIASES.get(single_species, single_species)
            for chain in chain_list:
                # Priority 0: Custom HMM directory (highest priority)
                if self._hmm_dir:
                    custom_hmm_path = self._hmm_dir / f"{canonical_species}_{chain}.hmm"
                    if custom_hmm_path.exists():
                        with pyhmmer.plan7.HMMFile(custom_hmm_path) as hmm_file:
                            hmms.append(next(hmm_file))
                        self._loaded_pairs.add((canonical_species, chain))
                        continue

                # Priority 1: Numbering (legacy) HMMs — specifically designed
                # for antibody numbering and produce correct IMGT gap patterns
                # (e.g. conserved deletion at position 10 for VH).
                numbering_species = single_species
                if numbering_species.strip() == "macaque":
                    numbering_species = "rhesus"
                if (numbering_species, chain) in self.numbering.species_chain_to_paths:
                    hmm_paths = self.numbering.species_chain_to_paths[(numbering_species, chain)]
                    for hmm_path in hmm_paths:
                        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                            hmm = next(hmm_file)
                            hmms.append(hmm)
                    self._loaded_pairs.add((canonical_species, chain))
                    continue

                # Priority 2: Local germlines HMM builder — fallback for
                # species/chains not covered by legacy Numbering HMMs.
                if use_local:
                    try:
                        if self._local_hmm_builder is None:
                            from sadie.germlines.renumbering_integration import (
                                LocalHMMBuilder,
                            )

                            self._local_hmm_builder = LocalHMMBuilder()

                        hmm = self._local_hmm_builder.get_hmm(species=canonical_species, chain=chain, source=source)
                        hmms.append(hmm)
                        self._loaded_pairs.add((canonical_species, chain))
                        continue
                    except Exception as e:
                        # Fall through to G3 on error
                        import logging

                        logging.warning(
                            f"Local HMM builder failed for {canonical_species} {chain}: {e}. " f"Falling back to G3."
                        )

                # Priority 3: Build G3 HMMs (final fallback)
                # g3.jordanrwillis.com may be unreachable; degrade gracefully.
                try:
                    if chain not in self.g3.chains or single_species not in self.g3.species or use_numbering_hmms is True:
                        continue  # No HMM available for this pair
                    hmm = self.g3.get_hmm(
                        source=source,
                        chain=chain,
                        species=single_species,
                        limit=None,
                    )
                    hmms.append(hmm)
                    self._loaded_pairs.add((canonical_species, chain))
                except Exception as e:
                    import logging
                    logging.warning(f"G3 fallback unavailable for {single_species} {chain}: {e}")

        # Guard rail: when both species and chains were explicitly provided
        # (not defaults), validate every requested pair has a loaded HMM.
        # This catches mixed requests where some pairs load successfully but
        # others are silently dropped (e.g., alpaca H+K where only H exists).
        if validate_pairs and species is not None and chains is not None:
            errors = []
            for single_species in species_list:
                canonical = _HMM_SPECIES_ALIASES.get(single_species, single_species)
                for chain in chain_list:
                    if (canonical, chain) not in self._loaded_pairs:
                        supported = self._get_supported_chains(canonical, source, use_local, use_numbering_hmms)
                        errors.append(
                            f"No HMM available for {canonical} {chain}. "
                            f"Supported chains for {canonical}: {supported}"
                        )
            if errors:
                raise ValueError("; ".join(errors))

        return hmms

    def _get_supported_chains(self, species: str, source: str, use_local: bool, use_numbering_hmms: bool) -> List[str]:
        """Determine which chains have HMMs available for a given species.

        Parameters
        ----------
        species : str
            Canonical species name (e.g., "alpaca", "chicken")
        source : str
            Data source (e.g., "imgt")
        use_local : bool
            Whether to check LocalHMMBuilder
        use_numbering_hmms : bool
            Whether legacy Numbering HMMs are forced

        Returns
        -------
        List[str]
            List of chain types with available HMMs (e.g., ["H"] or ["H", "K", "L"])
        """
        all_chains = ["H", "K", "L", "A", "B", "G", "D"]
        supported = []
        for chain in all_chains:
            # Check custom HMM directory
            if self._hmm_dir:
                if (self._hmm_dir / f"{species}_{chain}.hmm").exists():
                    supported.append(chain)
                    continue
            # Check local HMM builder
            if use_local and not use_numbering_hmms:
                try:
                    if self._local_hmm_builder is not None:
                        self._local_hmm_builder.get_hmm(species=species, chain=chain, source=source)
                        supported.append(chain)
                        continue
                except Exception:
                    pass
            # Check legacy Numbering paths
            if (species, chain) in self.numbering.species_chain_to_paths:
                supported.append(chain)
        return supported

    def digitize_seq(self, name: str | bytes, seq: Seq | str) -> pyhmmer.easel.DigitalSequence:
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
        try:
            return pyhmmer.easel.TextSequence(name=name, sequence=seq).digitize(self.alphabet)
        except TypeError:
            # pyhmmer >=0.12 expects str name, not bytes
            if isinstance(name, bytes):
                name = name.decode()
            return pyhmmer.easel.TextSequence(name=name, sequence=seq).digitize(self.alphabet)

    def transform_seqs(
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

        >>> transform_seqs(file)  # read fasta file with amino acids
        [<pyhmmer.easel.DigitalSequence at 0x1396614c0>]

        >>> transform_seqs(seq)  # read Biopython Seq object; auto-assigns id via a counter starting at 0
        [<pyhmmer.easel.DigitalSequence at 0x139661540>]

        >>> transform_seqs(seqrecord)  # read Biopython SeqRecord object
        [<pyhmmer.easel.DigitalSequence at 0x139de98c0>]

        >>> transform_seqs(loose_str)  # read string of amino acids; auto-assigns id via a counter starting at 0
        [<pyhmmer.easel.DigitalSequence at 0x139deb800>]

        >>> transform_seqs([file, seq, seqrecord, loose_str])  # read all of the above together as 1 list output
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
            if isinstance(seq_obj, str):
                if len(seq_obj) < 4096:  # max length of a path
                    if Path(seq_obj).is_file():
                        with pyhmmer.easel.SequenceFile(seq_obj, digital=True) as seq_file:
                            sequences.extend(list(seq_file))
                            continue
            # If sequence is a string, digitize it directly and add it to the list
            if isinstance(seq_obj, (Seq, str)):
                sequences.append(self.digitize_seq(name=str(sudo_name), seq=seq_obj))
                continue
            # If the sequence is a path, open it and digitize
            if isinstance(seq_obj, Path):
                with pyhmmer.easel.SequenceFile(seq_obj, digital=True) as seq_file:
                    sequences.extend(list(seq_file))
                    continue
            if isinstance(seq_obj, SeqRecord):
                sequences.append(self.digitize_seq(name=seq_obj.id, seq=seq_obj.seq))
                continue
            if isinstance(seq_obj, tuple):
                seq_id, seq = seq_obj
                sequences.append(self.digitize_seq(name=seq_id, seq=seq))
                continue
            raise ValueError(f"seq_obj {seq_obj} is not a valid sequence or path")

        return sequences

    # def load_stockholm(self, stockholm_path: Union[Path, str]) -> pyhmmer.easel.DigitalMSA:
    #     with pyhmmer.easel.MSAFile(stockholm_path, digital=True, alphabet=self.alphabet) as msa_file:
    #         msa = next(msa_file)
    #     return msa

    def hmmsearch(
        self,
        sequences: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str],
        bit_score_threshold: int = 80,
        limit: Optional[int] = 1,
        for_numbering: bool = False,
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
        bit_score_threshold : int
            Domain bit score threshold, default is 80.
        limit : int
            Number of domain hits to be returned, default is 1.
        for_numbering : bool
            If True, return the NUMBERING expected state_vector object.
        # for_j_region : bool
        #     If True, return the J-region expected state_vector object.
        #     Numbering reruns hmmsearch with the J-region to correct possible mistakes.

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
            List of HMMER hits for NUMBERING numbering.
        """
        # Convert sequences to Easel sequences
        sequences = self.transform_seqs(sequences)  # type: ignore[assignment]

        # Multiprocessing might scramble actual seq from order
        seq_name_2_seq = {_ensure_str(seq.name): seq.textize().sequence for seq in sequences}

        # Maintain order of sequences since pyhmmer is async
        results = {_ensure_str(seq.name): [] for seq in sequences}

        # HMMs not big enough to be worth multiprocessing
        for top_hits in pyhmmer.hmmsearch(self.hmms, sequences, cpus=1):
            for _i, hit in enumerate(top_hits):
                domain = hit.best_domain
                ali = hit.best_domain.alignment

                if domain.score < bit_score_threshold:
                    continue

                hit_name = _ensure_str(hit.name)
                hmm_name = _ensure_str(ali.hmm_name)
                results[hit_name].append(
                    {
                        "order": 0,  # best domain always has order 0
                        # "order": _i,  # no need for it unless we allow multiple domain hits
                        "n": len(hit.domains),
                        "query": hit_name,
                        "query_length": len(seq_name_2_seq[hit_name]),
                        "hmm_seq": ali.hmm_sequence,
                        "hmm_start": ali.hmm_from - 1,  # hmm seq starts with 1 and not 0
                        "hmm_end": ali.hmm_to,
                        "id": hmm_name,
                        "description": hit.description or "",
                        "evalue": float("{:.2e}".format(domain.c_evalue)),
                        # "bitscore": round(hit.score, 1),  # TODO: numbering doesnt use hit score, but domain score; maybe use this as an option later?
                        "bitscore": round(domain.score, 1),
                        "bias": round(hit.bias, 1),
                        "query_seq": ali.target_sequence.upper(),
                        "query_start": ali.target_from - 1,  # target seq starts on pythonic index
                        "query_end": ali.target_to,
                        "species": hmm_name.split("_")[0],
                        "chain_type": hmm_name.split("_")[1],
                    }
                )

        # Sort by bitscore and limit results
        best_results = []
        for query_id in results.keys():
            best_results.append(sorted(results[query_id], key=lambda x: x["bitscore"], reverse=True)[:limit])

        # TODO: This will be deprecated in the future and will be replaced direct formatting
        if for_numbering:
            default_out = ([["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"]], [], [])
            numbering_keys = [
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
            # we want to keep every empty results as well for Numbering
            # TODO: this is per sequence so in the future we would want to expand this to no limit once the check_for_j is rooted out.
            best_results = [
                result[0] if result else None for result in best_results
            ]  # Numbering only expects best result per query
            return [
                (
                    (
                        [numbering_keys, itemgetter(*numbering_keys)(result)],
                        [self.get_vector_state(**result)],
                        [result],
                    )
                    if result
                    else default_out
                )
                for result in best_results
            ]

        return best_results

    def get_vector_state(
        self,
        query: str,
        order: int,
        n: int,
        hmm_seq: str,
        query_seq: str,
        query_length: int,
        hmm_start: int,
        hmm_end: Optional[int] = None,
        query_start: Optional[int] = None,
        query_end: Optional[int] = None,
        **kwargs,
    ) -> List[Tuple[Tuple[int, str], int]]:
        """
        Get the NUMBERING state_vector object.

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
            List of NUMBERING state_vector objects.
        """
        assert len(hmm_seq) == len(query_seq), "The 2 seqs should be alignments of eachother"
        hmm_length = 128  # hardcoded since this is the length of the HMM for an antibody

        #
        # Allowing the user simple numbering if they already have alignments
        # if hmm_start is None:
        #     hmm_start = 0
        # if hmm_end is None:
        #     hmm_end = len(hmm_seq) - hmm_seq.count(".")
        # if query_start is None:
        #     query_start = 0
        # if query_end is None:
        #     query_end = len(query_seq) - query_seq.count("-")

        # print('0', hmm_length, hmm_start, hmm_end, query_length, query_start, query_end)

        if (order == 0) and (0 < hmm_start < 5):
            # print('prefix', hmm_length, hmm_start, hmm_end, query_length, query_start, query_end)
            n_extend = hmm_start
            if hmm_start > query_start:
                n_extend = min(query_start, hmm_start - query_start)
            query_seq = "8" * n_extend + query_seq
            hmm_seq = "x" * n_extend + hmm_seq
            query_start -= n_extend
            hmm_start -= n_extend
        # print('1', hmm_length, hmm_start, hmm_end, query_length, query_start, query_end)

        if n == 1 and query_end < query_length and (123 < hmm_end < hmm_length):  # Extend forwards
            # print('suffix', hmm_length, hmm_start, hmm_end, query_length, query_start, query_end)
            n_extend = min(hmm_length - hmm_end, query_length - query_end)
            query_seq += "8" * n_extend
            hmm_seq += "x" * n_extend
            query_end += n_extend
            hmm_end += n_extend
        # print('2', hmm_length, hmm_start, hmm_end, query_length, query_start, query_end)

        vector_state = []

        all_reference_states = list(range(1, 129))
        hmm_step = hmm_start  # real world index starting at 1 to match HMMER output
        query_step = query_start  # pythonic index starting at 0

        for i in range(len(hmm_seq)):
            # HMM seq insertion
            if hmm_seq[i] == ".":
                vector_state.append(((all_reference_states[hmm_step], "i"), query_step))
                query_step += 1
            # Query seq deletion
            elif query_seq[i] == "-":
                vector_state.append(((all_reference_states[hmm_step], "d"), None))
                hmm_step += 1
            # match or accepted missmatch
            else:
                vector_state.append(((all_reference_states[hmm_step], "m"), query_step))
                hmm_step += 1
                query_step += 1

        # print(vector_state)

        return vector_state

    def check_for_j(
        self,
        sequences: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str],
        alignments: List[List[Dict[str, Union[str, int]]]],
    ):
        """
        As the length of CDR3 gets long (over 30ish) an alignment that does not include the J region becomes more favourable.
        This leads to really long CDR3s not being numberable.

        To overcome this problem, when no J region is detected we try without the v region.
        """
        for i in range(len(sequences)):
            # Check the alignment for J region
            if len(alignments[i][1]) == 1:  # Only do for single domain chains.
                # Check whether a J region has been identified. If not check whether there is still a considerable amount of sequence
                # remaining.
                ali = alignments[i][1][0]

                # Find the last match position.
                last_state = ali[-1][0][0]  # type: ignore[index]
                last_si = ali[-1][1]  # type: ignore[index]
                if last_state < 120:  # No or very little J region
                    if last_si + 30 < len(
                        sequences[i][1]  # type: ignore[index]
                    ):  # Considerable amount of sequence left...suspicious of a long CDR3
                        # Find the position of the conserved cysteine (imgt 104).
                        cys_si = dict(ali).get((104, "m"), None)
                        if cys_si is not None:  # 104 found.
                            # Find the corresponding index in the alignment.
                            cys_ai = ali.index(((104, "m"), cys_si))  # type: ignore[attr-defined]

                            # Try to identify a J region in the remaining sequence after the 104. A low bit score threshold is used.
                            _, re_states, re_details = self.hmmsearch(
                                sequences=[(sequences[i][0], sequences[i][1][cys_si + 1 :])],  # type: ignore[index]
                                bit_score_threshold=10,
                                for_numbering=True,
                            )[0]

                            # Check if a J region was detected in the remaining sequence.
                            if re_states and re_states[0][-1][0][0] >= 126 and re_states[0][0][0][0] <= 117:  # type: ignore[index]
                                # Sandwich the presumed CDR3 region between the V and J regions.

                                vRegion = ali[: cys_ai + 1]  # type: ignore[index]
                                jRegion = [
                                    (state, index + cys_si + 1) for state, index in re_states[0] if state[0] >= 117
                                ]
                                cdrRegion = []
                                next = 105
                                for si in range(cys_si + 1, jRegion[0][1]):
                                    if next >= 116:
                                        cdrRegion.append(((116, "i"), si))
                                    else:
                                        cdrRegion.append(((next, "m"), si))
                                        next += 1

                                # Update the alignment entry.
                                alignments[i][1][0] = vRegion + cdrRegion + jRegion
                                alignments[i][2][0]["query_end"] = jRegion[-1][1] + 1  # type: ignore[index]

    def hmmersearch_with_j(
        self,
        sequences: Union[List[Union[Path, SeqRecord, str]], Path, SeqRecord, str],
        bit_score_threshold: int = 80,
        limit: Optional[int] = 1,
        for_numbering: bool = False,
        threads: int = 1,
    ) -> List[Tuple[str, List[List[Tuple[Tuple[int, str], int]]], List[Dict[str, Union[str, int]]]]]:
        """
        Run a HMMER search on the given sequences using the J region HMMs.
        """
        alignments = self.hmmsearch(
            sequences=sequences, bit_score_threshold=bit_score_threshold, limit=limit, for_numbering=for_numbering
        )
        # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
        # Check for J regions. Modifies sequence and alignment in place.
        if for_numbering:
            self.check_for_j(sequences=sequences, alignments=alignments)

        return alignments
