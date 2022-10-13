from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import List, Optional, Tuple

import pyhmmer
import requests as r
from pydantic import validate_arguments
from yarl import URL

from sadie.typing import Chain, Source, Species


class G3:
    """API Wrapper with OpenAPI found here https://g3.jordanrwillis.com/docs"""

    # TODO: most likely make this an import
    data_folder = Path(__file__).parent.parent / "data"
    segments = {"V", "D", "J"}
    chains = {"H", "K", "L"}

    def __init__(self):
        self.base_url = URL("https://g3.jordanrwillis.com/api/v1")
        self.not_usable_species = [
            "pig",
            "cow",
            "cat",  # missing L
            "alpaca",  # missing L and K
            "rhesus",  # TODO: breaks tests; fix and fall back on numbering for now
            "dog",  # TODO: viable but does not match. Need to check if diff species of dog from G3
        ]
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.builder = pyhmmer.plan7.Builder(self.alphabet, architecture="hand")
        self.background = pyhmmer.plan7.Background(self.alphabet)

    @property
    @lru_cache(maxsize=1)
    def sources(self):
        resp = r.get(self.base_url)
        resp.raise_for_status()
        return resp.json()["components"]["schemas"]["SourceName"]["enum"]

    @property
    @lru_cache(maxsize=1)
    def species(self):
        resp = r.get(self.base_url)
        resp.raise_for_status()
        species = resp.json()["components"]["schemas"]["CommonName"]["enum"]
        return [single_species for single_species in species if single_species not in self.not_usable_species]

    @lru_cache(maxsize=None)
    @validate_arguments
    def __get_gene_resp(
        self,
        source: Source = "imgt",
        species: Species = "human",
        segment: str = "V",
        limit: Optional[int] = None,
    ) -> str:
        params = {
            "source": source,
            "common": species,
            "segment": segment,
            "limit": limit if limit else "-1",
        }
        resp = r.get(self.base_url / "genes", params=params)
        resp.raise_for_status()
        return resp

    @validate_arguments
    def get_gene(
        self,
        source: Source = "imgt",
        species: Species = "human",
        chain: Chain = "H",
        segment: str = "V",
        limit: Optional[int] = None,
    ) -> str:
        resp = self.__get_gene_resp(source=source, species=species, segment=segment, limit=limit)
        return [x for x in resp.json() if x["gene"][2].lower() == chain.lower()]

    def get_stockholm_pairs(
        self,
        source: Source = "imgt",
        chain: Chain = "H",
        species: Species = "human",
        limit: Optional[int] = None,
    ) -> List[Tuple[str, str]]:

        sub_v = self.get_gene(source=source, species=species, chain=chain, segment="V", limit=limit)
        sub_j = self.get_gene(source=source, species=species, chain=chain, segment="J", limit=limit)

        stockholm_pairs = []
        for merge in product(sub_v, sub_j):

            v_seg = merge[0]
            j_seg = merge[1]

            functional = v_seg["imgt"]["imgt_functional"]
            v_part = v_seg["imgt"]["sequence_gapped_aa"].replace(".", "-")[:108].ljust(108).replace(" ", "-")
            cdr3_part = j_seg["imgt"]["cdr3_aa"]
            fwr4_part = j_seg["imgt"]["fwr4_aa"]
            v_name = v_seg["gene"]
            j_name = j_seg["gene"]

            name = f"{species}_{v_name}_{j_name}"

            # why?
            if functional != "F":
                continue

            # H rules
            # if chain.strip().lower() in "h":
            #     if len(cdr3_part[-3:] + fwr4_part) == 13:
            #         fwr4_part += "-"

            # K rules
            if chain.strip().lower() in "k":
                if len(cdr3_part[-3:] + fwr4_part) in [12, 13]:
                    fwr4_part += "-"

            # # L rules
            if chain.strip().lower() == "l":
                if len(cdr3_part[-3:] + fwr4_part) == 12:
                    fwr4_part += "-"

            # todo: alt fwr4_part based on it's size and who's askin
            multiplier = 128 - (len(v_part) + len(cdr3_part[-3:] + fwr4_part))

            align = v_part + "-" * multiplier + cdr3_part[-3:] + fwr4_part

            # sanity check if chains rules are working
            assert len(align) == 128

            stockholm_pairs.append((name, align))

        return stockholm_pairs

    # def get_msa(
    #     self,
    #     source: Source = "imgt",
    #     species: Species = "human",
    #     chain: Chain = "H",
    #     limit: Optional[int] = None,
    # ) -> str:
    #     stockholm_pairs = self.get_stockholm_pairs(source=source, chain=chain, species=species, limit=limit)
    #     sequences = []
    #     for name, align in stockholm_pairs:
    #         sequence = pyhmmer.easel.TextSequence(name=name.encode(), sequence=align)
    #         sequences.append(sequence)
    #     if not sequences:
    #         return None
    #     return pyhmmer.easel.TextMSA(name=f"{species}_{chain}".encode(), sequences=sequences).digitize(self.alphabet)

    @lru_cache(maxsize=None)
    def build_stockholm(
        self,
        source: Source = "imgt",
        species: Species = "human",
        chain: Chain = "H",
        limit: Optional[int] = None,
    ) -> Path:
        """
        Get a stockholm file in string format for the given species and chain.

        Parameters
        ----------
        source : str, optional
            Source of gene data, by default "imgt"
            options: 'imgt' or 'custom'
        species : str, optional
            species selected from avaliabe, by default "human"
        chain : str, optional
            chain for seq, by default "H"
            options: 'H', 'k', 'l' -> heavy, kappa light
        Returns
        -------
        str
            stockholm file in string format
        """
        if (self.data_folder / "stockholms").exists() is False:
            (self.data_folder / "stockholms").mkdir(parents=True, exist_ok=True)

        sto_path = self.data_folder / f"stockholms/{species}_{chain}.sto"

        sto_pairs = self.get_stockholm_pairs(source=source, chain=chain, species=species, limit=limit)

        head = f"# STOCKHOLM 1.0\n#=GF ID {species}_{chain}\n"
        body = "\n".join([f"{name}\t{ali}" for name, ali in sto_pairs])
        tail = "\n#=GC RF" + "\t" + "x" * 128 + "\n//\n"

        # TODO: hand arch needs a parsed file -- will be refactored to handle digital directly
        with open(sto_path, "w") as outfile:
            outfile.write(head + body + tail)

        return sto_path

    # @lru_cache(maxsize=None)
    def build_hmm(
        self,
        source: Source = "imgt",
        species: Species = "human",
        chain: Chain = "H",
        limit: Optional[int] = None,
    ) -> Path:
        sto_path = self.build_stockholm(source=source, chain=chain, species=species, limit=limit)

        if (self.data_folder / "hmms").exists() is False:
            (self.data_folder / "hmms").mkdir(parents=True, exist_ok=True)

        hmm_path = self.data_folder / f"hmms/{species}_{chain}.hmm"

        with pyhmmer.easel.MSAFile(sto_path, digital=True, alphabet=self.alphabet, format="stockholm") as msa_file:
            msa = next(msa_file)
            hmm, _, _ = self.builder.build_msa(msa, self.background)

        with open(hmm_path, "wb") as output_file:
            hmm.write(output_file)  # type: ignore

        return hmm_path

    # @lru_cache(maxsize=None)
    def get_hmm(
        self,
        source: Source = "imgt",
        species: Species = "human",
        chain: Chain = "H",
        limit: Optional[int] = None,
    ):
        hmm_path = self.data_folder / f"hmms/{species}_{chain}.hmm"

        if hmm_path.is_file() is False:
            hmm_path = self.build_hmm(source=source, chain=chain, species=species, limit=limit)

        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            hmm = next(hmm_file)

        return hmm  # type: ignore
