from functools import lru_cache
from itertools import product
from typing import Optional

import pyhmmer
import requests as r
from yarl import URL


class G3:
    """API Wrapper with OpenAPI found here https://g3.jordanrwillis.com/docs"""
    def __init__(self):
        self.base_url = URL("https://g3.jordanrwillis.com/api/v1")
        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self.builder = pyhmmer.plan7.Builder(self.alphabet)
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
        return resp.json()["components"]["schemas"]["CommonName"]["enum"]

    @lru_cache(maxsize=None)
    def __get_gene_resp(
        self,
        source: str = "imgt",
        species: str = "human",
        segment: str = "V",
        limit: Optional[int] = None,
    ) -> str:
        params = {
            "source": source.strip().lower(),
            "common": species.strip().lower(),
            "segment": segment.strip().upper(),
            "limit": limit if limit else "-1",
        }
        resp = r.get(self.base_url / "genes", params=params)
        resp.raise_for_status()
        return resp
        
    def get_gene(
        self,
        source: str = "imgt",
        species: str = "human",
        chain: str = "h",
        segment: str = "V",
        limit: Optional[int] = None,
    ) -> str:
        resp = self.__get_gene_resp(
            source=source, species=species, segment=segment, limit=limit
        )
        return [
            x
            for x in resp.json()
            if x["gene"][2].lower() == chain.lower()
            and x["gene"].split("*")[-1] == "01"
        ]

    def get_msa(
        self,
        source: str = "imgt",
        chain: str = "h",
        species: str = "human",
        limit: Optional[int] = None,
    ) -> str:
        sub_v = self.get_gene(
            source=source, species=species, chain=chain, segment="V", limit=limit
        )
        sub_j = self.get_gene(
            source=source, species=species, chain=chain, segment="J", limit=limit
        )

        sequences = []
        for merge in product(sub_v, sub_j):

            v_seg = merge[0]
            j_seg = merge[1]

            functional = v_seg["imgt"]["imgt_functional"]
            v_part = v_seg["imgt"]["sequence_gapped_aa"].replace(".", "-")
            cdr3_part = j_seg["imgt"]["cdr3_aa"]
            fwr4_part = j_seg["imgt"]["fwr4_aa"]
            v_name = v_seg["gene"]
            j_name = j_seg["gene"]

            name = f"{species}_{v_name}_{j_name}"

            # why?
            if functional != "F":
                continue

            # H rules are default

            # K rules
            if chain.strip().lower() in "k":
                if len(cdr3_part[-3:] + fwr4_part) == 12:
                    fwr4_part += "-"
                v_part = v_part[:-3]

            # L rules
            if chain.strip().lower() == "l":
                if len(cdr3_part[-3:] + fwr4_part) == 12:
                    fwr4_part += "-"
                v_part = v_part[:-5]

            # todo: alt fwr4_part based on it's size and who's askin
            multiplier = 128 - (len(v_part) + len(cdr3_part[-3:] + fwr4_part))

            align = v_part + "-" * multiplier + cdr3_part[-3:] + fwr4_part

            # sanity check if chains rules are working
            assert len(align) == 128

            sequence = pyhmmer.easel.TextSequence(name=name.encode(), sequence=align)
            sequences.append(sequence)

        if not sequences:
            return None
        return pyhmmer.easel.TextMSA(name=f"{species}_{chain}".encode(), sequences=sequences).digitize(self.alphabet)
    
    @lru_cache(maxsize=None)
    def get_hmm(
        self,
        source: str = "imgt",
        chain: str = "h",
        species: str = "human",
        limit: Optional[int] = None,
    ) -> str:
        if species in ['pig']:
            return None
        if species not in self.species:
            return None
        msa = self.get_msa(
            source=source, species=species, chain=chain, limit=limit
        )
        # Chain not found for that species
        if msa is None:
            return None
        hmm, _, _ = self.builder.build_msa(msa, self.background)
        return hmm