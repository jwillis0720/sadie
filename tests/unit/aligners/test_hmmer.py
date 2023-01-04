from pathlib import Path

import pyhmmer
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sadie.renumbering.aligners import HMMER


class TestHMMER:
    def setup_class(self):
        self.hmmer = HMMER()
        self.aa = "ACDEFGHIKLMNPQRSTVWY"

    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
            (Seq("ACDEFGHIKLMNPQRSTVWY"), "ACDEFGHIKLMNPQRSTVWY"),
        ],
    )
    def test_digitize_seq(self, seq: Seq | str, expected: str) -> None:
        assert self.hmmer.digitize_seq("name", seq).textize().sequence == expected

    def test_transform_seq(self, fixture_setup) -> None:
        fastaa = fixture_setup.get_catnap_heavy_aa()
        seqs = self.hmmer.transform_seqs(Path(fastaa))
        assert seqs is not None
        seqs = self.hmmer.transform_seqs(str(fastaa))
        assert seqs is not None
        seqs = self.hmmer.transform_seqs(self.aa)
        assert seqs[0].textize().sequence == self.aa
        seqs = self.hmmer.transform_seqs(SeqRecord(self.aa, id="test"))
        assert seqs[0].textize().sequence == self.aa
        seqs = self.hmmer.transform_seqs((self.aa, "test"))
        assert seqs[0].textize().sequence == self.aa
        with pytest.raises(ValueError):
            self.hmmer.transform_seqs(None)

    def test_get_hmm_models(self):
        hmms = self.hmmer.get_hmm_models(species=["macaque"], chains=["H"], source="imgt")
        assert hmms is not None
        hmms = self.hmmer.get_hmm_models(species=["pig"], chains=["H"], source="imgt")
        assert hmms is not None

    def test_hmmsearch(self):
        assert self.hmmer.hmmsearch("ACDEFGHIKLMNPQRSTVWY") == [[]]
        assert self.hmmer.hmmsearch("ACDEFGHIKLMNPQRSTVWY", for_numbering=True) == [
            ([["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"]], [], [])
        ]
