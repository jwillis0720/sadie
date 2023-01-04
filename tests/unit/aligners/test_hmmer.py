import pytest
from Bio.Seq import Seq

from sadie.renumbering.aligners import HMMER


class TestHMMER:
    def setup_class(self):
        self.hmmer = HMMER()

    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
            (Seq("ACDEFGHIKLMNPQRSTVWY"), "ACDEFGHIKLMNPQRSTVWY"),
        ],
    )
    def test_digitize_seq(self, seq: Seq | str, expected: str) -> None:
        assert self.hmmer.digitize_seq("name", seq).textize().sequence == expected
