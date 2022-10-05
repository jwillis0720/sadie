import math
from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning
from Bio.SubsMat import FreqTable as FreqTable
from typing import Any

log = math.log
NOTYPE: int
ACCREP: int
OBSFREQ: int
SUBS: int
EXPFREQ: int
LO: int
EPSILON: float

class SeqMat(dict):
    alphabet: Any
    ab_list: Any
    mat_name: Any
    sum_letters: Any
    relative_entropy: int
    def __init__(
        self, data: Any | None = ..., alphabet: Any | None = ..., mat_name: str = ..., build_later: int = ...
    ) -> None: ...
    entropy: int
    def make_entropy(self) -> None: ...
    def sum(self): ...
    def format(
        self,
        fmt: str = ...,
        letterfmt: str = ...,
        alphabet: Any | None = ...,
        non_sym: Any | None = ...,
        full: bool = ...,
    ): ...
    def __sub__(self, other): ...
    def __mul__(self, other): ...
    def __rmul__(self, other): ...
    def __add__(self, other): ...

class SubstitutionMatrix(SeqMat):
    def calculate_relative_entropy(self, obs_freq_mat): ...

class LogOddsMatrix(SeqMat):
    def calculate_relative_entropy(self, obs_freq_mat): ...

def make_log_odds_matrix(
    acc_rep_mat,
    exp_freq_table: Any | None = ...,
    logbase: int = ...,
    factor: float = ...,
    round_digit: int = ...,
    keep_nd: int = ...,
): ...
def observed_frequency_to_substitution_matrix(obs_freq_mat): ...
def read_text_matrix(data_file): ...

diagNO: int
diagONLY: int
diagALL: int

def two_mat_relative_entropy(mat_1, mat_2, logbase: int = ..., diag=...): ...
def two_mat_correlation(mat_1, mat_2): ...
def two_mat_DJS(mat_1, mat_2, pi_1: float = ..., pi_2: float = ...): ...
