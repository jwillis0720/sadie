from __future__ import annotations

import copy
import itertools
import os
from functools import lru_cache
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest
from hypothesis import HealthCheck, Phase, given, settings
from hypothesis import strategies as st
from Levenshtein import distance as lev_distance
from rapidfuzz.process import cdist as rapidfuzz_cdist
from sklearn.metrics import pairwise_distances

from sadie.airr import AirrTable, LinkedAirrTable
from sadie.cluster import Cluster

_DEFAULT_LOOKUP = ["cdr1_aa", "cdr2_aa", "cdr3_aa"]
_UNLINKED_LOOKUPS = [
    _DEFAULT_LOOKUP,
    ["cdr3_aa"],
    ["cdr3_aa", "cdr1_aa"],
    ["cdr1_aa", "cdr1_aa"],
]
_LINKED_LOOKUPS = [
    _DEFAULT_LOOKUP,
    ["cdr3_aa_heavy"],
    ["cdr3_aa_light", "cdr1_aa_heavy"],
    ["cdr1_aa_heavy", "cdr1_aa_heavy"],
]
_LOOKUP_COLUMNS = [
    "cdr1_aa_heavy",
    "cdr2_aa_heavy",
    "cdr3_aa_heavy",
    "cdr1_aa_light",
    "cdr2_aa_light",
    "cdr3_aa_light",
]
_GREEDY_ENABLED = os.environ.get("SADIE_CLUSTER_GREEDY") == "1"
_GREEDY_EXAMPLES = int(os.environ.get("SADIE_CLUSTER_GREEDY_EXAMPLES", "40000"))
_GREEDY_THREADED_EXAMPLES = int(os.environ.get("SADIE_CLUSTER_THREADED_EXAMPLES", "200"))
_GREEDY_PROGRESS = 0
_GREEDY_THREADED_PROGRESS = 0
pytestmark = pytest.mark.filterwarnings("ignore:DataFrame columns are not unique")


class _LegacyCluster(Cluster):
    """The exact distance kernel from main at 973395e4."""

    def _get_distance_df(self, df: pd.DataFrame) -> Any:
        if self.pad_somatic:
            _lookup = self.lookup + self.pad_somatic_values
            if self.include_only_v_gene:
                if self._type == "linked":
                    df["mutations_heavy"] = df["mutations_heavy"].apply(self._get_v_gene_only)
                    df["mutations_light"] = df["mutations_light"].apply(self._get_v_gene_only)
                else:
                    df["mutations"] = df["mutations"].apply(self._get_v_gene_only)
        else:
            _lookup = self.lookup
        df_lookup = df[_lookup].to_dict(orient="index")

        def calc_lev(x: npt.ArrayLike, y: npt.ArrayLike) -> float:
            dist = 0
            for metric in self.lookup:
                dist += lev_distance(str(df_lookup[x[0]][metric]), str(df_lookup[y[0]][metric]))  # type: ignore[index]
            if self.pad_somatic and x[0] != y[0]:  # type: ignore[index]
                if len(self.pad_somatic_values) == 2:
                    mutations_1_heavy = df_lookup[x[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    mutations_2_heavy = df_lookup[y[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    mutations_1_light = df_lookup[x[0]][self.pad_somatic_values[1]]  # type: ignore[index]
                    mutations_2_light = df_lookup[y[0]][self.pad_somatic_values[1]]  # type: ignore[index]
                    subtract_all = len(np.intersect1d(mutations_1_heavy, mutations_2_heavy)) + len(
                        np.intersect1d(mutations_1_light, mutations_2_light)
                    )
                else:
                    mutations_1 = df_lookup[x[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    mutations_2 = df_lookup[y[0]][self.pad_somatic_values[0]]  # type: ignore[index]
                    subtract_all = len(np.intersect1d(mutations_1, mutations_2))
                dist -= subtract_all
            return max(dist, 0)

        X: npt.ArrayLike = np.array(df.index).reshape(-1, 1)
        return pairwise_distances(X, metric=calc_lev, n_jobs=-1)


class _WeirdMutation(str):
    def __hash__(self) -> int:
        return 1

    def __eq__(self, other: object) -> bool:
        return True


class _ChangingLookup:
    def __init__(self):
        self.calls = 0

    def __str__(self) -> str:
        self.calls += 1
        return "A" if self.calls % 2 else "B"


@st.composite
def _mutation_token(draw: st.DrawFn) -> str:
    amino_acid = st.sampled_from(list("ACDEFGHIKLMNPQRSTVWYÅ"))
    position = st.one_of(st.sampled_from([0, 1, 92, 93, 94, 95, 120]), st.integers(min_value=0, max_value=150))
    return f"{draw(amino_acid)}{draw(position)}{draw(st.sampled_from(['', '.1']))}{draw(amino_acid)}"


_LOOKUP_VALUE = st.one_of(
    st.text(max_size=12),
    st.none(),
    st.booleans(),
    st.integers(min_value=-1000, max_value=1000),
    st.floats(width=32, allow_nan=True, allow_infinity=True),
    st.just(pd.NA),
)
_MUTATION_SPEC = st.tuples(
    st.sampled_from(["list", "array"]),
    st.lists(_mutation_token(), min_size=0, max_size=7).map(tuple),
)
_BOUNDARY_ATOM = st.one_of(
    st.text(max_size=6),
    st.binary(max_size=6),
    st.none(),
    st.booleans(),
    st.integers(min_value=-10, max_value=10),
    st.floats(width=32, allow_nan=True, allow_infinity=True),
    st.just(pd.NA),
)


@st.composite
def _boundary_mutation_spec(draw: st.DrawFn) -> tuple[str, Any]:
    kind = draw(st.sampled_from(["scalar", "list", "tuple", "set", "dict", "array0", "array1", "array2"]))
    if kind in {"scalar", "array0"}:
        value = draw(_BOUNDARY_ATOM)
    elif kind == "set":
        value = draw(st.sets(_BOUNDARY_ATOM, max_size=4))
    elif kind == "dict":
        value = draw(st.dictionaries(st.text(max_size=4), _BOUNDARY_ATOM, max_size=4))
    elif kind == "array2":
        value = draw(st.lists(_BOUNDARY_ATOM, min_size=4, max_size=4).map(tuple))
    else:
        value = draw(st.lists(_BOUNDARY_ATOM, max_size=5).map(tuple))
    return kind, value


@st.composite
def _cluster_case(draw: st.DrawFn, min_rows: int = 2, max_rows: int = 8) -> dict[str, Any]:
    row_count = draw(st.integers(min_value=min_rows, max_value=max_rows))
    lookup_rows = draw(
        st.lists(
            st.lists(_LOOKUP_VALUE, min_size=6, max_size=6).map(tuple),
            min_size=row_count,
            max_size=row_count,
        )
    )
    mutation_heavy = draw(st.lists(_MUTATION_SPEC, min_size=row_count, max_size=row_count))
    mutation_light = draw(st.lists(_MUTATION_SPEC, min_size=row_count, max_size=row_count))
    groups = draw(st.lists(st.sampled_from(["g0", "g1", "g2", "", "Å"]), min_size=row_count, max_size=row_count))
    batches = draw(st.lists(st.sampled_from(["b0", "b1", ""]), min_size=row_count, max_size=row_count))
    index_kind = draw(st.sampled_from(["range", "integer", "string"]))
    index_start = draw(st.integers(min_value=-1000, max_value=1000))
    if index_kind == "range":
        index = list(range(index_start, index_start + row_count))
    elif index_kind == "integer":
        step = draw(st.integers(min_value=2, max_value=9))
        index = [index_start + step * row for row in range(row_count)]
    else:
        prefix = draw(st.text(alphabet=list("abÅ_"), max_size=4))
        index = [f"{prefix}:{index_start + row}" for row in range(row_count)]

    return {
        "linked": draw(st.booleans()),
        "lookup_rows": lookup_rows,
        "mutation_heavy": mutation_heavy,
        "mutation_light": mutation_light,
        "groups": groups,
        "batches": batches,
        "index": index,
        "lookup_choice": draw(st.integers(min_value=0, max_value=3)),
        "pad_somatic": draw(st.booleans()),
        "include_only_v_gene": draw(st.booleans()),
        "groupby_choice": draw(st.integers(min_value=0, max_value=2)),
        "linkage": draw(st.sampled_from(["single", "average", "complete"])),
        "distance_threshold": draw(st.integers(min_value=0, max_value=15)),
    }


@st.composite
def _threaded_cluster_case(draw: st.DrawFn) -> dict[str, Any]:
    row_count = draw(st.sampled_from([149, 150, 151]))
    case = _fixed_case(draw(st.booleans()), row_count)
    lookup_palette = draw(st.lists(st.tuples(*[_LOOKUP_VALUE for _ in range(6)]), min_size=1, max_size=4))
    heavy_palette = draw(st.lists(_MUTATION_SPEC, min_size=1, max_size=4))
    light_palette = draw(st.lists(_MUTATION_SPEC, min_size=1, max_size=4))
    case.update(
        {
            "lookup_rows": [lookup_palette[row % len(lookup_palette)] for row in range(row_count)],
            "mutation_heavy": [heavy_palette[row % len(heavy_palette)] for row in range(row_count)],
            "mutation_light": [light_palette[row % len(light_palette)] for row in range(row_count)],
            "lookup_choice": draw(st.integers(min_value=0, max_value=3)),
            "pad_somatic": draw(st.booleans()),
            "include_only_v_gene": draw(st.booleans()),
        }
    )
    return case


def _materialize_mutation(spec: tuple[str, tuple[str, ...]]) -> list[str] | np.ndarray:
    kind, mutations = spec
    return list(mutations) if kind == "list" else np.array(mutations, dtype=object)


def _materialize_boundary_mutation(spec: tuple[str, Any]) -> Any:
    kind, value = spec
    if kind == "scalar":
        return copy.deepcopy(value)
    if kind == "list":
        return list(copy.deepcopy(value))
    if kind == "tuple":
        return tuple(copy.deepcopy(value))
    if kind == "set":
        return set(copy.deepcopy(value))
    if kind == "dict":
        return dict(copy.deepcopy(value))
    if kind == "array0":
        return np.array(copy.deepcopy(value), dtype=object)
    if kind == "array1":
        return np.array(copy.deepcopy(value), dtype=object)
    return np.array(copy.deepcopy(value), dtype=object).reshape(2, 2)


def _lookup(case: dict[str, Any]) -> list[str]:
    choices = _LINKED_LOOKUPS if case["linked"] else _UNLINKED_LOOKUPS
    return list(choices[case["lookup_choice"]])


def _resolved_lookup(case: dict[str, Any]) -> list[str]:
    lookup = _lookup(case)
    if case["linked"] and lookup == _DEFAULT_LOOKUP:
        return [f"{column}_heavy" for column in lookup] + [f"{column}_light" for column in lookup]
    return lookup


def _groupby(case: dict[str, Any]) -> list[str] | None:
    return [None, ["group"], ["group", "batch"]][case["groupby_choice"]]


def _distance_frame(case: dict[str, Any]) -> pd.DataFrame:
    linked = case["linked"]
    lookup_columns = _LOOKUP_COLUMNS if linked else _DEFAULT_LOOKUP
    lookup_indexes = range(6) if linked else range(3)
    data: dict[str, Any] = {
        column: [row[column_index] for row in case["lookup_rows"]]
        for column, column_index in zip(lookup_columns, lookup_indexes)
    }
    if linked:
        data["mutations_heavy"] = [_materialize_mutation(spec) for spec in case["mutation_heavy"]]
        data["mutations_light"] = [_materialize_mutation(spec) for spec in case["mutation_light"]]
        data["cellid"] = [f"cell-{row}" for row in range(len(case["lookup_rows"]))]
    else:
        data["mutations"] = [_materialize_mutation(spec) for spec in case["mutation_heavy"]]
        data["sequence_id"] = [f"sequence-{row}" for row in range(len(case["lookup_rows"]))]
    data["group"] = case["groups"]
    data["batch"] = case["batches"]
    return pd.DataFrame(data, index=case["index"], dtype=object)


@lru_cache(maxsize=1)
def _linked_template() -> pd.DataFrame:
    fixture = Path(__file__).resolve().parents[2] / "data/fixtures/airr_tables/joined_airrtable_with_mutational.feather"
    return pd.read_feather(fixture)


def _table(case: dict[str, Any]) -> AirrTable | LinkedAirrTable:
    generated = _distance_frame(case)
    if not case["linked"]:
        return AirrTable(generated, key_column="sequence_id")

    row_count = len(generated)
    source = _linked_template()
    table = source.iloc[np.arange(row_count) % len(source)].copy(deep=True)
    table.index = generated.index
    for column in generated.columns:
        values = np.empty(row_count, dtype=object)
        for row, value in enumerate(generated[column]):
            values[row] = copy.deepcopy(value)
        table[column] = values
    return LinkedAirrTable(table, key_column="cellid")


def _distance_cluster(cluster_type: type[Cluster], case: dict[str, Any]) -> Cluster:
    cluster = cluster_type.__new__(cluster_type)
    cluster.lookup = _resolved_lookup(case)
    cluster.pad_somatic = case["pad_somatic"]
    cluster.include_only_v_gene = case["include_only_v_gene"]
    cluster._type = "linked" if case["linked"] else "unlinked"
    cluster.pad_somatic_values = ["mutations_heavy", "mutations_light"] if case["linked"] else ["mutations"]
    return cluster


def _assert_distance_equivalent(case: dict[str, Any]) -> None:
    legacy_frame = _distance_frame(case)
    current_frame = _distance_frame(case)
    legacy = _distance_cluster(_LegacyCluster, case)._get_distance_df(legacy_frame)
    current = _distance_cluster(Cluster, case)._get_distance_df(current_frame)

    np.testing.assert_array_equal(current, legacy)
    pd.testing.assert_frame_equal(current_frame, legacy_frame)
    assert current.dtype == np.float64
    assert current.shape == (len(current_frame), len(current_frame))
    np.testing.assert_array_equal(current, current.T)
    np.testing.assert_array_equal(np.diag(current), np.zeros(len(current_frame)))
    assert np.all(current >= 0)
    assert np.all(current == np.floor(current))


def _run_cluster(
    cluster_type: type[Cluster], case: dict[str, Any]
) -> tuple[AirrTable | LinkedAirrTable, Cluster, AirrTable | LinkedAirrTable]:
    input_table = _table(case)
    cluster = cluster_type(
        input_table,
        linkage=case["linkage"],
        groupby=_groupby(case),
        lookup=_lookup(case),
        pad_somatic=case["pad_somatic"],
        include_only_v_gene=case["include_only_v_gene"],
    )
    result = cluster.cluster(case["distance_threshold"])
    return input_table, cluster, result


def _assert_model_equivalent(current: Any, legacy: Any) -> None:
    assert (current is None) == (legacy is None)
    if current is None:
        return
    assert type(current) is type(legacy)
    for attribute in ("labels_", "children_", "distances_"):
        np.testing.assert_array_equal(getattr(current, attribute), getattr(legacy, attribute))
    assert current.n_leaves_ == legacy.n_leaves_
    assert current.n_connected_components_ == legacy.n_connected_components_


def _assert_cluster_equivalent(case: dict[str, Any]) -> None:
    legacy_input, legacy_cluster, legacy_result = _run_cluster(_LegacyCluster, case)
    current_input, current_cluster, current_result = _run_cluster(Cluster, case)

    assert type(current_result) is type(legacy_result)
    assert current_result.key_column == legacy_result.key_column
    assert current_result.verified == legacy_result.verified
    if isinstance(current_result, LinkedAirrTable):
        assert current_result.suffixes == legacy_result.suffixes
    pd.testing.assert_frame_equal(pd.DataFrame(current_result), pd.DataFrame(legacy_result), check_exact=True)
    pd.testing.assert_frame_equal(pd.DataFrame(current_input), pd.DataFrame(legacy_input), check_exact=True)
    pd.testing.assert_frame_equal(
        pd.DataFrame(current_cluster.airrtable), pd.DataFrame(legacy_cluster.airrtable), check_exact=True
    )
    np.testing.assert_array_equal(current_cluster.distance_df, legacy_cluster.distance_df)
    _assert_model_equivalent(current_cluster.model, legacy_cluster.model)
    assert (current_result is current_cluster.airrtable) == (legacy_result is legacy_cluster.airrtable)
    assert (current_result is current_input) == (legacy_result is legacy_input)


def _fixed_case(linked: bool, row_count: int = 6) -> dict[str, Any]:
    return {
        "linked": linked,
        "lookup_rows": [
            (f"A{row}", f"B{row % 3}", f"C{row % 2}", f"D{row}", f"E{row % 3}", f"F{row % 2}")
            for row in range(row_count)
        ],
        "mutation_heavy": [
            ("array" if row % 2 else "list", ("A0T", "C93G", "G94A", "A0T")[: row % 5]) for row in range(row_count)
        ],
        "mutation_light": [
            ("list" if row % 2 else "array", ("T1C", "G93A", "C120T")[: row % 4]) for row in range(row_count)
        ],
        "groups": [f"g{row % 3}" for row in range(row_count)],
        "batches": [f"b{row % 2}" for row in range(row_count)],
        "index": [10 + row * 3 for row in range(row_count)],
        "lookup_choice": 0,
        "pad_somatic": False,
        "include_only_v_gene": False,
        "groupby_choice": 0,
        "linkage": "complete",
        "distance_threshold": 3,
    }


def _assert_raw_mutation_outcome_equivalent(mutations: list[Any]) -> None:
    case = _fixed_case(False, len(mutations))
    case["lookup_choice"] = 1
    case["pad_somatic"] = True

    def frame() -> pd.DataFrame:
        values = np.empty(len(mutations), dtype=object)
        for row, mutation in enumerate(mutations):
            values[row] = copy.deepcopy(mutation)
        return pd.DataFrame(
            {"cdr3_aa": ["AAAA" if row % 2 == 0 else "TTTT" for row in range(len(mutations))], "mutations": values},
            index=[10 + row * 3 for row in range(len(mutations))],
            dtype=object,
        )

    def outcome(cluster_type: type[Cluster], source: pd.DataFrame) -> tuple[str, Any]:
        try:
            return "result", _distance_cluster(cluster_type, case)._get_distance_df(source)
        except Exception as error:
            return "error", error

    legacy_frame = frame()
    current_frame = frame()
    legacy_outcome = outcome(_LegacyCluster, legacy_frame)
    current_outcome = outcome(Cluster, current_frame)

    assert current_outcome[0] == legacy_outcome[0]
    if current_outcome[0] == "result":
        np.testing.assert_array_equal(current_outcome[1], legacy_outcome[1])
        pd.testing.assert_frame_equal(current_frame, legacy_frame)
    else:
        # NumPy can reverse operand names in equivalent TypeError messages.
        assert type(current_outcome[1]) is type(legacy_outcome[1])


@given(case=_cluster_case(min_rows=1, max_rows=12))
@settings(max_examples=500, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_generated_distance_kernel_matches_legacy(case: dict[str, Any]):
    _assert_distance_equivalent(case)


@given(case=_cluster_case(min_rows=2, max_rows=8))
@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_generated_cluster_output_matches_legacy(case: dict[str, Any]):
    _assert_cluster_equivalent(case)


@given(left=_boundary_mutation_spec(), right=_boundary_mutation_spec(), singleton=st.booleans())
@settings(max_examples=500, deadline=None, suppress_health_check=[HealthCheck.too_slow])
def test_generated_mutation_boundaries_match_legacy(left: tuple[str, Any], right: tuple[str, Any], singleton: bool):
    mutations = [_materialize_boundary_mutation(left)]
    if not singleton:
        mutations.append(_materialize_boundary_mutation(right))
    _assert_raw_mutation_outcome_equivalent(mutations)


def test_bounded_exhaustive_lookup_values_match_legacy():
    values = ["".join(value) for length in range(5) for value in itertools.product("ABÅ", repeat=length)]
    values += ["é", "e\u0301", "🙂", "\0", "", " ", None, np.nan, pd.NA, -1, 0, 1]
    case = _fixed_case(False, len(values))
    case["lookup_choice"] = 1
    case["lookup_rows"] = [("", "", value, "", "", "") for value in values]
    _assert_distance_equivalent(case)


@pytest.mark.parametrize("representation", ["list", "array"])
@pytest.mark.parametrize("include_only_v_gene", [False, True])
def test_bounded_exhaustive_mutation_lists_match_legacy(representation: str, include_only_v_gene: bool):
    tokens = ("A0T", "C93G", "G94A", "T120C")
    mutation_lists = [mutations for length in range(4) for mutations in itertools.product(tokens, repeat=length)]
    row_count = len(mutation_lists) * 2
    case = _fixed_case(False, row_count)
    case["lookup_choice"] = 1
    case["lookup_rows"] = [
        ("", "", "AAAA" if block == 0 else "TTTT", "", "", "") for block in range(2) for _ in mutation_lists
    ]
    case["mutation_heavy"] = [(representation, mutations) for _ in range(2) for mutations in mutation_lists]
    case["pad_somatic"] = True
    case["include_only_v_gene"] = include_only_v_gene
    _assert_distance_equivalent(case)


def test_bounded_exhaustive_linked_chain_padding_matches_legacy():
    tokens = ("A0T", "C94G")
    mutation_lists = [mutations for length in range(3) for mutations in itertools.product(tokens, repeat=length)]
    combinations = list(itertools.product(mutation_lists, repeat=2))
    case = _fixed_case(True, len(combinations))
    case["lookup_choice"] = 1
    case["lookup_rows"] = [("", "", "AAAA" if row % 2 else "TTTT", "", "", "") for row in range(len(combinations))]
    case["mutation_heavy"] = [("list", heavy) for heavy, _ in combinations]
    case["mutation_light"] = [("array", light) for _, light in combinations]
    case["pad_somatic"] = True
    _assert_distance_equivalent(case)


@pytest.mark.parametrize("linked", [False, True])
@pytest.mark.greedy
@pytest.mark.skipif(not _GREEDY_ENABLED, reason="set SADIE_CLUSTER_GREEDY=1 to run exhaustive option combinations")
def test_cluster_option_cross_product_matches_legacy(linked: bool):
    for linkage, groupby_choice, lookup_choice, pad_somatic, include_only_v_gene, threshold in itertools.product(
        ["single", "average", "complete"],
        range(3),
        range(4),
        [False, True],
        [False, True],
        [0, 3, 12],
    ):
        case = _fixed_case(linked)
        case.update(
            {
                "linkage": linkage,
                "groupby_choice": groupby_choice,
                "lookup_choice": lookup_choice,
                "pad_somatic": pad_somatic,
                "include_only_v_gene": include_only_v_gene,
                "distance_threshold": threshold,
            }
        )
        _assert_cluster_equivalent(case)


@pytest.mark.parametrize("row_count", [149, 150, 151])
@pytest.mark.parametrize("pad_somatic", [False, True])
def test_thread_boundary_matches_legacy_and_is_deterministic(row_count: int, pad_somatic: bool):
    case = _fixed_case(False, row_count)
    case["pad_somatic"] = pad_somatic
    legacy_frame = _distance_frame(case)
    expected = _distance_cluster(_LegacyCluster, case)._get_distance_df(legacy_frame)

    results = []
    with patch("sadie.cluster.cluster.cdist", wraps=rapidfuzz_cdist) as cdist_spy:
        for _ in range(10):
            result = _distance_cluster(Cluster, case)._get_distance_df(_distance_frame(case))
            np.testing.assert_array_equal(result, expected)
            results.append(result)
    assert {call.kwargs["workers"] for call in cdist_spy.call_args_list} == ({-1} if row_count >= 150 else {1})
    for result in results[1:]:
        np.testing.assert_array_equal(result, results[0])


@pytest.mark.parametrize("linked", [False, True])
def test_grouped_threaded_cluster_matches_legacy(linked: bool):
    case = _fixed_case(linked, 300)
    case.update(
        {
            "groups": ["g0"] * 150 + ["g1"] * 150,
            "lookup_choice": 2,
            "pad_somatic": True,
            "include_only_v_gene": True,
            "groupby_choice": 1,
        }
    )
    with patch("sadie.cluster.cluster.cdist", wraps=rapidfuzz_cdist) as cdist_spy:
        _assert_cluster_equivalent(case)
    assert {call.kwargs["workers"] for call in cdist_spy.call_args_list} == {-1}


@pytest.mark.parametrize(
    ("left", "right"),
    [
        ("A1T", "A1T"),
        (b"A1T", b"A1T"),
        (1, 1),
        (np.nan, np.nan),
        (None, None),
        ({"A1T"}, {"A1T"}),
        ({"mutation": "A1T"}, {"mutation": "A1T"}),
        ([["A1T"]], [["A1T"]]),
        ([1, "1"], ["1"]),
        ([np.nan], [np.nan]),
        ([""], ["\0"]),
        (pd.NA, pd.NA),
        ([_WeirdMutation("A")], [_WeirdMutation("B")]),
    ],
)
def test_malformed_mutation_cells_preserve_legacy_outcome(left: Any, right: Any):
    _assert_raw_mutation_outcome_equivalent([left, right])


@pytest.mark.parametrize("mutation", ["A1T", b"A1T", 1, np.nan, None, {"A1T"}, [["A1T"]], pd.NA])
def test_singleton_malformed_mutation_cells_match_legacy(mutation: Any):
    _assert_raw_mutation_outcome_equivalent([mutation])


@pytest.mark.parametrize("pad_somatic_values", [[], ["mutations", "extra", "third"]])
def test_nonstandard_padding_columns_match_legacy(pad_somatic_values: list[str]):
    case = _fixed_case(False, 2)
    case.update({"lookup_choice": 1, "pad_somatic": True})
    frame = _distance_frame(case)
    frame["extra"] = [["A0T"], ["A0T"]]
    frame["third"] = [["A0T"], ["A0T"]]

    def outcome(cluster_type: type[Cluster]) -> tuple[str, Any]:
        cluster = _distance_cluster(cluster_type, case)
        cluster.pad_somatic_values = pad_somatic_values
        try:
            return "result", cluster._get_distance_df(frame.copy(deep=True))
        except Exception as error:
            return "error", error

    legacy = outcome(_LegacyCluster)
    current = outcome(Cluster)
    assert current[0] == legacy[0]
    if current[0] == "result":
        np.testing.assert_array_equal(current[1], legacy[1])
    else:
        assert type(current[1]) is type(legacy[1])


def test_stateful_lookup_objects_match_legacy():
    case = _fixed_case(False, 2)
    case["lookup_choice"] = 1

    def frame() -> pd.DataFrame:
        result = _distance_frame(case)
        result["cdr3_aa"] = [_ChangingLookup(), _ChangingLookup()]
        return result

    legacy = _distance_cluster(_LegacyCluster, case)._get_distance_df(frame())
    with patch("sadie.cluster.cluster.cdist", wraps=rapidfuzz_cdist) as cdist_spy:
        current = _distance_cluster(Cluster, case)._get_distance_df(frame())
    np.testing.assert_array_equal(current, legacy)
    cdist_spy.assert_not_called()


def test_empty_distance_matrix_matches_legacy():
    case = _fixed_case(False, 0)
    case["lookup_choice"] = 1
    try:
        legacy = _distance_cluster(_LegacyCluster, case)._get_distance_df(_distance_frame(case))
    except Exception as legacy_error:
        with pytest.raises(type(legacy_error)):
            _distance_cluster(Cluster, case)._get_distance_df(_distance_frame(case))
    else:
        current = _distance_cluster(Cluster, case)._get_distance_df(_distance_frame(case))
        np.testing.assert_array_equal(current, legacy)


@pytest.mark.parametrize(
    "index",
    [
        [0, np.nan],
        [0, np.inf],
        [0, pd.NaT],
        [0, 0],
        pd.Index([1, "1"], dtype=object),
        pd.Index([("a", 1), ("b", 2)]),
        pd.MultiIndex.from_tuples([("a", 1), ("b", 2)]),
        pd.date_range("2020-01-01", periods=2),
        pd.DatetimeIndex(["2020-01-01", pd.NaT]),
        pd.to_timedelta([1, 2], unit="D"),
        pd.TimedeltaIndex([pd.Timedelta(days=1), pd.NaT]),
    ],
)
def test_index_result_or_error_matches_legacy(index: Any):
    row_count = len(index)
    case = _fixed_case(False, row_count)
    case["lookup_choice"] = 1
    case["index"] = index
    legacy_frame = _distance_frame(case)
    current_frame = _distance_frame(case)

    def outcome(cluster_type: type[Cluster], source: pd.DataFrame) -> tuple[str, Any]:
        try:
            return "result", _distance_cluster(cluster_type, case)._get_distance_df(source)
        except Exception as error:
            return "error", error

    legacy_outcome = outcome(_LegacyCluster, legacy_frame)
    current_outcome = outcome(Cluster, current_frame)
    assert current_outcome[0] == legacy_outcome[0]
    if current_outcome[0] == "result":
        np.testing.assert_array_equal(current_outcome[1], legacy_outcome[1])
    else:
        assert type(current_outcome[1]) is type(legacy_outcome[1])


@pytest.mark.greedy
@pytest.mark.skipif(not _GREEDY_ENABLED, reason="set SADIE_CLUSTER_GREEDY=1 to fuzz the threaded boundary")
@given(case=_threaded_cluster_case())
@settings(
    max_examples=_GREEDY_THREADED_EXAMPLES,
    deadline=None,
    derandomize=True,
    database=None,
    phases=(Phase.generate, Phase.shrink),
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.large_base_example],
)
def test_greedy_threaded_distance_matches_legacy(case: dict[str, Any]):
    global _GREEDY_THREADED_PROGRESS
    _assert_distance_equivalent(case)
    _GREEDY_THREADED_PROGRESS += 1
    if _GREEDY_THREADED_PROGRESS % 10 == 0:
        print(
            f"greedy threaded equivalence: {_GREEDY_THREADED_PROGRESS}/{_GREEDY_THREADED_EXAMPLES}",
            flush=True,
        )


@pytest.mark.greedy
@pytest.mark.skipif(not _GREEDY_ENABLED, reason="set SADIE_CLUSTER_GREEDY=1 to run the multi-hour suite")
@given(
    case=_cluster_case(min_rows=2, max_rows=24),
    left=_boundary_mutation_spec(),
    right=_boundary_mutation_spec(),
)
@settings(
    max_examples=_GREEDY_EXAMPLES,
    deadline=None,
    derandomize=True,
    database=None,
    phases=(Phase.generate, Phase.shrink),
    suppress_health_check=[HealthCheck.too_slow],
)
def test_greedy_generated_cluster_and_distance_match_legacy(
    case: dict[str, Any], left: tuple[str, Any], right: tuple[str, Any]
):
    global _GREEDY_PROGRESS
    _assert_distance_equivalent(case)
    _assert_cluster_equivalent(case)
    _assert_raw_mutation_outcome_equivalent(
        [_materialize_boundary_mutation(left), _materialize_boundary_mutation(right)]
    )
    _GREEDY_PROGRESS += 1
    if _GREEDY_PROGRESS % 250 == 0:
        print(f"greedy cluster equivalence: {_GREEDY_PROGRESS}/{_GREEDY_EXAMPLES}", flush=True)
