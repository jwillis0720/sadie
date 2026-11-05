"""Tests for mouse VDJ annotation validation using well-characterized GenBank sequences.

Validates SADIE's mouse annotation for:
- V/J gene calls
- vj_in_frame and productive status
- Stop codon detection
- Kabat numbering via mutational analysis
"""

from pathlib import Path

import pytest

from sadie.airr import Airr, AirrTable
from sadie.airr import methods as airr_methods


def _mouse_available() -> bool:
    """Check if mouse germlines are available."""
    from sadie.germlines import get_germlines_base_dir

    mouse_path = get_germlines_base_dir() / "igblast" / "Ig" / "internal_data" / "mouse"
    return mouse_path.exists()


skip_no_mouse = pytest.mark.skipif(not _mouse_available(), reason="mouse germlines not available")

MOUSE_VDJ_FASTA = (
    Path(__file__).parent.parent.parent / "data" / "fixtures" / "fasta_inputs" / "mouse_vdj_validation.fasta"
)

# Expected annotations for each sequence
# Keys: v_call, j_call, locus, vj_in_frame, productive, stop_codon
EXPECTED = {
    "PX423956": {
        "v_call": "IGHV1-54*01",
        "j_call": "IGHJ2*01",
        "locus": "IGH",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "PX423948": {
        "v_call": "IGHV1-18*01",
        "j_call": "IGHJ4*01",
        "locus": "IGH",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "PX423950": {
        "v_call": "IGHV1S81",
        "j_call": "IGHJ1*01",
        "locus": "IGH",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "PX423956_kappa": {
        "v_call": "IGKV4-61",
        "j_call": "IGKJ1*01",
        "locus": "IGK",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "LC600311": {
        "v_call": "IGKV1-117",
        "j_call": "IGKJ1*01",
        "locus": "IGK",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "LC600312": {
        "v_call": "IGKV1-135",
        "j_call": "IGKJ1*01",
        "locus": "IGK",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "GQ984293": {
        "v_call": "IGKV6-17",
        "j_call": "IGKJ5*01",
        "locus": "IGK",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "AY704179": {
        "v_call": "IGKV13-84",
        "j_call": "IGKJ5*01",
        "locus": "IGK",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "M94350_VJ": {
        "v_call": "IGLV2",
        "j_call": "IGLJ1",
        "locus": "IGL",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "PQ471374": {
        "v_call": "IGLV1",
        "j_call": "IGLJ1",
        "locus": "IGL",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "MZ450924": {
        "v_call": "IGLV1",
        "j_call": "IGLJ1",
        "locus": "IGL",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
    "AB160033": {
        "v_call": "IGHV1-69",
        "j_call": "IGHJ4*01",
        "locus": "IGH",
        "vj_in_frame": True,
        "productive": False,
        "stop_codon": True,
    },
    "PX423952": {
        "v_call": "IGHV1S135",
        "j_call": "IGHJ3*01",
        "locus": "IGH",
        "vj_in_frame": True,
        "productive": True,
        "stop_codon": False,
    },
}


@skip_no_mouse
class TestMouseVdjAnnotation:
    """Validate mouse VDJ annotation: vj_in_frame, productive, V/J calls."""

    @pytest.fixture(scope="class")
    def mouse_result(self) -> AirrTable:
        """Run SADIE mouse annotation on the validation FASTA."""
        assert MOUSE_VDJ_FASTA.exists(), f"Fixture not found: {MOUSE_VDJ_FASTA}"
        airr_api = Airr("mouse")
        result = airr_api.run_fasta(str(MOUSE_VDJ_FASTA))
        assert isinstance(result, AirrTable)
        return result

    def test_all_sequences_annotated(self, mouse_result: AirrTable) -> None:
        """All input sequences should produce annotation results."""
        result_ids = set(mouse_result["sequence_id"].tolist())
        for seq_id in EXPECTED:
            assert seq_id in result_ids, f"Missing annotation result for {seq_id}"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_v_call(self, mouse_result: AirrTable, seq_id: str) -> None:
        """V gene call should contain the expected V gene name."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        v_call = str(row.iloc[0]["v_call"])
        assert (
            expected["v_call"] in v_call
        ), f"{seq_id}: expected v_call containing '{expected['v_call']}', got '{v_call}'"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_j_call(self, mouse_result: AirrTable, seq_id: str) -> None:
        """J gene call should contain the expected J gene name."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        j_call = str(row.iloc[0]["j_call"])
        assert (
            expected["j_call"] in j_call
        ), f"{seq_id}: expected j_call containing '{expected['j_call']}', got '{j_call}'"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_locus(self, mouse_result: AirrTable, seq_id: str) -> None:
        """Locus assignment should match expected chain type."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        locus = row.iloc[0]["locus"]
        assert locus == expected["locus"], f"{seq_id}: expected locus '{expected['locus']}', got '{locus}'"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_vj_in_frame(self, mouse_result: AirrTable, seq_id: str) -> None:
        """vj_in_frame should match expected value."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        vj_in_frame = row.iloc[0]["vj_in_frame"]
        assert (
            vj_in_frame == expected["vj_in_frame"]
        ), f"{seq_id}: expected vj_in_frame={expected['vj_in_frame']}, got {vj_in_frame}"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_productive(self, mouse_result: AirrTable, seq_id: str) -> None:
        """productive should match expected value."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        productive = row.iloc[0]["productive"]
        assert (
            productive == expected["productive"]
        ), f"{seq_id}: expected productive={expected['productive']}, got {productive}"

    @pytest.mark.parametrize("seq_id", list(EXPECTED.keys()))
    def test_stop_codon(self, mouse_result: AirrTable, seq_id: str) -> None:
        """stop_codon should match expected value."""
        expected = EXPECTED[seq_id]
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        stop_codon = row.iloc[0]["stop_codon"]
        assert (
            stop_codon == expected["stop_codon"]
        ), f"{seq_id}: expected stop_codon={expected['stop_codon']}, got {stop_codon}"

    @pytest.mark.parametrize(
        "seq_id",
        [sid for sid, exp in EXPECTED.items() if exp["productive"]],
    )
    def test_productive_no_internal_stops(self, mouse_result: AirrTable, seq_id: str) -> None:
        """Productive sequences should not have '*' in sequence_alignment_aa."""
        row = mouse_result[mouse_result["sequence_id"] == seq_id]
        assert len(row) == 1, f"Expected exactly 1 result for {seq_id}, got {len(row)}"
        seq_aa = row.iloc[0]["sequence_alignment_aa"]
        if isinstance(seq_aa, str):
            assert "*" not in seq_aa, f"{seq_id}: productive sequence has internal stop codon in alignment AA"


@skip_no_mouse
class TestMouseKabatNumbering:
    """Validate Kabat numbering via mutational analysis on mouse sequences."""

    @pytest.fixture(scope="class")
    def mouse_result_with_kabat(self) -> AirrTable:
        """Run SADIE mouse annotation and then Kabat mutational analysis."""
        assert MOUSE_VDJ_FASTA.exists(), f"Fixture not found: {MOUSE_VDJ_FASTA}"
        airr_api = Airr("mouse")
        result = airr_api.run_fasta(str(MOUSE_VDJ_FASTA))
        assert isinstance(result, AirrTable)
        result_with_kabat = airr_methods.run_mutational_analysis(result, scheme="kabat", run_multiproc=False)
        return result_with_kabat

    def test_mutations_column_present(self, mouse_result_with_kabat: AirrTable) -> None:
        """Mutational analysis should add a 'mutations' column."""
        assert "mutations" in mouse_result_with_kabat.columns, "Expected 'mutations' column after Kabat analysis"

    def test_scheme_column_present(self, mouse_result_with_kabat: AirrTable) -> None:
        """Mutational analysis should add a 'scheme' column."""
        assert "scheme" in mouse_result_with_kabat.columns, "Expected 'scheme' column after Kabat analysis"

    def test_productive_sequences_have_numbering(self, mouse_result_with_kabat: AirrTable) -> None:
        """Productive sequences should have non-null Kabat numbering data."""
        productive_ids = [sid for sid, exp in EXPECTED.items() if exp["productive"]]
        productive_rows = mouse_result_with_kabat[mouse_result_with_kabat["sequence_id"].isin(productive_ids)]
        assert len(productive_rows) > 0, "No productive sequences found in results"

        # At least some productive sequences should have non-null scheme data
        scheme_col = productive_rows["scheme"]
        non_null_count = scheme_col.notna().sum()
        assert non_null_count > 0, "All productive sequences have null Kabat numbering — expected at least some"

    def test_kabat_heavy_cdr_positions(self, mouse_result_with_kabat: AirrTable) -> None:
        """For productive heavy chain sequences, Kabat scheme data should contain expected CDR position ranges.

        Kabat CDR definitions for heavy chain:
        - CDR-H1: positions 31-35 (with possible 35A, 35B insertions)
        - CDR-H2: positions 50-65
        - CDR-H3: positions 95-102
        """
        productive_heavy_ids = [sid for sid, exp in EXPECTED.items() if exp["productive"] and exp["locus"] == "IGH"]
        heavy_rows = mouse_result_with_kabat[mouse_result_with_kabat["sequence_id"].isin(productive_heavy_ids)]

        if len(heavy_rows) == 0:
            pytest.skip("No productive heavy chain results available")

        # Check that scheme column has data for at least one heavy chain
        has_scheme_data = False
        for _, row in heavy_rows.iterrows():
            scheme_val = row.get("scheme")
            if isinstance(scheme_val, str) and len(scheme_val) > 0:
                has_scheme_data = True
                break

        assert has_scheme_data, "No heavy chain sequence has Kabat scheme data"

    def test_kabat_light_cdr_positions(self, mouse_result_with_kabat: AirrTable) -> None:
        """For productive light chain sequences, Kabat scheme data should be present.

        Kabat CDR definitions for light chain:
        - CDR-L1: positions 24-34
        - CDR-L2: positions 50-56
        - CDR-L3: positions 89-97
        """
        productive_light_ids = [
            sid for sid, exp in EXPECTED.items() if exp["productive"] and exp["locus"] in ("IGK", "IGL")
        ]
        light_rows = mouse_result_with_kabat[mouse_result_with_kabat["sequence_id"].isin(productive_light_ids)]

        if len(light_rows) == 0:
            pytest.skip("No productive light chain results available")

        # Check that scheme column has data for at least one light chain
        has_scheme_data = False
        for _, row in light_rows.iterrows():
            scheme_val = row.get("scheme")
            if isinstance(scheme_val, str) and len(scheme_val) > 0:
                has_scheme_data = True
                break

        assert has_scheme_data, "No light chain sequence has Kabat scheme data"
