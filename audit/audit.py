#!/usr/bin/env python
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pandas as pd

STATISTICAL_COLS = frozenset({
    "v_support", "d_support", "j_support",
})

ALLELE_DEPENDENT_COLS = frozenset({
    "v_call", "v_call_top", "d_call", "d_call_top", "j_call", "j_call_top",
    "v_germline_alignment", "v_germline_alignment_aa",
    "d_germline_alignment", "d_germline_alignment_aa",
    "j_germline_alignment", "j_germline_alignment_aa",
    "germline_alignment", "germline_alignment_aa",
    "v_identity", "v_mutation", "v_mutation_aa",
    "d_identity", "d_mutation", "d_mutation_aa",
    "j_identity", "j_mutation", "j_mutation_aa",
    "v_score", "d_score", "j_score",
    "v_cigar", "d_cigar", "j_cigar",
    "complete_vdj",
})

C_REGION_COLS = frozenset({
    "c_call", "c_score", "c_cigar", "c_support", "c_identity",
    "c_sequence_start", "c_sequence_end", "c_germline_start", "c_germline_end",
    "c_alignment_start", "c_alignment_end",
    "c_sequence_alignment", "c_sequence_alignment_aa",
    "c_germline_alignment", "c_germline_alignment_aa",
})


@dataclass
class AuditResult:
    germlines_df: pd.DataFrame
    g3_df: pd.DataFrame
    common_cols: list
    differences: dict = field(default_factory=dict)

    @property
    def c_region_cols(self) -> list:
        return [c for c in self.common_cols if c in C_REGION_COLS]

    @property
    def structural_cols(self) -> list:
        return [c for c in self.common_cols
                if c not in STATISTICAL_COLS
                and c not in C_REGION_COLS
                and c not in ALLELE_DEPENDENT_COLS]

    @property
    def total_values(self) -> int:
        return len(self.common_cols) * len(self.germlines_df)

    @property
    def total_diffs(self) -> int:
        return sum(self.differences.values())

    @property
    def parity_pct(self) -> float:
        if self.total_values == 0:
            return 100.0
        return (1 - self.total_diffs / self.total_values) * 100

    @property
    def structural_parity(self) -> float:
        structural_values = len(self.structural_cols) * len(self.germlines_df)
        if structural_values == 0:
            return 100.0
        structural_diff_count = sum(
            v for k, v in self.differences.items()
            if k not in C_REGION_COLS
            and k not in STATISTICAL_COLS
            and k not in ALLELE_DEPENDENT_COLS
        )
        return (1 - structural_diff_count / structural_values) * 100

    def get_diffs_by_category(self) -> tuple[dict, dict, dict, dict]:
        c_region = {k: v for k, v in self.differences.items() if k in C_REGION_COLS}
        statistical = {k: v for k, v in self.differences.items() if k in STATISTICAL_COLS}
        allele = {k: v for k, v in self.differences.items()
                  if k in ALLELE_DEPENDENT_COLS and k not in STATISTICAL_COLS}
        structural = {k: v for k, v in self.differences.items()
                      if k not in C_REGION_COLS
                      and k not in STATISTICAL_COLS
                      and k not in ALLELE_DEPENDENT_COLS}
        return c_region, statistical, allele, structural


def patch_germline_manager_for_imgt_only():
    from sadie.germlines.manager import GermlineManager
    original_init = GermlineManager.__init__

    def imgt_only_init(self, providers=None, data_dir=None):
        original_init(self, providers=["imgt"], data_dir=data_dir)

    GermlineManager.__init__ = imgt_only_init


def rebuild_databases(base_path: Path):
    print("=" * 60)
    print("REBUILDING GERMLINE DATABASES WITH IMGT ONLY")
    print("=" * 60)
    patch_germline_manager_for_imgt_only()
    from sadie.germlines.pipeline import GermlinePipeline
    pipeline = GermlinePipeline(base_path / "src" / "sadie" / "germlines")
    print("Rebuilding with IMGT-only sources...")
    pipeline.force_rebuild("human")
    print("Database rebuild complete!\n")


def load_sequences(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} sequences")
    sequences = df[["sequence_id_heavy", "sequence_heavy"]].dropna()
    sequences = sequences.rename(columns={
        "sequence_id_heavy": "sequence_id",
        "sequence_heavy": "sequence"
    })
    if sequences["sequence_id"].duplicated().any():
        dup_count = sequences["sequence_id"].duplicated().sum()
        print(f"Found {dup_count} duplicate sequence IDs, keeping first occurrence")
        sequences = sequences.drop_duplicates(subset="sequence_id", keep="first")
    print(f"Prepared {len(sequences)} unique sequences for annotation")
    return sequences


def run_backend(sequences: pd.DataFrame, use_germlines: bool) -> pd.DataFrame:
    from sadie.airr import Airr
    backend_name = "Germlines" if use_germlines else "G3"
    print(f"\n=== Running with {backend_name} Backend ===")
    os.environ["SADIE_USE_GERMLINES_MODULE"] = str(use_germlines).lower()
    airr = Airr("human")
    result = airr.run_dataframe(sequences, seq_id_field="sequence_id", seq_field="sequence")
    if "source" in result.columns:
        result = result.drop(columns=["source"])
    print(f"{backend_name} backend: {len(result)} results, {len(result.columns)} columns")
    return result


def normalize_series(series: pd.Series) -> pd.Series:
    return series.astype(str).replace("nan", "").replace("<NA>", "")


def find_differences(germlines_df: pd.DataFrame, g3_df: pd.DataFrame,
                     common_cols: list) -> dict:
    differences = {}
    for col in common_cols:
        g_vals = normalize_series(germlines_df[col])
        g3_vals = normalize_series(g3_df[col])
        diff_count = (g_vals != g3_vals).sum()
        if diff_count > 0:
            differences[col] = diff_count
    return differences


def compare_results(result_germlines: pd.DataFrame,
                    result_g3: pd.DataFrame) -> AuditResult:
    print("\n=== Comparing Results ===")
    germlines_cols = set(result_germlines.columns)
    g3_cols = set(result_g3.columns)

    print(f"Germlines columns: {len(germlines_cols)}")
    print(f"G3 columns: {len(g3_cols)}")
    print(f"Only in germlines: {germlines_cols - g3_cols}")
    print(f"Only in G3: {g3_cols - germlines_cols}")

    common_cols = sorted(germlines_cols & g3_cols)
    c_cols = [c for c in common_cols if c in C_REGION_COLS]
    print(f"C region columns present: {len(c_cols)}")
    print(f"Non-C region columns: {len(common_cols) - len(c_cols)}")
    print(f"Total common columns for comparison: {len(common_cols)}")

    germlines_aligned = result_germlines[common_cols].sort_values("sequence_id").reset_index(drop=True)
    g3_aligned = result_g3[common_cols].sort_values("sequence_id").reset_index(drop=True)

    differences = find_differences(germlines_aligned, g3_aligned, common_cols)
    return AuditResult(germlines_aligned, g3_aligned, common_cols, differences)


def print_sample_differences(audit: AuditResult, max_cols: int = 5, max_samples: int = 3):
    if not audit.differences:
        print("\nNo differences found - backends produce identical results")
        return

    print(f"\nFound differences in {len(audit.differences)} columns:")
    for col, count in sorted(audit.differences.items(), key=lambda x: -x[1]):
        print(f"  {col}: {count} differences")

    print("\n=== Sample Differences ===")
    for col in list(audit.differences.keys())[:max_cols]:
        print(f"\n--- {col} ---")
        g_vals = normalize_series(audit.germlines_df[col])
        g3_vals = normalize_series(audit.g3_df[col])
        mask = g_vals != g3_vals
        sample_idx = mask[mask].index[:max_samples]
        for idx in sample_idx:
            seq_id = audit.germlines_df.loc[idx, "sequence_id"]
            print(f"  Sequence: {seq_id}")
            print(f"    Germlines: {audit.germlines_df.loc[idx, col]}")
            print(f"    G3:        {audit.g3_df.loc[idx, col]}")


def print_summary(audit: AuditResult, num_sequences: int):
    c_diffs, stat_diffs, allele_diffs, struct_diffs = audit.get_diffs_by_category()

    print("\n" + "=" * 60)
    print("AUDIT SUMMARY (IMGT-ONLY MODE)")
    print("=" * 60)
    print(f"Sequences tested: {num_sequences}")
    print(f"Total common columns: {len(audit.common_cols)}")
    print(f"  - C region columns: {len(audit.c_region_cols)}")
    print(f"  - Statistical columns (E-values): {len([c for c in audit.common_cols if c in STATISTICAL_COLS])}")
    print(f"  - Allele-dependent columns: {len([c for c in audit.common_cols if c in ALLELE_DEPENDENT_COLS])}")
    print(f"  - Pure structural columns: {len(audit.structural_cols)}")
    print()
    print(f"OVERALL PARITY: {audit.parity_pct:.2f}%")
    print(f"  - Total values compared: {audit.total_values}")
    print(f"  - Values with differences: {audit.total_diffs}")
    print()
    print(f"PURE STRUCTURAL PARITY: {audit.structural_parity:.2f}%")
    print(f"  - Structural values compared: {len(audit.structural_cols) * len(audit.germlines_df)}")
    print(f"  - Structural differences: {sum(struct_diffs.values())}")
    if struct_diffs:
        print("  - Columns with differences:")
        for col, count in sorted(struct_diffs.items(), key=lambda x: -x[1]):
            print(f"      {col}: {count}")
    print()
    print("ALLELE-DEPENDENT DIFFERENCES (expected, IMGT versions differ):")
    print(f"  - Count: {sum(allele_diffs.values())} ({len(allele_diffs)} columns)")
    print("STATISTICAL DIFFERENCES (E-values, expected):")
    print(f"  - Count: {sum(stat_diffs.values())} ({len(stat_diffs)} columns)")
    print("=" * 60)

    if audit.structural_parity >= 99.0:
        print("PASS: Pure structural parity is excellent (>=99%)")
    elif audit.structural_parity >= 95.0:
        print("PASS: Pure structural parity is good (>=95%)")
    else:
        print(f"NEEDS REVIEW: Pure structural parity is {audit.structural_parity:.2f}%")


def main(csv_path: Optional[Path] = None):
    base_path = Path(__file__).parent.parent
    if csv_path is None:
        csv_path = base_path / "audit" / "20260112_HCV_DB_example.csv"

    rebuild_databases(base_path)
    sequences = load_sequences(csv_path)

    result_germlines = run_backend(sequences, use_germlines=True)
    result_g3 = run_backend(sequences, use_germlines=False)

    audit = compare_results(result_germlines, result_g3)
    print_sample_differences(audit)
    print_summary(audit, len(sequences))


if __name__ == "__main__":
    main()
