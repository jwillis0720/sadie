#!/usr/bin/env python
"""AIRR Backend Parity Audit - Validate germlines backend vs G3 backend."""

import os
import pandas as pd
from sadie.airr import Airr

# Load data
df = pd.read_csv("audit/20260112_HCV_DB_example.csv")
print(f"Loaded {len(df)} sequences")

# Extract heavy chain sequences
sequences = df[["sequence_id_heavy", "sequence_heavy"]].dropna()
sequences = sequences.rename(columns={"sequence_id_heavy": "sequence_id", "sequence_heavy": "sequence"})

# Handle duplicates by keeping first occurrence
if sequences["sequence_id"].duplicated().any():
    dup_count = sequences["sequence_id"].duplicated().sum()
    print(f"Found {dup_count} duplicate sequence IDs, keeping first occurrence")
    sequences = sequences.drop_duplicates(subset="sequence_id", keep="first")

print(f"Prepared {len(sequences)} unique sequences for annotation")

# C region columns to track separately
C_REGION_COLS = [
    "c_call", "c_score", "c_cigar", "c_support", "c_identity",
    "c_sequence_start", "c_sequence_end", "c_germline_start", "c_germline_end",
    "c_alignment_start", "c_alignment_end",
    "c_sequence_alignment", "c_sequence_alignment_aa",
    "c_germline_alignment", "c_germline_alignment_aa"
]

# Run with Germlines Backend
print("\n=== Running with Germlines Backend ===")
os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"
airr_germlines = Airr("human")
result_germlines = airr_germlines.run_dataframe(sequences, seq_id_field="sequence_id", seq_field="sequence")
if "source" in result_germlines.columns:
    result_germlines = result_germlines.drop(columns=["source"])
print(f"Germlines backend: {len(result_germlines)} results, {len(result_germlines.columns)} columns")

# Run with G3 Backend
print("\n=== Running with G3 Backend ===")
os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"
airr_g3 = Airr("human")
result_g3 = airr_g3.run_dataframe(sequences, seq_id_field="sequence_id", seq_field="sequence")
if "source" in result_g3.columns:
    result_g3 = result_g3.drop(columns=["source"])
print(f"G3 backend: {len(result_g3)} results, {len(result_g3.columns)} columns")

# Compare Results
print("\n=== Comparing Results ===")
germlines_cols = set(result_germlines.columns)
g3_cols = set(result_g3.columns)

print(f"Germlines columns: {len(germlines_cols)}")
print(f"G3 columns: {len(g3_cols)}")
print(f"Only in germlines: {germlines_cols - g3_cols}")
print(f"Only in G3: {g3_cols - germlines_cols}")

# Include all common columns (including C region)
common_cols = sorted(germlines_cols & g3_cols)
c_region_cols_present = [c for c in common_cols if c in C_REGION_COLS]
non_c_region_cols = [c for c in common_cols if c not in C_REGION_COLS]
print(f"C region columns present: {len(c_region_cols_present)}")
print(f"Non-C region columns: {len(non_c_region_cols)}")
print(f"Total common columns for comparison: {len(common_cols)}")

# Align data
result_germlines_aligned = result_germlines[common_cols].sort_values("sequence_id").reset_index(drop=True)
result_g3_aligned = result_g3[common_cols].sort_values("sequence_id").reset_index(drop=True)

# Find differences
differences = {}
for col in common_cols:
    # Convert to string first, then compare (handles nullable int types)
    g_vals = result_germlines_aligned[col].astype(str).replace("nan", "").replace("<NA>", "")
    g3_vals = result_g3_aligned[col].astype(str).replace("nan", "").replace("<NA>", "")
    mask = g_vals != g3_vals
    diff_count = mask.sum()
    if diff_count > 0:
        differences[col] = diff_count

if differences:
    print(f"\nFound differences in {len(differences)} columns:")
    for col, count in sorted(differences.items(), key=lambda x: -x[1]):
        print(f"  {col}: {count} differences")
    
    # Show sample differences
    print("\n=== Sample Differences ===")
    for col in list(differences.keys())[:5]:
        print(f"\n--- {col} ---")
        g_vals = result_germlines_aligned[col].astype(str).replace("nan", "").replace("<NA>", "")
        g3_vals = result_g3_aligned[col].astype(str).replace("nan", "").replace("<NA>", "")
        mask = g_vals != g3_vals
        sample_idx = mask[mask].index[:3]
        for idx in sample_idx:
            seq_id = result_germlines_aligned.loc[idx, "sequence_id"]
            print(f"  Sequence: {seq_id}")
            print(f"    Germlines: {result_germlines_aligned.loc[idx, col]}")
            print(f"    G3:        {result_g3_aligned.loc[idx, col]}")
else:
    print("\nNo differences found - backends produce identical results")

# Summary
total_values = len(common_cols) * len(result_germlines_aligned)
total_diffs = sum(differences.values()) if differences else 0
parity_pct = (1 - total_diffs / total_values) * 100 if total_values > 0 else 100

print("\n" + "=" * 50)
print("AUDIT SUMMARY")
print("=" * 50)
print(f"Sequences tested: {len(sequences)}")
print(f"Common columns (excluding C region): {len(common_cols)}")
print(f"Total values compared: {total_values}")
print(f"Values with differences: {total_diffs}")
print(f"Parity: {parity_pct:.2f}%")
print("=" * 50)

if parity_pct == 100:
    print("PASS: Backends produce identical results")
else:
    print(f"FAIL: {len(differences)} columns have differences")
