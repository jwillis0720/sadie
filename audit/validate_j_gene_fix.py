#!/usr/bin/env python
"""
Validate J Gene Fix
===================

Tests that J gene matching and CDR3 annotation work after aux file fix.
"""

import os
from pathlib import Path

import pandas as pd

# Force germlines backend
os.environ["SADIE_USE_GERMLINES_MODULE"] = "true"

from sadie.airr import Airr


def main():
    # Test sequence (known productive heavy chain)
    test_fasta = Path("audit/test_sequences.fasta")

    # Create test FASTA if not exists
    if not test_fasta.exists():
        test_fasta.write_text(
            ">test_seq1\n"
            "CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTC"
            "TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGGCATGCACTGGGTCCGCCAGGCT"
            "CCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTAC"
            "GCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTAT"
            "CTGCAAATGAACAGCCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGA"
            "CGGTTTGCTTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG\n"
        )

    # Run annotation with germlines backend
    airr = Airr("human")

    # Read the sequence
    with open(test_fasta) as f:
        lines = f.readlines()
        seq_id = lines[0].strip().lstrip(">")
        seq = "".join(l.strip() for l in lines[1:])

    # Create a single-row dataframe
    df = pd.DataFrame([{"sequence_id": seq_id, "sequence": seq}])
    result_df = airr.run_dataframe(df, seq_id_field="sequence_id", seq_field="sequence")

    # Get the first (only) result
    result = result_df.iloc[0].to_dict()

    # Check critical fields
    print("=== J Gene Fix Validation ===")
    print(f"v_call: {result.get('v_call', 'N/A')}")
    print(f"d_call: {result.get('d_call', 'N/A')}")
    print(f"j_call: {result.get('j_call', 'N/A')}")  # Should NOT be NaN
    print(f"junction: {result.get('junction', 'N/A')}")  # Should NOT be NaN
    print(f"junction_aa: {result.get('junction_aa', 'N/A')}")
    print(f"cdr3: {result.get('cdr3', 'N/A')}")  # Should NOT be NaN
    print(f"cdr3_aa: {result.get('cdr3_aa', 'N/A')}")
    print(f"fwr4: {result.get('fwr4', 'N/A')}")  # Should NOT be NaN
    print(f"complete_vdj: {result.get('complete_vdj', 'N/A')}")  # Should be True

    # Validation checks
    success = True
    if pd.isna(result.get('j_call')) or result.get('j_call') is None:
        print("\n❌ FAIL: j_call is NaN/None")
        success = False
    else:
        print("\n✓ PASS: j_call is populated")

    if pd.isna(result.get('junction')) or result.get('junction') is None:
        print("❌ FAIL: junction is NaN/None")
        success = False
    else:
        print("✓ PASS: junction is populated")

    if pd.isna(result.get('cdr3')) or result.get('cdr3') is None:
        print("❌ FAIL: cdr3 is NaN/None")
        success = False
    else:
        print("✓ PASS: cdr3 is populated")

    if pd.isna(result.get('fwr4')) or result.get('fwr4') is None:
        print("❌ FAIL: fwr4 is NaN/None")
        success = False
    else:
        print("✓ PASS: fwr4 is populated")

    print(f"\n=== Overall: {'SUCCESS' if success else 'FAILURE'} ===")
    return success


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
