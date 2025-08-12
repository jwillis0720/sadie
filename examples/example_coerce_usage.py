#!/usr/bin/env python
"""Example of using the coerce parameter in SADIE AIRR annotation"""

from pathlib import Path

from sadie.airr import Airr

# Example FASTA file with sequences that may have alleles not in auxiliary files
fasta_file = Path("tests/data/fixtures/cdr3_known_bugs.fasta")

print("Running AIRR annotation without coercion...")
# Standard annotation - CDR3 may be null if exact allele match not found in auxiliary files
airr_standard = Airr("human", coerce=False)
result_standard = airr_standard.run_fasta(fasta_file)

# Check for sequences with missing CDR3
missing_cdr3 = result_standard['cdr3'].isna().sum()
print(f"Sequences with missing CDR3 (no coerce): {missing_cdr3}")

print("\nRunning AIRR annotation with coercion...")
# Annotation with coercion - will use closest available allele for CDR3 boundaries
airr_coerce = Airr("human", coerce=True)
result_coerce = airr_coerce.run_fasta(fasta_file)

# Check for sequences with missing CDR3 after coercion
missing_cdr3_coerce = result_coerce['cdr3'].isna().sum()
print(f"Sequences with missing CDR3 (with coerce): {missing_cdr3_coerce}")

# Show which sequences were coerced
coerced_mask = result_coerce['comments'].str.contains('FR/CDR boundaries derived from', na=False)
if coerced_mask.any():
    print("\nCoerced sequences:")
    for idx, row in result_coerce[coerced_mask].iterrows():
        print(f"  {row['sequence_id']}: {row['comments']}")
else:
    print("\nNo sequences required coercion.")

# Example of the improved annotation
print("\nExample comparison for IGKV3-15*02 sequences:")
igkv3_15_mask = result_standard['v_call'].str.contains('IGKV3-15', na=False)
if igkv3_15_mask.any():
    for idx in result_standard[igkv3_15_mask].index:
        print(f"\nSequence {result_standard.loc[idx, 'sequence_id']}:")
        print(f"  Without coerce - V call: {result_standard.loc[idx, 'v_call']}, CDR3: {result_standard.loc[idx, 'cdr3']}")
        print(f"  With coerce - V call: {result_coerce.loc[idx, 'v_call']}, CDR3: {result_coerce.loc[idx, 'cdr3']}")
