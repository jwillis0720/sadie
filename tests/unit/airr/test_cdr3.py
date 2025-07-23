"""Unit tests for CDR3 field population in AIRR annotation."""

from pathlib import Path

import pandas as pd
import pytest

from sadie.airr import Airr, AirrTable
from tests.conftest import SadieFixture


def process_permutation(args):
    """Process a single parameter permutation for CDR3 testing."""
    i, perm, seq5_fasta_str, param_names = args
    params = dict(zip(param_names, perm))

    try:
        airr_api = Airr("human", **params)
        airr_table = airr_api.run_fasta(seq5_fasta_str)

        if not airr_table.empty:
            row = airr_table.iloc[0]
            result = {
                "test_type": f"permutation_{i}",
                **params,
                "cdr3": row.get("cdr3", None),
                "cdr3_aa": row.get("cdr3_aa", None),
                "cdr3_start": row.get("cdr3_start", None),
                "cdr3_end": row.get("cdr3_end", None),
                "productive": row.get("productive", None),
                "v_call": row.get("v_call", None),
                "d_call": row.get("d_call", None),
                "j_call": row.get("j_call", None),
                "error": None,
            }
        else:
            result = {
                "test_type": f"permutation_{i}",
                **params,
                "cdr3": None,
                "cdr3_aa": None,
                "cdr3_start": None,
                "cdr3_end": None,
                "productive": None,
                "v_call": None,
                "d_call": None,
                "j_call": None,
                "error": "Empty result",
            }
    except Exception as e:
        result = {
            "test_type": f"permutation_{i}",
            **params,
            "cdr3": None,
            "cdr3_aa": None,
            "cdr3_start": None,
            "cdr3_end": None,
            "productive": None,
            "v_call": None,
            "d_call": None,
            "j_call": None,
            "error": str(e),
        }

    return result


class TestCDR3Population:
    """Test suite to ensure CDR3 fields are properly populated during AIRR annotation."""

    @pytest.fixture
    def cdr3_known_bugs_fasta(self, fixture_setup: SadieFixture) -> Path:
        """Get the path to the CDR3 known bugs fasta file."""
        return Path("tests/data/fixtures/cdr3_known_bugs.fasta")

    def test_cdr3_field_not_none(self, cdr3_known_bugs_fasta: Path) -> None:
        """Test that CDR3 field is not None for sequences in cdr3_known_bugs.fasta."""
        # Initialize Airr with human reference (assuming these are human sequences)
        airr_api = Airr("human")

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Check that we got results
        assert not airr_table.empty, "AIRR annotation returned empty table"

        # Check that CDR3 field is not None for any sequence
        cdr3_nulls = airr_table["cdr3"].isna()
        sequences_with_null_cdr3 = airr_table[cdr3_nulls]["sequence_id"].tolist()

        assert not cdr3_nulls.any(), (
            f"CDR3 field is None for sequences: {sequences_with_null_cdr3}. "
            f"Total sequences with null CDR3: {cdr3_nulls.sum()} out of {len(airr_table)}"
        )

    def test_cdr3_aa_field_not_none(self, cdr3_known_bugs_fasta: Path) -> None:
        """Test that CDR3_aa field is not None for productive sequences."""
        # Initialize Airr with human reference
        airr_api = Airr("human")

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Filter for productive sequences (CDR3_aa should be present for these)
        productive_sequences = airr_table[airr_table["productive"] == True]

        if not productive_sequences.empty:
            cdr3_aa_nulls = productive_sequences["cdr3_aa"].isna()
            sequences_with_null_cdr3_aa = productive_sequences[cdr3_aa_nulls]["sequence_id"].tolist()

            assert not cdr3_aa_nulls.any(), (
                f"CDR3_aa field is None for productive sequences: {sequences_with_null_cdr3_aa}. "
                f"Total productive sequences with null CDR3_aa: {cdr3_aa_nulls.sum()} out of {len(productive_sequences)}"
            )

    def test_cdr3_start_end_positions(self, cdr3_known_bugs_fasta: Path) -> None:
        """Test that CDR3 start and end positions are valid when CDR3 is present."""
        # Initialize Airr with human reference
        airr_api = Airr("human")

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Check sequences that have CDR3
        sequences_with_cdr3 = airr_table[airr_table["cdr3"].notna()]

        for idx, row in sequences_with_cdr3.iterrows():
            cdr3_start = row["cdr3_start"]
            cdr3_end = row["cdr3_end"]

            # Check that start and end positions are not null
            assert pd.notna(cdr3_start), f"CDR3_start is null for sequence {row['sequence_id']}"
            assert pd.notna(cdr3_end), f"CDR3_end is null for sequence {row['sequence_id']}"

            # Check that positions are valid integers
            assert isinstance(cdr3_start, (int, float)), f"CDR3_start is not numeric for sequence {row['sequence_id']}"
            assert isinstance(cdr3_end, (int, float)), f"CDR3_end is not numeric for sequence {row['sequence_id']}"

            # Check that end > start
            assert cdr3_end > cdr3_start, f"CDR3_end <= CDR3_start for sequence {row['sequence_id']}"

            # Check that the CDR3 sequence length matches the positions
            expected_length = int(cdr3_end - cdr3_start + 1)
            actual_length = len(row["cdr3"])
            assert actual_length == expected_length, (
                f"CDR3 length mismatch for sequence {row['sequence_id']}: "
                f"expected {expected_length} based on positions, got {actual_length}"
            )

    def test_cdr3_extraction_from_sequence(self, cdr3_known_bugs_fasta: Path) -> None:
        """Test that CDR3 sequence is correctly extracted from the full sequence."""
        # Initialize Airr with human reference
        airr_api = Airr("human")

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Check sequences that have CDR3
        sequences_with_cdr3 = airr_table[airr_table["cdr3"].notna()]

        for idx, row in sequences_with_cdr3.iterrows():
            if pd.notna(row["cdr3_start"]) and pd.notna(row["cdr3_end"]):
                # Extract CDR3 from the full sequence using the positions
                # Note: positions are 1-based in AIRR format
                start = int(row["cdr3_start"]) - 1  # Convert to 0-based
                end = int(row["cdr3_end"])  # End is inclusive in AIRR

                extracted_cdr3 = row["sequence"][start:end]

                # The CDR3 field might have gaps, so remove them for comparison
                cdr3_no_gaps = row["cdr3"].replace("-", "")

                assert extracted_cdr3 == cdr3_no_gaps, (
                    f"CDR3 extraction mismatch for sequence {row['sequence_id']}: "
                    f"extracted '{extracted_cdr3}' but CDR3 field has '{cdr3_no_gaps}'"
                )

    def test_all_sequences_processed(self, cdr3_known_bugs_fasta: Path) -> None:
        """Test that all sequences in the input file are processed."""
        from Bio import SeqIO

        # Count sequences in input file
        input_sequences = list(SeqIO.parse(str(cdr3_known_bugs_fasta), "fasta"))
        input_count = len(input_sequences)

        # Initialize Airr with human reference
        airr_api = Airr("human")

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Check that all sequences were processed
        assert len(airr_table) == input_count, (
            f"Not all sequences were processed: " f"input had {input_count} sequences, but output has {len(airr_table)}"
        )

        # Check that all sequence IDs match
        input_ids = {seq.id for seq in input_sequences}
        output_ids = set(airr_table["sequence_id"])

        assert input_ids == output_ids, (
            f"Sequence ID mismatch: "
            f"missing from output: {input_ids - output_ids}, "
            f"extra in output: {output_ids - input_ids}"
        )

    def test_cdr3_with_debug_enabled(self, cdr3_known_bugs_fasta: Path, caplog) -> None:
        """Test CDR3 annotation with debug mode enabled to show IgBLAST command."""
        import logging

        # Set logger to INFO level to capture debug output
        logging.getLogger("IgBLAST").setLevel(logging.INFO)

        # Initialize Airr with human reference and debug enabled
        airr_api = Airr("human", debug=True)

        # Run AIRR annotation on the test file
        airr_table = airr_api.run_fasta(str(cdr3_known_bugs_fasta))

        # Check that we got results
        assert not airr_table.empty, "AIRR annotation returned empty table"

        # Check that the IgBLAST command was logged
        assert "IgBLAST command:" in caplog.text, "IgBLAST command was not logged"
        assert "IGDATA environment variable:" in caplog.text, "IGDATA environment variable was not logged"

        # Print the captured log for visibility
        print("\n--- IgBLAST Debug Output ---")
        for record in caplog.records:
            if record.name == "IgBLAST" and record.levelname == "INFO":
                print(record.message)
        print("--- End Debug Output ---\n")

        # Still check that CDR3 is populated correctly
        cdr3_nulls = airr_table["cdr3"].isna()
        assert not cdr3_nulls.any(), f"CDR3 field is None for some sequences even with debug enabled"

    def test_parameter_permutations_seq5(self, tmp_path: Path) -> None:
        """Test Airr annotation with different parameter permutations for seq5."""
        import csv
        from functools import partial
        from itertools import product
        from multiprocessing import Pool, cpu_count

        from Bio import SeqIO

        # Extract seq5 from the fasta file
        fasta_path = Path("tests/data/fixtures/cdr3_known_bugs.fasta")
        sequences = list(SeqIO.parse(str(fasta_path), "fasta"))
        seq5 = [seq for seq in sequences if seq.id == "Seq5"][0]

        # Create a temporary fasta file with only seq5
        seq5_fasta = tmp_path / "seq5.fasta"
        SeqIO.write([seq5], str(seq5_fasta), "fasta")

        # Define specific parameter ranges optimized for antibody VDJ analysis
        # Reduced ranges for faster testing while still covering key variations
        param_ranges = {
            # Gene penalties: More negative = stricter matching
            "v_gene_penalty": [-2, -1],  # Default -1
            "d_gene_penalty": [-3, -2],  # Default -2
            "j_gene_penalty": [-3, -2],  # Default -3
            # Number of alignments to report
            "num_alignments_v": [2, 3],  # Default 2
            "num_alignments_d": [2, 3],  # Default 2
            "num_alignments_j": [3, 4],  # Default 3
            # D gene matching parameters
            "min_d_match": [6, 7, 8],  # Default 7
            # Word size for initial seed finding
            "word_size": [4, 5, 7],  # Default 5
            # Gap penalties
            "gap_open": [3, 4],  # Default 4
            "gap_extend": [1],  # Default 1
        }

        # Default parameters
        default_params = {
            "v_gene_penalty": -1,
            "d_gene_penalty": -2,
            "j_gene_penalty": -3,
            "num_alignments_v": 2,
            "num_alignments_d": 2,
            "num_alignments_j": 3,
            "min_d_match": 7,
            "word_size": 5,
            "gap_open": 4,
            "gap_extend": 1,
        }

        # Prepare results list
        results = []

        # Test 1: Default parameters
        print("\nTesting default parameters...")
        try:
            airr_api = Airr("human", **default_params)
            airr_table = airr_api.run_fasta(str(seq5_fasta))

            if not airr_table.empty:
                row = airr_table.iloc[0]
                result = {
                    "test_type": "default",
                    **default_params,
                    "cdr3": row.get("cdr3", None),
                    "cdr3_aa": row.get("cdr3_aa", None),
                    "cdr3_start": row.get("cdr3_start", None),
                    "cdr3_end": row.get("cdr3_end", None),
                    "productive": row.get("productive", None),
                    "v_call": row.get("v_call", None),
                    "d_call": row.get("d_call", None),
                    "j_call": row.get("j_call", None),
                    "error": None,
                }
            else:
                result = {
                    "test_type": "default",
                    **default_params,
                    "cdr3": None,
                    "cdr3_aa": None,
                    "cdr3_start": None,
                    "cdr3_end": None,
                    "productive": None,
                    "v_call": None,
                    "d_call": None,
                    "j_call": None,
                    "error": "Empty result",
                }
        except Exception as e:
            result = {
                "test_type": "default",
                **default_params,
                "cdr3": None,
                "cdr3_aa": None,
                "cdr3_start": None,
                "cdr3_end": None,
                "productive": None,
                "v_call": None,
                "d_call": None,
                "j_call": None,
                "error": str(e),
            }
        results.append(result)

        # Test 2: All permutations using parallel processing
        print("Testing all parameter permutations...")

        # Count total permutations
        total_perms = 1
        for param_name, values in param_ranges.items():
            total_perms *= len(values)

        print(f"Total permutations to test: {total_perms:,}")

        # Create all permutations list
        all_permutations = list(
            product(
                param_ranges["v_gene_penalty"],
                param_ranges["d_gene_penalty"],
                param_ranges["j_gene_penalty"],
                param_ranges["num_alignments_v"],
                param_ranges["num_alignments_d"],
                param_ranges["num_alignments_j"],
                param_ranges["min_d_match"],
                param_ranges["word_size"],
                param_ranges["gap_open"],
                param_ranges["gap_extend"],
            )
        )

        # Limit testing to a reasonable number of permutations
        MAX_PERMUTATIONS = 2000  # Limit for test runtime
        if len(all_permutations) > MAX_PERMUTATIONS:
            print(f"Limiting test to {MAX_PERMUTATIONS} out of {len(all_permutations)} permutations")
            import random

            random.seed(42)  # For reproducibility
            all_permutations = random.sample(all_permutations, MAX_PERMUTATIONS)

        param_names = [
            "v_gene_penalty",
            "d_gene_penalty",
            "j_gene_penalty",
            "num_alignments_v",
            "num_alignments_d",
            "num_alignments_j",
            "min_d_match",
            "word_size",
            "gap_open",
            "gap_extend",
        ]

        # Prepare arguments for parallel processing
        seq5_fasta_str = str(seq5_fasta)
        args_list = [(i, perm, seq5_fasta_str, param_names) for i, perm in enumerate(all_permutations)]

        # Use multiprocessing to run permutations in parallel
        # Limit to available CPUs or 8, whichever is smaller
        n_workers = cpu_count() * 2
        print(f"Using {n_workers} parallel workers...")

        with Pool(processes=n_workers) as pool:
            # Process in chunks to show progress
            chunk_size = 5000
            all_results = []

            for i in range(0, len(args_list), chunk_size):
                chunk = args_list[i : i + chunk_size]
                chunk_results = pool.map(process_permutation, chunk)
                all_results.extend(chunk_results)
                print(f"Progress: {min(i+chunk_size, len(args_list))}/{total_perms} permutations tested...")

            results.extend(all_results)

        # Sort results by quality metrics
        # Priority: 1) Has CDR3, 2) Productive, 3) Has J gene, 4) CDR3 length, 5) J gene specificity
        def sort_key(r):
            # Scoring system (higher is better)
            score = 0

            # Major criteria
            if r["error"] is None:
                score += 10000  # No error
            if r["cdr3"] is not None:
                score += 1000  # Has CDR3
            if r["productive"] == True:
                score += 500  # Is productive
            if r["j_call"] is not None and r["j_call"] != "":
                score += 100  # Has J gene assignment
                # Bonus for specific J gene calls (not ambiguous)
                if "," not in str(r["j_call"]):  # Single J gene assignment
                    score += 50

            # CDR3 quality metrics
            if r["cdr3"] is not None:
                cdr3_len = len(str(r["cdr3"]))
                # Typical CDR3 lengths are 24-66 nt (8-22 aa)
                if 24 <= cdr3_len <= 66:
                    score += 20
                if 30 <= cdr3_len <= 48:  # Most common range
                    score += 10

            # V and D gene presence
            if r["v_call"] is not None and r["v_call"] != "":
                score += 10
            if r["d_call"] is not None and r["d_call"] != "":
                score += 5

            return score

        # Sort results (highest score first)
        sorted_results = sorted(results, key=sort_key, reverse=True)

        # Save sorted results to CSV
        output_csv = tmp_path / "parameter_permutation_results_sorted.csv"
        fieldnames = (
            ["test_type"]
            + param_names
            + ["cdr3", "cdr3_aa", "cdr3_start", "cdr3_end", "productive", "v_call", "d_call", "j_call", "error"]
        )

        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sorted_results)

        print(f"\nResults saved to: {output_csv}")
        print(f"Total tests run: {len(results)}")

        # Also save top 100 results separately for easy analysis
        top_results_csv = tmp_path / "parameter_permutation_top100.csv"
        with open(top_results_csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sorted_results[:100])

        print(f"Top 100 results saved to: {top_results_csv}")

        # Print summary of best results
        print("\n--- Top 5 Parameter Combinations ---")
        for i, result in enumerate(sorted_results[:5]):
            print(f"\n{i+1}. Test: {result['test_type']}")
            print(f"   CDR3: {result['cdr3']}")
            print(f"   J gene: {result['j_call']}")
            print(f"   Productive: {result['productive']}")
            print(
                f"   Parameters: v_penalty={result['v_gene_penalty']}, d_penalty={result['d_gene_penalty']}, "
                f"j_penalty={result['j_gene_penalty']}, word_size={result['word_size']}"
            )

        # Verify the CSV was created
        assert output_csv.exists(), "Output CSV file was not created"

        # Read and verify the CSV content
        with open(output_csv, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            rows = list(reader)
            assert len(rows) == len(results), "CSV row count doesn't match results"

        # Check that at least some permutations produced valid CDR3
        # Need to check for both None and NaN values
        import pandas as pd

        valid_cdr3_count = sum(
            1 for r in results if r["cdr3"] is not None and pd.notna(r["cdr3"]) and r["error"] is None
        )
        print(f"\nValid CDR3 results: {valid_cdr3_count} out of {len(results)} tested permutations")

        # Analyze J gene assignment success
        j_gene_count = sum(1 for r in results if r["j_call"] is not None and r["j_call"] != "" and r["error"] is None)
        print(f"Valid J gene assignments: {j_gene_count} out of {len(results)} tested permutations")

        assert valid_cdr3_count > 0, "No valid CDR3 results were produced"
