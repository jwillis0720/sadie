#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

SOURCES_DIR = Path(__file__).parent.parent / "sources"

VALID_NUCLEOTIDES = set("ACGTNacgtn.")
IUPAC_AMBIGUOUS = set("RYSWKMBDHVryswkmbdhv")


def validate_fasta_header(header: str) -> tuple[bool, str]:
    if not header.startswith(">"):
        return False, "Header must start with '>'"
    if len(header) < 2:
        return False, "Header too short"
    return True, ""


def validate_sequence(seq: str) -> tuple[bool, str]:
    if not seq:
        return False, "Empty sequence"
    valid_chars = VALID_NUCLEOTIDES | IUPAC_AMBIGUOUS
    invalid = set(seq) - valid_chars
    if invalid:
        return False, f"Invalid characters: {invalid}"
    return True, ""


def validate_fasta_file(path: Path) -> tuple[bool, list[str]]:
    errors = []
    if not path.exists():
        return False, [f"File not found: {path}"]
    if path.stat().st_size == 0:
        return False, [f"Empty file: {path}"]

    try:
        content = path.read_text()
    except Exception as e:
        return False, [f"Read error: {e}"]

    lines = content.strip().split('\n')
    if not lines:
        return False, ["No content"]

    i = 0
    seq_count = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        header = lines[i].strip()
        valid, msg = validate_fasta_header(header)
        if not valid:
            errors.append(f"Line {i+1}: {msg}")
            i += 1
            continue

        seq_lines = []
        i += 1
        while i < len(lines) and not lines[i].startswith(">"):
            seq_lines.append(lines[i].strip())
            i += 1

        seq = "".join(seq_lines)
        valid, msg = validate_sequence(seq)
        if not valid:
            errors.append(f"Sequence for {header[:50]}: {msg}")
        seq_count += 1

    if seq_count == 0:
        errors.append("No sequences found")

    return len(errors) == 0, errors


def validate_structure(species_dir: Path) -> tuple[bool, list[str]]:
    errors = []
    if not species_dir.exists():
        return False, [f"Species directory not found: {species_dir}"]
    if not species_dir.is_dir():
        return False, [f"Not a directory: {species_dir}"]

    required_heavy = ["IGHV.fasta", "IGHD.fasta", "IGHJ.fasta"]
    for req in required_heavy:
        fpath = species_dir / req
        gapped = species_dir / req.replace(".fasta", "_gapped.fasta")
        if not fpath.exists() and not gapped.exists():
            errors.append(f"Missing required file: {req}")

    for fasta in species_dir.glob("*.fasta"):
        if not fasta.stat().st_mode & 0o444:
            errors.append(f"File not readable: {fasta.name}")

    return len(errors) == 0, errors


def validate_provider(provider: str, species: Optional[str] = None) -> tuple[bool, dict]:
    provider_dir = SOURCES_DIR / provider
    if not provider_dir.exists():
        return False, {"error": f"Provider directory not found: {provider_dir}"}

    results = {"provider": provider, "species": {}, "total_files": 0, "valid_files": 0}
    species_dirs = [provider_dir / species] if species else list(provider_dir.iterdir())

    for sp_dir in species_dirs:
        if not sp_dir.is_dir() or sp_dir.name.startswith("."):
            continue
        sp_name = sp_dir.name
        results["species"][sp_name] = {"files": {}, "valid": True}

        for fasta in sp_dir.glob("*.fasta"):
            results["total_files"] += 1
            valid, errors = validate_fasta_file(fasta)
            results["species"][sp_name]["files"][fasta.name] = {
                "valid": valid,
                "errors": errors
            }
            if valid:
                results["valid_files"] += 1
            else:
                results["species"][sp_name]["valid"] = False

    all_valid = results["valid_files"] == results["total_files"]
    return all_valid, results


def validate_all() -> tuple[bool, dict]:
    all_results = {"providers": {}, "summary": {"total": 0, "valid": 0}}
    all_valid = True

    for provider_dir in SOURCES_DIR.iterdir():
        if not provider_dir.is_dir() or provider_dir.name.startswith("."):
            continue
        provider = provider_dir.name
        valid, results = validate_provider(provider)
        all_results["providers"][provider] = results
        all_results["summary"]["total"] += results.get("total_files", 0)
        all_results["summary"]["valid"] += results.get("valid_files", 0)
        if not valid:
            all_valid = False

    return all_valid, all_results


def print_results(results: dict, verbose: bool = False):
    for provider, data in results.get("providers", {}).items():
        for sp_name, sp_data in data.get("species", {}).items():
            status = "✓" if sp_data["valid"] else "✗"
            print(f"{status} {provider}/{sp_name}")
            if verbose or not sp_data["valid"]:
                for fname, fdata in sp_data.get("files", {}).items():
                    fstatus = "✓" if fdata["valid"] else "✗"
                    print(f"    {fstatus} {fname}")
                    for err in fdata.get("errors", []):
                        print(f"        - {err}")

    summary = results.get("summary", {})
    print(f"\nSummary: {summary.get('valid', 0)}/{summary.get('total', 0)} files valid")


def main():
    parser = argparse.ArgumentParser(description="Validate germline FASTA files")
    parser.add_argument("--provider", help="Validate specific provider")
    parser.add_argument("--species", help="Validate specific species")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--file", help="Validate specific file")
    args = parser.parse_args()

    if args.file:
        path = Path(args.file)
        valid, errors = validate_fasta_file(path)
        if valid:
            print(f"✓ {path.name} is valid")
            sys.exit(0)
        else:
            print(f"✗ {path.name} has errors:")
            for err in errors:
                print(f"  - {err}")
            sys.exit(1)

    if args.provider:
        valid, results = validate_provider(args.provider, args.species)
        print_results({"providers": {args.provider: results}}, args.verbose)
    else:
        valid, results = validate_all()
        print_results(results, args.verbose)

    sys.exit(0 if valid else 1)


if __name__ == "__main__":
    main()
