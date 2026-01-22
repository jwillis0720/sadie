#!/usr/bin/env python3
"""
Download OGRDB Germline Data
============================

Script to download OGRDB germline sequences from the Zenodo archive.

Data Source: https://zenodo.org/records/18145568
Archive URL: https://zenodo.org/records/18145568/files/ogrdb_archive.tgz?download=1

The OGRDB archive contains:
- SQL dump with gene_description table containing:
  - sequence: ungapped nucleotide sequence
  - coding_seq_imgt: IMGT-gapped nucleotide sequence
- Genotype CSV files with additional metadata

Usage:
    python download_ogrdb.py --species human
    python download_ogrdb.py --species human mouse --output-dir ./data
    python download_ogrdb.py --help

Output Directory Structure:
    sources/ogrdb/
    ├── human/
    │   ├── IGHV.fasta          # Ungapped sequences
    │   ├── IGHV_gapped.fasta   # IMGT-gapped sequences
    │   ├── IGHD.fasta
    │   ├── IGHJ.fasta
    │   ├── IGHJ_gapped.fasta
    │   └── ...
    └── mouse/
        └── ...
"""

import argparse
import json
import logging
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.request import urlretrieve
from urllib.error import URLError

logger = logging.getLogger(__name__)

CHECKPOINT_FILE = ".download_progress.json"


def _load_checkpoint(checkpoint_path: Path) -> dict:
    if checkpoint_path.exists():
        try:
            return json.loads(checkpoint_path.read_text())
        except Exception:
            return {}
    return {}


def _save_checkpoint(checkpoint_path: Path, data: dict):
    checkpoint_path.write_text(json.dumps(data, indent=2))

# Zenodo archive URL
ZENODO_ARCHIVE_URL = "https://zenodo.org/records/18145568/files/ogrdb_archive.tgz?download=1"

# Species name mapping (OGRDB names to internal names)
SPECIES_MAP = {
    "Homo sapiens": "human",
    "Mus musculus": "mouse",
    "Macaca mulatta": "rhesus_macaque",
    "Rattus norvegicus": "rat",
    "Canis lupus familiaris": "dog",
    "Oryctolagus cuniculus": "rabbit",
    # Add more as needed
}

SPECIES_MAP_REVERSE = {v: k for k, v in SPECIES_MAP.items()}

# Segment type detection from gene names
SEGMENT_PATTERNS = {
    "V": re.compile(r"IG[HKL]V|TR[ABGD]V", re.IGNORECASE),
    "D": re.compile(r"IG[H]D|TR[BD]D", re.IGNORECASE),
    "J": re.compile(r"IG[HKL]J|TR[ABGD]J", re.IGNORECASE),
    "C": re.compile(r"IG[HKL]C|TR[ABGD]C", re.IGNORECASE),
}

# Chain detection from gene names
CHAIN_PATTERNS = {
    "H": re.compile(r"IGH[VDJC]", re.IGNORECASE),
    "K": re.compile(r"IGK[VJC]", re.IGNORECASE),
    "L": re.compile(r"IGL[VJC]", re.IGNORECASE),
}


class OGRDBDownloader:
    """Download and process OGRDB archive from Zenodo."""
    
    def __init__(
        self,
        output_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None
    ):
        """
        Initialize OGRDB downloader.
        
        Parameters
        ----------
        output_dir : Path, optional
            Output directory for FASTA files.
            Defaults to sources/ogrdb/
        cache_dir : Path, optional
            Cache directory for downloaded archive.
            Defaults to ~/.cache/sadie/ogrdb/
        """
        if output_dir is None:
            output_dir = Path(__file__).parent.parent / "sources" / "ogrdb"
        
        if cache_dir is None:
            cache_dir = Path.home() / ".cache" / "sadie" / "ogrdb"
        
        self.output_dir = Path(output_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def download(
        self,
        species: List[str],
        force: bool = False
    ) -> None:
        """
        Download OGRDB data for specified species.
        
        Parameters
        ----------
        species : List[str]
            Species to extract (e.g., ["human", "mouse"])
        force : bool
            Force re-download even if cached
        """
        start_time = time.time()
        
        # Download archive
        archive_path = self._download_archive(force=force)
        
        # Extract and process archive
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Extract archive
            logger.info(f"Extracting archive to {temp_path}...")
            self._extract_archive(archive_path, temp_path)
            
            # Find and process SQL dump
            sql_files = list(temp_path.rglob("*.sql"))
            if not sql_files:
                # Try looking for SQLite database
                db_files = list(temp_path.rglob("*.db")) + list(temp_path.rglob("*.sqlite"))
                if db_files:
                    for db_file in db_files:
                        self._process_sqlite_db(db_file, species)
                else:
                    raise FileNotFoundError(
                        "No SQL dump or SQLite database found in archive. "
                        f"Contents: {list(temp_path.rglob('*'))}"
                    )
            else:
                for sql_file in sql_files:
                    logger.info(f"Processing SQL file: {sql_file}")
                    self._process_sql_dump(sql_file, species)
        
        duration_ms = int((time.time() - start_time) * 1000)
        logger.info(
            f"operation=download provider=ogrdb "
            f"species={','.join(species)} duration_ms={duration_ms} status=success"
        )
    
    def _download_archive(self, force: bool = False) -> Path:
        """
        Download OGRDB archive from Zenodo.
        
        Parameters
        ----------
        force : bool
            Force re-download
            
        Returns
        -------
        Path
            Path to downloaded archive
        """
        archive_path = self.cache_dir / "ogrdb_archive.tgz"
        
        if archive_path.exists() and not force:
            logger.info(f"Using cached archive: {archive_path}")
            return archive_path
        
        logger.info(f"Downloading OGRDB archive from Zenodo...")
        logger.info(f"URL: {ZENODO_ARCHIVE_URL}")
        
        try:
            def progress_hook(count, block_size, total_size):
                if total_size > 0:
                    percent = int(count * block_size * 100 / total_size)
                    if count % 100 == 0:  # Log every 100 blocks
                        logger.info(f"Download progress: {percent}%")
            
            urlretrieve(ZENODO_ARCHIVE_URL, archive_path, reporthook=progress_hook)
            logger.info(f"Downloaded archive to {archive_path}")
            
        except URLError as e:
            raise RuntimeError(
                f"Failed to download OGRDB archive: {e}. "
                "Check your internet connection or download manually from "
                f"{ZENODO_ARCHIVE_URL}"
            )
        
        return archive_path
    
    def _extract_archive(self, archive_path: Path, extract_dir: Path) -> None:
        """
        Extract tgz archive.
        
        Parameters
        ----------
        archive_path : Path
            Path to archive
        extract_dir : Path
            Directory to extract to
        """
        try:
            with tarfile.open(archive_path, "r:gz") as tar:
                tar.extractall(extract_dir)
        except tarfile.TarError as e:
            raise RuntimeError(f"Failed to extract archive: {e}")
    
    def _process_sql_dump(
        self,
        sql_path: Path,
        species: List[str]
    ) -> None:
        """
        Process SQL dump to extract sequences from gene_description table.
        
        The gene_description table has columns:
        - id (0), sequence_name (1), description_id (2), lab_address (3)
        - release_version (4), release_date (5), release_description (6)
        - alt_names (7), locus (8), sequence_type (9), inference_type (10)
        - affirmation_level (11), status (12), gene_subgroup (13)
        - subgroup_designation (14), allele_designation (15)
        - sequence (16), coding_seq_imgt (17), ... species is at index ~42
        
        Parameters
        ----------
        sql_path : Path
            Path to SQL dump file
        species : List[str]
            Species to extract
        """
        logger.info(f"Parsing SQL dump: {sql_path}")
        
        # Convert species to OGRDB format for matching
        species_set = set()
        for sp in species:
            ogrdb_name = SPECIES_MAP_REVERSE.get(sp, sp)
            species_set.add(ogrdb_name.lower())
            species_set.add(sp.lower())
        
        # Sequences organized by species/segment/chain
        sequences: Dict[str, Dict[str, Dict[str, List[Tuple[str, str, str]]]]] = {}
        
        # Read SQL file and find gene_description INSERT statements
        with open(sql_path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        
        # Find ALL gene_description INSERT blocks
        # Format: INSERT INTO `gene_description` VALUES\n(row1),\n(row2),...;
        gene_desc_matches = re.findall(
            r"INSERT INTO [`'\"]?gene_description[`'\"]?\s+VALUES\s*\n?(.*?);",
            content,
            re.IGNORECASE | re.DOTALL
        )
        
        if not gene_desc_matches:
            logger.warning("No gene_description INSERT statement found")
            return
        
        logger.info(f"Found {len(gene_desc_matches)} INSERT blocks for gene_description")
        
        # Parse all INSERT blocks
        all_rows = []
        for values_block in gene_desc_matches:
            rows = self._parse_multi_row_insert(values_block)
            all_rows.extend(rows)
        
        logger.info(f"Found {len(all_rows)} total rows in gene_description table")
        
        count = 0
        for row in all_rows:
            result = self._parse_gene_description_row(row, species_set)
            if result:
                gene_name, sp, locus, seq_type, ungapped, gapped = result
                self._add_sequence_to_dict(
                    sequences, gene_name, sp, locus, seq_type, ungapped, gapped
                )
                count += 1
        
        logger.info(f"Extracted {count} sequences for species: {species}")
        
        # Write FASTA files
        self._write_fasta_files(sequences)
    
    def _parse_multi_row_insert(self, values_block: str) -> List[str]:
        """
        Parse a multi-row INSERT VALUES block into individual rows.
        
        Handles: (val1, val2, ...),\n(val1, val2, ...),...
        
        Parameters
        ----------
        values_block : str
            The VALUES portion of an INSERT statement
            
        Returns
        -------
        List[str]
            List of row strings (without outer parentheses)
        """
        rows = []
        depth = 0
        current_row = []
        in_string = False
        string_char = None
        i = 0
        
        while i < len(values_block):
            char = values_block[i]
            
            # Handle string escaping
            if in_string:
                current_row.append(char)
                if char == string_char:
                    # Check for escaped quote
                    if i + 1 < len(values_block) and values_block[i + 1] == string_char:
                        current_row.append(values_block[i + 1])
                        i += 1
                    else:
                        in_string = False
                elif char == '\\' and i + 1 < len(values_block):
                    # Escape sequence
                    current_row.append(values_block[i + 1])
                    i += 1
            elif char in ("'", '"'):
                in_string = True
                string_char = char
                current_row.append(char)
            elif char == '(':
                if depth == 0:
                    current_row = []  # Start new row
                else:
                    current_row.append(char)
                depth += 1
            elif char == ')':
                depth -= 1
                if depth == 0:
                    # End of row
                    rows.append(''.join(current_row))
                else:
                    current_row.append(char)
            elif depth > 0:
                current_row.append(char)
            
            i += 1
        
        return rows
    
    def _parse_gene_description_row(
        self,
        row: str,
        species_set: set
    ) -> Optional[Tuple[str, str, str, str, str, str]]:
        """
        Parse a single gene_description row.
        
        Column indices (0-indexed):
        - 1: sequence_name (gene name like IGHV1-2*i01)
        - 8: locus (IGH, IGK, IGL)
        - 9: sequence_type (V, D, J, C)
        - 16: sequence (ungapped)
        - 17: coding_seq_imgt (gapped)
        - ~42: species
        
        Parameters
        ----------
        row : str
            Comma-separated values string
        species_set : set
            Set of species to include (lowercase)
            
        Returns
        -------
        Tuple or None
            (gene_name, species, locus, seq_type, ungapped, gapped) or None
        """
        values = self._split_sql_row(row)
        
        if len(values) < 43:
            return None
        
        # Extract fields
        sequence_name = self._clean_sql_string(values[1])
        locus = self._clean_sql_string(values[8])
        seq_type = self._clean_sql_string(values[9])
        ungapped = self._clean_sql_string(values[16])
        gapped = self._clean_sql_string(values[17])
        
        # Species is around index 42 (after many NULL fields)
        # Look for it by finding 'Homo sapiens' or similar
        species_val = None
        for i in range(40, min(50, len(values))):
            val = self._clean_sql_string(values[i])
            if val and any(sp in val.lower() for sp in species_set):
                species_val = val
                break
        
        if not species_val:
            return None
        
        # Validate we have required data
        if not sequence_name or not locus or not seq_type:
            return None
        
        if not ungapped and not gapped:
            return None
        
        # Map species to internal name
        # SPECIES_MAP: {"Homo sapiens": "human", "Mus musculus": "mouse", ...}
        species_internal = None
        for ogrdb_name, internal_name in SPECIES_MAP.items():
            if ogrdb_name.lower() in species_val.lower():
                species_internal = internal_name
                break
        
        if not species_internal:
            # Try matching internal name directly (e.g., "human", "mouse")
            for ogrdb_name, internal_name in SPECIES_MAP.items():
                if internal_name.lower() in species_val.lower():
                    species_internal = internal_name
                    break
        
        if not species_internal:
            return None  # Skip if species not recognized
        
        return (sequence_name, species_internal, locus, seq_type, ungapped, gapped)
    
    def _split_sql_row(self, row: str) -> List[str]:
        """
        Split a SQL row into values, handling quoted strings and NULLs.
        
        Parameters
        ----------
        row : str
            Comma-separated values
            
        Returns
        -------
        List[str]
            List of values
        """
        values = []
        current = []
        in_string = False
        string_char = None
        i = 0
        
        while i < len(row):
            char = row[i]
            
            if in_string:
                if char == string_char:
                    # Check for escaped quote
                    if i + 1 < len(row) and row[i + 1] == string_char:
                        current.append(char)
                        i += 1
                    else:
                        in_string = False
                        current.append(char)
                elif char == '\\' and i + 1 < len(row):
                    # Escape sequence - skip the backslash for common escapes
                    next_char = row[i + 1]
                    if next_char in ('n', 'r', 't'):
                        current.append({'n': '\n', 'r': '\r', 't': '\t'}[next_char])
                    else:
                        current.append(next_char)
                    i += 1
                else:
                    current.append(char)
            elif char in ("'", '"'):
                in_string = True
                string_char = char
                current.append(char)
            elif char == ',':
                values.append(''.join(current).strip())
                current = []
            else:
                current.append(char)
            
            i += 1
        
        if current:
            values.append(''.join(current).strip())
        
        return values
    
    def _clean_sql_string(self, value: str) -> Optional[str]:
        """
        Clean a SQL string value.
        
        Parameters
        ----------
        value : str
            SQL value (may be quoted or NULL)
            
        Returns
        -------
        str or None
            Cleaned string value
        """
        if not value or value.upper() == 'NULL':
            return None
        
        # Remove quotes
        if (value.startswith("'") and value.endswith("'")) or \
           (value.startswith('"') and value.endswith('"')):
            value = value[1:-1]
        
        # Handle empty strings
        if not value.strip():
            return None
        
        return value
    
    def _add_sequence_to_dict(
        self,
        sequences: Dict,
        gene_name: str,
        species: str,
        locus: str,
        seq_type: str,
        ungapped: Optional[str],
        gapped: Optional[str]
    ) -> None:
        """
        Add a sequence to the sequences dictionary.
        
        Parameters
        ----------
        sequences : Dict
            Output dictionary
        gene_name : str
            Gene name (e.g., IGHV1-2*i01)
        species : str
            Internal species name
        locus : str
            Locus (IGH, IGK, IGL)
        seq_type : str
            Sequence type (V, D, J)
        ungapped : str
            Ungapped sequence
        gapped : str
            Gapped sequence
        """
        # Determine chain from locus
        chain_map = {"IGH": "H", "IGK": "K", "IGL": "L", "TRA": "A", "TRB": "B", "TRD": "D", "TRG": "G"}
        chain = chain_map.get(locus)
        
        if not chain:
            return
        
        # Only process V, D, J segments
        segment = seq_type.upper()
        if segment not in ("V", "D", "J"):
            return
        
        # Initialize nested structure
        if species not in sequences:
            sequences[species] = {}
        if segment not in sequences[species]:
            sequences[species][segment] = {}
        if chain not in sequences[species][segment]:
            sequences[species][segment][chain] = []
        
        sequences[species][segment][chain].append((gene_name, ungapped, gapped))
    
    def _parse_insert_values(
        self,
        values_str: str,
        species_patterns: List[str],
        sequences: Dict
    ) -> None:
        """
        Parse INSERT VALUES string to extract sequence data.
        
        Parameters
        ----------
        values_str : str
            SQL VALUES string
        species_patterns : List[str]
            Species patterns to match
        sequences : Dict
            Output sequences dictionary
        """
        # Split by comma, handling quoted strings
        # This is a simplified parser - may need refinement
        values = self._split_sql_values(values_str)
        
        # Expected columns (may vary):
        # sequence_name, species, sequence, coding_seq_imgt, ...
        # We'll look for patterns
        
        gene_name = None
        species_name = None
        ungapped_seq = None
        gapped_seq = None
        
        for i, val in enumerate(values):
            val = val.strip().strip("'\"")
            
            # Detect gene name (e.g., IGHV1-69*01)
            if re.match(r"^IG[HKL][VDJC]", val, re.IGNORECASE) or \
               re.match(r"^TR[ABGD][VDJC]", val, re.IGNORECASE):
                gene_name = val
            
            # Detect species
            for sp_pattern in species_patterns:
                if sp_pattern in val.lower():
                    species_name = val
                    break
            
            # Detect sequences (long strings of ACGT)
            if len(val) > 50 and re.match(r"^[ACGTN.]+$", val, re.IGNORECASE):
                if "." in val or "-" in val:
                    gapped_seq = val.upper()
                else:
                    ungapped_seq = val.upper()
        
        if gene_name and (ungapped_seq or gapped_seq):
            self._add_sequence(
                sequences, gene_name, species_name,
                ungapped_seq, gapped_seq, species_patterns
            )
    
    def _split_sql_values(self, values_str: str) -> List[str]:
        """
        Split SQL VALUES string handling quoted strings.
        
        Parameters
        ----------
        values_str : str
            SQL VALUES string
            
        Returns
        -------
        List[str]
            Split values
        """
        values = []
        current = []
        in_quote = False
        quote_char = None
        
        for char in values_str:
            if char in ("'", '"') and not in_quote:
                in_quote = True
                quote_char = char
                current.append(char)
            elif char == quote_char and in_quote:
                in_quote = False
                current.append(char)
                quote_char = None
            elif char == "," and not in_quote:
                values.append("".join(current))
                current = []
            else:
                current.append(char)
        
        if current:
            values.append("".join(current))
        
        return values
    
    def _parse_sql_alternative(
        self,
        sql_path: Path,
        species: List[str],
        sequences: Dict
    ) -> None:
        """
        Alternative SQL parsing for different dump formats.
        
        Parameters
        ----------
        sql_path : Path
            Path to SQL file
        species : List[str]
            Species to extract
        sequences : Dict
            Output sequences dictionary
        """
        # Convert SQL dump to SQLite and query
        logger.info("Converting SQL dump to SQLite for parsing...")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = Path(temp_dir) / "ogrdb.db"
            
            try:
                # Create SQLite database from SQL dump
                self._sql_to_sqlite(sql_path, db_path)
                self._process_sqlite_db(db_path, species, sequences)
            except Exception as e:
                logger.warning(f"SQLite conversion failed: {e}")
                # Fall back to regex-based parsing
                self._parse_sql_regex(sql_path, species, sequences)
    
    def _sql_to_sqlite(self, sql_path: Path, db_path: Path) -> None:
        """
        Convert SQL dump to SQLite database.
        
        Parameters
        ----------
        sql_path : Path
            Path to SQL dump
        db_path : Path
            Path for SQLite database
        """
        # Read SQL and convert MySQL/PostgreSQL syntax to SQLite
        with open(sql_path, "r", encoding="utf-8", errors="replace") as f:
            sql_content = f.read()
        
        # Basic MySQL to SQLite conversions
        sql_content = re.sub(r"AUTO_INCREMENT", "", sql_content, flags=re.IGNORECASE)
        sql_content = re.sub(r"ENGINE\s*=\s*\w+", "", sql_content, flags=re.IGNORECASE)
        sql_content = re.sub(r"DEFAULT CHARSET\s*=\s*\w+", "", sql_content, flags=re.IGNORECASE)
        sql_content = re.sub(r"UNSIGNED", "", sql_content, flags=re.IGNORECASE)
        sql_content = re.sub(r"`", '"', sql_content)
        
        conn = sqlite3.connect(db_path)
        try:
            conn.executescript(sql_content)
            conn.commit()
        finally:
            conn.close()
    
    def _process_sqlite_db(
        self,
        db_path: Path,
        species: List[str],
        sequences: Optional[Dict] = None
    ) -> None:
        """
        Process SQLite database to extract sequences.
        
        Parameters
        ----------
        db_path : Path
            Path to SQLite database
        species : List[str]
            Species to extract
        sequences : Dict, optional
            Output sequences dictionary (created if None)
        """
        if sequences is None:
            sequences = {}
        
        species_patterns = []
        for sp in species:
            ogrdb_name = SPECIES_MAP_REVERSE.get(sp, sp)
            species_patterns.append(ogrdb_name.lower())
            species_patterns.append(sp.lower())
        
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        
        try:
            cursor = conn.cursor()
            
            # Get table names
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = [row[0] for row in cursor.fetchall()]
            logger.info(f"Found tables: {tables}")
            
            # Look for gene_description or similar table
            gene_table = None
            for table in tables:
                if "gene" in table.lower() or "sequence" in table.lower():
                    gene_table = table
                    break
            
            if not gene_table:
                logger.warning(f"No gene table found. Tables: {tables}")
                return
            
            # Get column names
            cursor.execute(f"PRAGMA table_info({gene_table})")
            columns = [row[1] for row in cursor.fetchall()]
            logger.info(f"Columns in {gene_table}: {columns}")
            
            # Build query to extract sequences
            cursor.execute(f"SELECT * FROM {gene_table}")
            
            for row in cursor.fetchall():
                row_dict = dict(row)
                self._process_db_row(row_dict, species_patterns, sequences)
        
        finally:
            conn.close()
        
        self._write_fasta_files(sequences)
    
    def _process_db_row(
        self,
        row: Dict,
        species_patterns: List[str],
        sequences: Dict
    ) -> None:
        """
        Process a database row to extract sequence.
        
        Parameters
        ----------
        row : Dict
            Database row
        species_patterns : List[str]
            Species patterns to match
        sequences : Dict
            Output sequences dictionary
        """
        # Find gene name column
        gene_name = None
        for col in ["sequence_name", "gene_name", "name", "allele_name"]:
            if col in row and row[col]:
                gene_name = row[col]
                break
        
        if not gene_name:
            return
        
        # Find species column
        species_name = None
        for col in ["species", "organism", "species_subgroup"]:
            if col in row and row[col]:
                species_name = row[col]
                break
        
        # Find ungapped sequence
        ungapped_seq = None
        for col in ["sequence", "nt_sequence", "nucleotide_sequence"]:
            if col in row and row[col]:
                seq = row[col].upper()
                if re.match(r"^[ACGTN]+$", seq):
                    ungapped_seq = seq
                    break
        
        # Find gapped sequence
        gapped_seq = None
        for col in ["coding_seq_imgt", "nt_sequence_gapped", "gapped_sequence", "imgt_sequence"]:
            if col in row and row[col]:
                seq = row[col].upper()
                if re.match(r"^[ACGTN.]+$", seq):
                    gapped_seq = seq
                    break
        
        if gene_name and (ungapped_seq or gapped_seq):
            self._add_sequence(
                sequences, gene_name, species_name,
                ungapped_seq, gapped_seq, species_patterns
            )
    
    def _parse_sql_regex(
        self,
        sql_path: Path,
        species: List[str],
        sequences: Dict
    ) -> None:
        """
        Parse SQL dump using regex patterns.
        
        Parameters
        ----------
        sql_path : Path
            Path to SQL file
        species : List[str]
            Species to extract
        sequences : Dict
            Output sequences dictionary
        """
        logger.info("Using regex-based SQL parsing...")
        
        species_patterns = []
        for sp in species:
            ogrdb_name = SPECIES_MAP_REVERSE.get(sp, sp)
            species_patterns.append(ogrdb_name.lower())
            species_patterns.append(sp.lower())
        
        with open(sql_path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        
        # Pattern to match gene sequences in various formats
        # Looking for: gene_name, species, sequence patterns
        gene_pattern = re.compile(
            r"'((?:IG|TR)[HKL]?[VDJC][^']*\*\d+[^']*)'"  # Gene name like IGHV1-69*01
            r"[^']*"
            r"'([^']{50,})'"  # Sequence (at least 50 chars)
        )
        
        for match in gene_pattern.finditer(content):
            gene_name = match.group(1)
            sequence = match.group(2).upper()
            
            if re.match(r"^[ACGTN.]+$", sequence):
                is_gapped = "." in sequence
                self._add_sequence(
                    sequences, gene_name, None,
                    None if is_gapped else sequence,
                    sequence if is_gapped else None,
                    species_patterns
                )
    
    def _add_sequence(
        self,
        sequences: Dict,
        gene_name: str,
        species_name: Optional[str],
        ungapped_seq: Optional[str],
        gapped_seq: Optional[str],
        species_patterns: List[str]
    ) -> None:
        """
        Add sequence to sequences dictionary.
        
        Parameters
        ----------
        sequences : Dict
            Output sequences dictionary
        gene_name : str
            Gene name
        species_name : str, optional
            Species name from database
        ungapped_seq : str, optional
            Ungapped sequence
        gapped_seq : str, optional
            Gapped sequence
        species_patterns : List[str]
            Allowed species patterns
        """
        # Determine species
        species = "human"  # Default
        if species_name:
            for internal, ogrdb in SPECIES_MAP_REVERSE.items():
                if species_name.lower() in ogrdb.lower() or \
                   ogrdb.lower() in species_name.lower():
                    species = internal
                    break
        
        # Check if this species is requested
        if species.lower() not in [p.lower() for p in species_patterns]:
            # Species not in requested list, still include if matches
            for pattern in species_patterns:
                if pattern in species.lower():
                    break
            else:
                return  # Skip this species
        
        # Determine segment and chain from gene name
        segment = None
        chain = None
        
        for seg, pattern in SEGMENT_PATTERNS.items():
            if pattern.search(gene_name):
                segment = seg
                break
        
        for ch, pattern in CHAIN_PATTERNS.items():
            if pattern.search(gene_name):
                chain = ch
                break
        
        if not segment or not chain:
            logger.debug(f"Could not determine segment/chain for {gene_name}")
            return
        
        # Initialize nested dict structure
        if species not in sequences:
            sequences[species] = {}
        if segment not in sequences[species]:
            sequences[species][segment] = {}
        if chain not in sequences[species][segment]:
            sequences[species][segment][chain] = []
        
        # Add sequence (gene_name, ungapped, gapped)
        sequences[species][segment][chain].append((gene_name, ungapped_seq, gapped_seq))
    
    def _write_fasta_files(self, sequences: Dict) -> None:
        """
        Write sequences to FASTA files.
        
        Parameters
        ----------
        sequences : Dict
            Sequences organized by species/segment/chain
        """
        for species, segments in sequences.items():
            species_dir = self.output_dir / species
            species_dir.mkdir(parents=True, exist_ok=True)
            
            for segment, chains in segments.items():
                for chain, seqs in chains.items():
                    if not seqs:
                        continue
                    
                    # Deduplicate by gene name
                    seen = set()
                    unique_seqs = []
                    for gene_name, ungapped, gapped in seqs:
                        if gene_name not in seen:
                            seen.add(gene_name)
                            unique_seqs.append((gene_name, ungapped, gapped))
                    
                    # Write ungapped FASTA
                    ungapped_path = species_dir / f"IG{chain}{segment}.fasta"
                    ungapped_count = 0
                    with open(ungapped_path, "w") as f:
                        for gene_name, ungapped, gapped in unique_seqs:
                            seq = ungapped or (gapped.replace(".", "") if gapped else None)
                            if seq:
                                f.write(f">{gene_name}\n")
                                f.write(f"{seq}\n")
                                ungapped_count += 1
                    
                    if ungapped_count > 0:
                        logger.info(f"Wrote {ungapped_count} sequences to {ungapped_path}")
                    
                    # Write gapped FASTA (if any gapped sequences)
                    gapped_seqs = [(n, g) for n, u, g in unique_seqs if g]
                    if gapped_seqs:
                        gapped_path = species_dir / f"IG{chain}{segment}_gapped.fasta"
                        with open(gapped_path, "w") as f:
                            for gene_name, gapped in gapped_seqs:
                                f.write(f">{gene_name}\n")
                                f.write(f"{gapped}\n")
                        logger.info(f"Wrote {len(gapped_seqs)} gapped sequences to {gapped_path}")


def main():
    """Main entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        description="Download OGRDB germline data from Zenodo archive"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=["human"],
        help="Species to download (e.g., human mouse)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for FASTA files"
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Cache directory for downloaded archive"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if cached"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    downloader = OGRDBDownloader(
        output_dir=args.output_dir,
        cache_dir=args.cache_dir
    )
    
    try:
        downloader.download(args.species, force=args.force)
        print(f"\nOGRDB data downloaded successfully to {downloader.output_dir}")
    except Exception as e:
        logger.error(f"Download failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
