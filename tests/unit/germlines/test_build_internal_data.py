"""Tests for build_internal_data script and combined VDJC structure.

These tests verify that:
1. build_internal_data creates combined VDJC FASTA files (not symlinks)
2. The BLAST database is built from the combined file
3. No symlinks exist in internal_data directories
4. GermlineData paths point to the correct locations
"""

import pytest
from pathlib import Path

from sadie.germlines import get_germlines_base_dir


class TestInternalDataStructure:
    """Tests for the internal_data directory structure."""

    @pytest.fixture
    def internal_data_dir(self) -> Path:
        """Get the internal_data directory."""
        return get_germlines_base_dir() / "igblast" / "Ig" / "internal_data"

    @pytest.fixture
    def database_dir(self) -> Path:
        """Get the database directory."""
        return get_germlines_base_dir() / "igblast" / "database"

    def test_no_symlinks_in_internal_data(self, internal_data_dir: Path) -> None:
        """Test that no symlinks exist in internal_data directories."""
        if not internal_data_dir.exists():
            pytest.skip("internal_data directory not built")

        symlinks_found = []
        for species_dir in internal_data_dir.iterdir():
            if species_dir.is_dir():
                for item in species_dir.iterdir():
                    if item.is_symlink():
                        symlinks_found.append(str(item))

        assert (
            len(symlinks_found) == 0
        ), f"Symlinks found in internal_data: {symlinks_found}"

    def test_combined_fasta_exists_for_human(self, internal_data_dir: Path) -> None:
        """Test that combined FASTA file exists for human."""
        human_internal = internal_data_dir / "human"
        if not human_internal.exists():
            pytest.skip("Human internal_data not built")

        combined_fasta = human_internal / "human_V.fasta"
        assert combined_fasta.exists(), "Combined FASTA should exist"

        # Should NOT have separate D, J, C files in internal_data
        assert not (human_internal / "human_D.fasta").exists(), (
            "Separate D file should not exist in internal_data"
        )
        assert not (human_internal / "human_J.fasta").exists(), (
            "Separate J file should not exist in internal_data"
        )
        assert not (human_internal / "human_C.fasta").exists(), (
            "Separate C file should not exist in internal_data"
        )

    def test_combined_fasta_contains_multiple_segments(
        self, internal_data_dir: Path
    ) -> None:
        """Test that combined FASTA contains V, D, J, and C genes."""
        human_internal = internal_data_dir / "human"
        combined_fasta = human_internal / "human_V.fasta"

        if not combined_fasta.exists():
            pytest.skip("Human internal_data not built")

        # Read and parse FASTA
        content = combined_fasta.read_text()
        headers = [line for line in content.split("\n") if line.startswith(">")]

        # Check for different segment types
        has_v = any("IGHV" in h or "IGKV" in h or "IGLV" in h for h in headers)
        has_d = any("IGHD" in h for h in headers)
        has_j = any("IGHJ" in h or "IGKJ" in h or "IGLJ" in h for h in headers)
        has_c = any("IGHA" in h or "IGHG" in h or "IGHM" in h or "IGKC" in h for h in headers)

        assert has_v, "Combined FASTA should contain V genes"
        assert has_d, "Combined FASTA should contain D genes"
        assert has_j, "Combined FASTA should contain J genes"
        assert has_c, "Combined FASTA should contain C genes"

        # Should have significant number of sequences
        assert len(headers) > 100, f"Expected 100+ sequences, got {len(headers)}"

    def test_blast_database_files_exist(self, internal_data_dir: Path) -> None:
        """Test that BLAST database files are created from combined FASTA."""
        human_internal = internal_data_dir / "human"
        if not human_internal.exists():
            pytest.skip("Human internal_data not built")

        # Check for BLAST database files
        expected_extensions = [".nhr", ".nin", ".nsq"]
        for ext in expected_extensions:
            db_file = human_internal / f"human_V{ext}"
            assert db_file.exists(), f"BLAST database file {db_file.name} should exist"

    def test_ndm_file_exists(self, internal_data_dir: Path) -> None:
        """Test that NDM file is generated correctly."""
        human_internal = internal_data_dir / "human"
        if not human_internal.exists():
            pytest.skip("Human internal_data not built")

        ndm_file = human_internal / "human.ndm.imgt"
        assert ndm_file.exists(), "NDM file should exist"

        # NDM file should have entries
        content = ndm_file.read_text()
        lines = [l for l in content.strip().split("\n") if l]
        assert len(lines) > 0, "NDM file should have entries"


class TestGermlineDataPaths:
    """Tests for GermlineData path configuration."""

    def test_germline_data_paths_point_to_database(self) -> None:
        """Test that V/D/J/C dirs point to database/, not internal_data/."""
        from sadie.airr.igblast.germline import GermlineData

        try:
            gd = GermlineData("human")
        except (FileNotFoundError, ValueError) as e:
            pytest.skip(f"Human germline data not available: {e}")

        # V/D/J/C should point to database/
        assert "database" in str(gd.v_gene_dir), (
            f"v_gene_dir should point to database/, got {gd.v_gene_dir}"
        )
        assert "database" in str(gd.d_gene_dir), (
            f"d_gene_dir should point to database/, got {gd.d_gene_dir}"
        )
        assert "database" in str(gd.j_gene_dir), (
            f"j_gene_dir should point to database/, got {gd.j_gene_dir}"
        )
        assert "database" in str(gd.c_gene_dir), (
            f"c_gene_dir should point to database/, got {gd.c_gene_dir}"
        )

    def test_germline_data_igdata_points_to_ig(self) -> None:
        """Test that igdata points to Ig/ directory (contains internal_data)."""
        from sadie.airr.igblast.germline import GermlineData

        try:
            gd = GermlineData("human")
        except (FileNotFoundError, ValueError) as e:
            pytest.skip(f"Human germline data not available: {e}")

        # igdata should point to Ig/ (not internal_data directly)
        assert gd.igdata.name == "Ig", f"igdata should point to Ig/, got {gd.igdata}"
        assert (gd.igdata / "internal_data").exists(), (
            "igdata/internal_data should exist"
        )


class TestBuildInternalDataScript:
    """Tests for the build_internal_data script functions."""

    def test_build_combined_fasta(self, tmp_path: Path) -> None:
        """Test build_combined_fasta creates combined file with deduplication."""
        from sadie.germlines.scripts.build_internal_data import build_combined_fasta

        # Create mock database directory with test FASTAs
        database_dir = tmp_path / "database" / "test"
        database_dir.mkdir(parents=True)
        internal_data_dir = tmp_path / "internal_data" / "test"
        internal_data_dir.mkdir(parents=True)

        # Create test V FASTA
        v_fasta = database_dir / "test_V.fasta"
        v_fasta.write_text(">IGHV1-1*01\nATGC\n>IGHV1-2*01\nGCTA\n")

        # Create test D FASTA
        d_fasta = database_dir / "test_D.fasta"
        d_fasta.write_text(">IGHD1-1*01\nAAA\n")

        # Create test J FASTA
        j_fasta = database_dir / "test_J.fasta"
        j_fasta.write_text(">IGHJ1*01\nTTT\n")

        # Create test C FASTA with duplicate D entry
        c_fasta = database_dir / "test_C.fasta"
        c_fasta.write_text(">IGHA1*01\nCCC\n>IGHD1-1*01\nAAA\n")  # Duplicate!

        # Build combined FASTA
        combined = build_combined_fasta("test", database_dir, internal_data_dir)

        assert combined.exists()
        content = combined.read_text()

        # Should have 5 unique sequences (duplicate IGHD1-1*01 should be skipped)
        # V: 2, D: 1, J: 1, C: 1 (IGHA1*01 only, IGHD1-1*01 is duplicate)
        headers = [l for l in content.split("\n") if l.startswith(">")]
        assert len(headers) == 5, f"Expected 5 sequences, got {len(headers)}: {headers}"

        # Verify sequences are present
        assert ">IGHV1-1*01" in content
        assert ">IGHV1-2*01" in content
        assert ">IGHD1-1*01" in content
        assert ">IGHJ1*01" in content
        assert ">IGHA1*01" in content

        # Verify IGHD1-1*01 appears only once (from D file, not C file)
        assert content.count(">IGHD1-1*01") == 1, "Duplicate should have been removed"

    def test_build_blast_db(self, tmp_path: Path) -> None:
        """Test build_blast_db creates BLAST database files."""
        from sadie.germlines.scripts.build_internal_data import build_blast_db
        import shutil

        # Skip if makeblastdb not available
        if not shutil.which("makeblastdb"):
            # Check for bundled binary
            import platform
            system = platform.system().lower()
            bundled = Path(__file__).parent.parent.parent.parent / "src/sadie/reference/bin" / system / "makeblastdb"
            if not bundled.exists():
                pytest.skip("makeblastdb not available")

        # Create test FASTA
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(">seq1\nATGCATGC\n>seq2\nGCTAGCTA\n")

        # Build database
        db_prefix = tmp_path / "test_db"
        result = build_blast_db(fasta_path, db_prefix)

        assert result is True, "build_blast_db should return True on success"

        # Check database files exist
        assert (tmp_path / "test_db.nhr").exists() or (tmp_path / "test_db.nin").exists()
