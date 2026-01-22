import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock


class TestGermlineDataLegacyAPI:
    def test_germline_data_import(self):
        from sadie.airr.igblast.germline import GermlineData

        assert GermlineData is not None

    def test_germline_data_has_required_attributes(self, tmp_path):
        (tmp_path / "Ig/blastdb/human").mkdir(parents=True)
        (tmp_path / "Ig/blastdb/human/human_V.nhr").write_text("")
        (tmp_path / "Ig/blastdb/human/human_V.nin").write_text("")
        (tmp_path / "Ig/blastdb/human/human_D.nhr").write_text("")
        (tmp_path / "Ig/blastdb/human/human_D.nin").write_text("")
        (tmp_path / "Ig/blastdb/human/human_J.nhr").write_text("")
        (tmp_path / "Ig/blastdb/human/human_J.nin").write_text("")
        (tmp_path / "aux_db/imgt").mkdir(parents=True, exist_ok=True)
        (tmp_path / "aux_db/imgt/human_gl.aux").write_text("dummy aux data")
        (tmp_path / "Ig/internal_data/human").mkdir(parents=True)

        with patch("sadie.airr.igblast.germline._use_germlines_module", return_value=False):
            pass

    def test_get_available_datasets_returns_set(self):
        from sadie.airr.igblast.germline import GermlineData

        datasets = GermlineData.get_available_datasets()

        assert isinstance(datasets, set)

    def test_feature_flag_deprecation_warning(self):
        from sadie.airr.igblast.germline import _use_germlines_module
        import os

        original = os.environ.get("SADIE_USE_GERMLINES_MODULE")
        try:
            os.environ["SADIE_USE_GERMLINES_MODULE"] = "false"

            import importlib
            import sadie.airr.igblast.germline as germ_module
            importlib.reload(germ_module)

            result = germ_module._use_germlines_module()
            assert result is False
        finally:
            if original is not None:
                os.environ["SADIE_USE_GERMLINES_MODULE"] = original
            elif "SADIE_USE_GERMLINES_MODULE" in os.environ:
                del os.environ["SADIE_USE_GERMLINES_MODULE"]

    def test_feature_flag_default_true(self):
        from sadie.airr.igblast.germline import _use_germlines_module
        import os

        original = os.environ.get("SADIE_USE_GERMLINES_MODULE")
        try:
            if "SADIE_USE_GERMLINES_MODULE" in os.environ:
                del os.environ["SADIE_USE_GERMLINES_MODULE"]

            import importlib
            import sadie.airr.igblast.germline as germ_module
            importlib.reload(germ_module)

            result = germ_module._use_germlines_module()
            assert result is True
        finally:
            if original is not None:
                os.environ["SADIE_USE_GERMLINES_MODULE"] = original


class TestGermlineDataPaths:
    def test_v_gene_dir_attribute_exists(self, tmp_path):
        from sadie.airr.igblast.germline import GermlineData

        (tmp_path / "Ig/blastdb/human").mkdir(parents=True)
        for seg in ["V", "D", "J"]:
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nhr").write_text("")
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nin").write_text("")
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nsq").write_text("")
        (tmp_path / "aux_db/imgt").mkdir(parents=True)
        (tmp_path / "aux_db/imgt/human_gl.aux").write_text("dummy")
        (tmp_path / "Ig/internal_data/human").mkdir(parents=True)

        gd = GermlineData("human", database_dir=tmp_path)

        assert hasattr(gd, "v_gene_dir")
        assert hasattr(gd, "d_gene_dir")
        assert hasattr(gd, "j_gene_dir")

    def test_aux_path_attribute_exists(self, tmp_path):
        from sadie.airr.igblast.germline import GermlineData

        (tmp_path / "Ig/blastdb/human").mkdir(parents=True)
        for seg in ["V", "D", "J"]:
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nhr").write_text("")
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nin").write_text("")
            (tmp_path / f"Ig/blastdb/human/human_{seg}.nsq").write_text("")
        (tmp_path / "aux_db/imgt").mkdir(parents=True)
        (tmp_path / "aux_db/imgt/human_gl.aux").write_text("dummy")
        (tmp_path / "Ig/internal_data/human").mkdir(parents=True)

        gd = GermlineData("human", database_dir=tmp_path)

        assert hasattr(gd, "aux_path")
