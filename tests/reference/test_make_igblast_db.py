import glob
import os
import shutil
import tempfile

import pandas as pd
import pytest
from click.testing import CliRunner
from sadie.reference import app


def test_make_igblast_reference():
    """Confirm the CLI works as expecte"""
    runner = CliRunner(echo_stdin=True)
    with tempfile.TemporaryDirectory(suffix="igblast_dir") as tmpdir:
        result = runner.invoke(app.make_igblast_reference, ["--outpath", tmpdir], catch_exceptions=True)
        if result.exit_code != 0:
            print(result)
        assert result.exit_code == 0
        assert os.path.exists(tmpdir)

        directories_created = glob.glob(tmpdir + "/*")
        assert sorted(directories_created) == sorted(
            [
                f"{tmpdir}/blastdb",
                f"{tmpdir}/aux_data",
                f"{tmpdir}/internal_data",
            ]
        )
        blast_dbs_created = glob.glob(f"{tmpdir}/blastdb/*")
        assert sorted(blast_dbs_created) == sorted([f"{tmpdir}/blastdb/TCR", f"{tmpdir}/blastdb/Ig"])
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    assert not os.path.exists(tmpdir)
