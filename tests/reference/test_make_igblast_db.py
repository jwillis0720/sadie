import glob
import os
import shutil

import pandas as pd
import pytest
from click.testing import CliRunner
from pybody.reference import app


def test_make_igblast_reference():
    """Confirm the CLI works as expecte"""
    runner = CliRunner(echo_stdin=True)
    # with tempfile.TemporaryDirectory() as tmpdir:
    directory = "/tmp/igblast_dirs"
    if os.path.exists(directory):
        shutil.rmtree(directory)

    os.makedirs(directory)
    result = runner.invoke(app.make_igblast_reference, ["--outpath", directory], catch_exceptions=True)
    if result.exit_code != 0:
        print(result)
    assert result.exit_code == 0
    assert os.path.exists(directory)

    directories_created = glob.glob(directory + "/*")
    assert sorted(directories_created) == sorted(
        [
            "/tmp/igblast_dirs/blastdb",
            "/tmp/igblast_dirs/aux_data",
            "/tmp/igblast_dirs/internal_data",
        ]
    )
    blast_dbs_created = glob.glob("/tmp/igblast_dirs/blastdb/*")
    assert sorted(blast_dbs_created) == sorted(["/tmp/igblast_dirs/blastdb/TCR", "/tmp/igblast_dirs/blastdb/Ig"])
    shutil.rmtree(directory)
    assert not os.path.exists(directory)
