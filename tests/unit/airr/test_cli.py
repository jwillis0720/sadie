import os
from pathlib import Path

from click.testing import CliRunner

from sadie.app import airr, sadie
from tests.conftest import SadieFixture


def test_airr_cli(fixture_setup: SadieFixture):
    tmp_path = fixture_setup.tmp_path
    runner = CliRunner()

    # check that sadie can be invoked alone
    results = runner.invoke(sadie)
    assert results.exit_code == 0

    input_file = fixture_setup.get_catnap_heavy_nt()

    # this will be the default if we don't specify
    output_path = input_file.parent / Path(input_file.stem + ".tsv.gz")
    results = runner.invoke(airr, [str(input_file)])
    assert results.exit_code == 0
    os.remove(output_path)

    # explicitly specify the output file
    results = runner.invoke(airr, ["--skip-igl", "--skip-mutation", str(input_file), str(output_path)])
    assert results.exit_code == 0
    os.remove(output_path)

    # can we do feather
    tmp_out = tmp_path / "airr.feather"
    results = runner.invoke(airr, ["--skip-igl", "--skip-mutation", str(input_file), str(tmp_out)])
    assert Path(tmp_out).exists()
    assert results.exit_code == 0

    # can we do csv
    tmp_out = tmp_path / "airr.csv"
    results = runner.invoke(airr, ["--skip-igl", "--skip-mutation", str(input_file), str(tmp_out)])
    assert Path(tmp_out).exists()
    assert results.exit_code == 0

    # can we do gb
    tmp_out = tmp_path / "airr.gb"
    results = runner.invoke(airr, ["--skip-igl", "--skip-mutation", str(input_file), str(tmp_out)])
    assert Path(tmp_out).exists()
    assert results.exit_code == 0
