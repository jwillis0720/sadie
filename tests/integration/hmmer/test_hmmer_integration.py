import logging
import glob
import os
import tempfile
from itertools import product
from pathlib import Path

from click.testing import CliRunner
from sadie.app import renumbering
from sadie.renumbering import Renumbering

logger = logging.getLogger()


# cli helper
def _run_cli(p_tuple):
    with tempfile.NamedTemporaryFile() as tmpfile:
        cli_input = [
            "--query",
            str(p_tuple[0]),
            "-o",
            tmpfile.name,
            "-a",
            p_tuple[1],
            "-f",
            p_tuple[2],
            "-s",
            p_tuple[3],
            "-r",
            p_tuple[4],
            "-vvvvv",
        ]
        if not Renumbering.check_combination(p_tuple[3], p_tuple[4]):
            logger.info(f"skipping {p_tuple[4]}-{p_tuple[3]}")
            return True

        logger.info(f"CLI input {' '.join(cli_input)}")
        runner = CliRunner()
        result = runner.invoke(renumbering, cli_input)
        if result.exit_code != 0:
            logger.info(f"!!!STDERR {result}")

        assert result.exit_code == 0

        # numbering appends reults and alignment so glob
        out_file = Path(tmpfile.name).stem + "*"
    for f in glob.glob(out_file):
        os.remove(f)
    return True


def test_cli(fixture_setup):
    """Confirm the CLI works as expecte"""
    test_file_heavy = fixture_setup.get_catnap_heavy_aa()
    test_file_light = fixture_setup.get_catnap_light_aa()
    species = ["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "dog", "cat"]
    species = ",".join(species)
    ft = ["csv", "json", "feather"]
    schemes = ["imgt", "kabat", "chothia"]
    regions = ["imgt", "kabat", "chothia", "abm", "contact", "scdr"]
    products = product([test_file_heavy, test_file_light], [species], ft, schemes, regions)

    # pool = Pool()
    results = list(map(_run_cli, products))
    print(results)
