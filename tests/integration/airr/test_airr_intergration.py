import pandas as pd
from sadie.airr import Airr


def test_vs_oas_sampling():
    """
    Test airr vs the OAS airr datasets. I have chosen 15K heavy and 10K light to compare against
    """
    end_point = "https://sadie.s3.us-east-2.amazonaws.com/integration/OAS_sample_subsample.bz2"

    # uzing bzip for max compression to minimize egress S3 expense
    df_bz2 = pd.read_csv(end_point, index_col=0).reset_index()

    # Airr api
    airr_api = Airr(species="human", database="imgt", functional="all")
    airr_api.run_dataframe(
        df_bz2,
    )


# def test_cli():
#     """Confirm the CLI works as expecte"""
#     runner = CliRunner()
#     fastas = get_fasta_inputs("*")
#     for query_file in fastas:
#         with tempfile.NamedTemporaryFile() as tmpfile:
#             logger.debug(tmpfile.name)
#             result = runner.invoke(app.run_airr, ["--query", query_file, "-o", tmpfile.name])
#             assert result.exit_code == 0
#             assert os.path.exists(tmpfile.name)
#     assert not os.path.exists(tmpfile.name)
#     for query_file in fastas:
#         with tempfile.NamedTemporaryFile() as tmpfile:
#             logger.debug(tmpfile.name)
#             result = runner.invoke(app.run_airr, ["--query", query_file, "-s", "dog", "-o", tmpfile.name])
#             assert result.exit_code == 0
#             assert os.path.exists(tmpfile.name)
#         assert not os.path.exists(tmpfile.name)
