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
