from sadie.reference import Reference

reference_path = "reference.yml"
reference_object = Reference.parse_yaml(reference_path)

outpath = "my_output_database_path"
germline_path = reference_object.make_airr_database(outpath)
