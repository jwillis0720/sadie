from sadie.reference import References

reference_path = "reference.yml"
references_object = References.from_yaml(reference_path)

outpath = "my_output_database_path"
germline_path = reference_object.make_airr_database(outpath)
