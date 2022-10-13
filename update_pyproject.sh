poetry add $( cat requirements-core.txt )
# poetry add --optional $( cat requirements-dev.txt )  # could replicate pip [] with this
poetry add --group dev  $( cat requirements-dev.txt )
