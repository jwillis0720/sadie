#!/usr/bin/sh
pip install .
pip install wheel twine setuptools
rm -rf dist/*
python setup.py sdist bdist_wheel
twine upload --verbose --repository pypi  dist/*.whl
