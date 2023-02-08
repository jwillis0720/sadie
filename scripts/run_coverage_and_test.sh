#!/usr/bin/sh
coverage run --source=sadie -m pytest tests/
coverage html
coverage xml
coverage report
# python.any for wheel
