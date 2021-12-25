#!/usr/bin/sh
coverage run --source=sadie -m pytest tests/unit/
coverage html
coverage xml
coverage report
pytest -x -n0 -v tests/integration
