"""The setup script."""
import sys
import subprocess

from setuptools import setup, find_packages
from setuptools.command.test import test


with open("README.md") as readme_file:
    readme = readme_file.read()
requirements = open("requirements.txt").readlines()


class PyTest(test):
    """Define what happens when we run python setup.py test"""

    def initialize_options(self):
        test.initialize_options(self)
        self.pytest_args = "tests/"

    def run_tests(self):
        # Import here, because outside the eggs aren't loaded
        # import shlex
        import pytest

        err_number = pytest.main(["-vvv"])
        sys.exit(err_number)


setup(
    name="pibody",
    version="0.0.1",
    python_requires=">=3.6",
    description="PiBody: The Complete Antibody Library",
    author="Jordan R Willis",
    author_email="jwillis0720@gmail.com",
    url="https://github.com/jwillis0720/",
    packages=find_packages(include=["pibody*"], exclude=["tests*", "docs"]),
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords=["airr", "annotating", "antibodies"],
    entry_points={
        "console_scripts": ["airr=pibody.airr.app:run_airr"],
    },
    test_suite="tests",
    cmdclass={"tests": PyTest},
)
