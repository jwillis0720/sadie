"""The setup script."""
import sys
from setuptools import setup
from setuptools.command.test import test


# From https://stackoverflow.com/questions/45150304/how-to-force-a-python-wheel-to-be-platform-specific-when-building-it
# Tell bdsit wheel to be platform specific even though we arent compiling extensions
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False

except ImportError:
    bdist_wheel = None


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


# great info
# https://stackoverflow.com/questions/50155464/using-pytest-with-a-src-layer
setup(
    test_suite="tests",
    cmdclass={"test": PyTest, "bdist_wheel": bdist_wheel},
)
