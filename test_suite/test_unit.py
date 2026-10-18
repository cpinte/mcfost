import os
import subprocess
import pytest

# Fortran unit tests in test_suite/unit/, built by "make unit_tests" in src/
# against libmcfost.a, i.e. with the same compiler (SYSTEM) and flags as mcfost.
# In CI they are built in the compile step (where the compiler environment is
# loaded), so the make call below has nothing to do; locally it builds them.

_test_dir = os.path.dirname(os.path.abspath(__file__))
_src_dir = os.path.join(_test_dir, "..", "src")
_unit_dir = os.path.join(_test_dir, "unit")

unit_tests = sorted(
    f[:-4] for f in os.listdir(_unit_dir) if f.startswith("test_") and f.endswith(".f90")
)


@pytest.fixture(scope="module")
def build_unit_tests():
    result = subprocess.run(["make", "-C", _src_dir, "unit_tests"], capture_output=True, text=True)
    assert result.returncode == 0, "make unit_tests failed:\n" + result.stdout + result.stderr


@pytest.mark.parametrize("name", unit_tests)
def test_unit(build_unit_tests, name):
    result = subprocess.run([os.path.join(_unit_dir, name)], capture_output=True, text=True)
    print(result.stdout)
    assert result.returncode == 0, result.stdout + result.stderr
