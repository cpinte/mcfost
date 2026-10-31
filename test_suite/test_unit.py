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
    # run from test_suite/: test_mcfost2phantom reads unit/mcfost2phantom.para
    result = subprocess.run([os.path.join(_unit_dir, name)], capture_output=True, text=True, cwd=_test_dir)
    print(result.stdout)
    assert result.returncode == 0, result.stdout + result.stderr


def test_mcfost2phantom_overflow(build_unit_tests):
    # library called without max_particles, then with more particles than on the first
    # call: mcfost must stop with an error instead of writing out of its arrays
    result = subprocess.run([os.path.join(_unit_dir, "test_mcfost2phantom"), "overflow"],
                            capture_output=True, text=True, cwd=_test_dir)
    print(result.stdout)
    assert result.returncode != 0, result.stdout
    assert "more particles than allocated by mcfost" in result.stdout
