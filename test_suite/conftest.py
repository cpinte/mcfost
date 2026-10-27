# Long tests (marked @pytest.mark.slow) only run with "pytest --run-slow", e.g. the
# optically thick Pascucci et al. 2004 benchmarks: the default run (CI, pull requests)
# skips them.

import pytest


def pytest_addoption(parser):
    parser.addoption("--run-slow", action="store_true", default=False,
                     help="also run the long tests (marked slow)")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long test, only run with --run-slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="long test: use --run-slow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
