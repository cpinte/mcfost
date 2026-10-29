# Long tests (marked @pytest.mark.slow), e.g. the optically thick Pascucci et al. 2004
# benchmarks, run locally but are skipped on CI (pull requests). Skip them locally with
# "pytest -m 'not slow'".

import os
import pytest

# True on CI (GitHub Actions sets CI=true): skips the long tests and those needing too much memory
ON_CI = os.environ.get('CI', None) == 'true'


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long test, skipped on CI")


def pytest_collection_modifyitems(config, items):
    if not ON_CI:
        return
    skip_slow = pytest.mark.skip(reason="long test: skipped on CI")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
