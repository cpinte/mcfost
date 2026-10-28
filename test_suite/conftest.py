# Long tests (marked @pytest.mark.slow), e.g. the optically thick Pascucci et al. 2004
# benchmarks, run locally but are skipped on CI (pull requests). Skip them locally with
# "pytest -m 'not slow'".

import os
import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: long test, skipped on CI")


def pytest_collection_modifyitems(config, items):
    if os.environ.get('CI', None) != 'true':
        return
    skip_slow = pytest.mark.skip(reason="long test: skipped on CI")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
