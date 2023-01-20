"""Pytest configuration and fixtures for the Diofant test suite."""

import random
import sys

import pytest

import diofant


def pytest_report_header(config):
    return f"""\nDiofant version: {diofant.__version__}
cache: {diofant.core.cache.USE_CACHE}
ground types: {diofant.core.compatibility.GROUND_TYPES}\n"""


def pytest_configure(config):
    config.addinivalue_line('markers', 'slow: marks tests as slow')
    config.addinivalue_line('markers', 'regression: marks a regression test')


def pytest_collection_modifyitems(items):
    for item in items:
        if 'issue' in item.nodeid:
            item.add_marker(pytest.mark.regression)


@pytest.fixture(autouse=True, scope='module')
def _file_clear_cache():
    diofant.core.cache.clear_cache()


@pytest.fixture(autouse=True, scope='session')
def _set_displayhook():
    sys.__displayhook__ = sys.displayhook  # https://bugs.python.org/26092


@pytest.fixture(autouse=True)
def _add_np(doctest_namespace):
    for sym in (diofant.symbols('a:d t x:z') +
                diofant.symbols('k m n', integer=True) +
                diofant.symbols('f:h', cls=diofant.Function)):
        doctest_namespace[str(sym)] = sym
    for name in dir(diofant):
        doctest_namespace[name] = getattr(diofant, name)
    random.seed(0)
