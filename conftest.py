import sys
import warnings

import pytest

from diofant.core.cache import clear_cache, USE_CACHE
from diofant.core.compatibility import GROUND_TYPES


def pytest_report_header(config):
    return """
cache: %s
ground types: %s
""" % (USE_CACHE, GROUND_TYPES)


@pytest.fixture(autouse=True, scope='module')
def file_clear_cache():
    clear_cache()


@pytest.fixture(autouse=True, scope='module')
def check_disabled(request):
    if getattr(request.module, 'disabled', False):
        pytest.skip("test requirements not met.")


@pytest.fixture(autouse=True, scope='session')
def set_displayhook():
    sys.__displayhook__ = sys.displayhook  # https://bugs.python.org/26092


@pytest.fixture(autouse=True, scope='session')
def enable_deprecationwarnings():
    warnings.simplefilter('error', DeprecationWarning)
