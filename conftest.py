import re
import sys
import warnings

import hypothesis
import pytest

import diofant


collect_ignore = ["setup.py"]

sp = re.compile(r'([0-9]+)/([1-9][0-9]*)')

hypothesis.settings.register_profile("default",
                                     hypothesis.settings(max_examples=100))
hypothesis.settings.register_profile("debug",
                                     hypothesis.settings(max_examples=100,
                                                         verbosity=hypothesis.Verbosity.verbose))


def process_split(session, config, items):
    split = config.getoption("--split")
    if not split:
        return
    m = sp.match(split)
    if not m:
        raise ValueError("split must be a string of the form a/b "
                         "where a and b are ints.")
    i, t = map(int, m.groups())
    start, end = (i - 1)*len(items)//t, i*len(items)//t

    if i < t:
        del items[end:]
    del items[:start]


def pytest_report_header(config):
    return """
cache: %s
ground types: %s
""" % (diofant.core.cache.USE_CACHE, diofant.core.compatibility.GROUND_TYPES)


def pytest_addoption(parser):
    parser.addoption("--split", action="store", default="", help="split tests")


def pytest_collection_modifyitems(session, config, items):
    process_split(session, config, items)


@pytest.fixture(autouse=True, scope='module')
def file_clear_cache():
    diofant.core.cache.clear_cache()


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


@pytest.fixture(autouse=True, scope='session')
def enable_mpl_agg_backend():
    try:
        import matplotlib as mpl
        mpl.use('Agg')
        del mpl
    except ImportError:
        pass


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    for sym in (diofant.symbols('a b c d x y z t') +
                diofant.symbols('k m n', integer=True) +
                diofant.symbols('f g h', cls=diofant.Function)):
        doctest_namespace[str(sym)] = sym
    for name in dir(diofant):
        doctest_namespace[name] = getattr(diofant, name)
