import re

import pytest

from sympy.core.cache import clear_cache


def pytest_report_header(config):
    from sympy.utilities.misc import ARCH
    s = "architecture: %s\n" % ARCH
    from sympy.core.cache import USE_CACHE
    s += "cache:        %s\n" % USE_CACHE
    from sympy.core.compatibility import GROUND_TYPES, HAS_GMPY
    version = ''
    if GROUND_TYPES == 'gmpy':
        if HAS_GMPY == 1:
            import gmpy
        elif HAS_GMPY == 2:
            import gmpy2 as gmpy
        version = gmpy.version()
    s += "ground types: %s %s\n" % (GROUND_TYPES, version)
    return s


def pytest_terminal_summary(terminalreporter):
    if (terminalreporter.stats.get('error', None) or
            terminalreporter.stats.get('failed', None)):
        terminalreporter.write_sep(
            ' ', 'DO *NOT* COMMIT!', red=True, bold=True)


@pytest.fixture(autouse=True, scope='module')
def file_clear_cache():
    clear_cache()


@pytest.fixture(autouse=True, scope='module')
def check_disabled(request):
    if getattr(request.module, 'disabled', False):
        pytest.skip("test requirements not met.")


@pytest.fixture(autouse=True, scope='session')
def set_displayhook():
    import sys
    from sympy import init_printing

    # hook our nice, hash-stable strprinter
    init_printing(pretty_print=False, use_unicode=False)

    # doctest restore sys.displayhook from __displayhook__,
    # see https://bugs.python.org/issue26092.
    sys.__displayhook__ = sys.displayhook


sp = re.compile(r'([0-9]+)/([1-9][0-9]*)')


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
        # remove elements from end of list first
        del items[end:]
    del items[:start]


def pytest_addoption(parser):
    parser.addoption("--split", action="store", default="",
                     help="split tests")


def pytest_collection_modifyitems(session, config, items):
    """ pytest hook. """
    # handle splits
    process_split(session, config, items)
