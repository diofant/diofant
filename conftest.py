import pytest

from diofant.core.cache import clear_cache


def pytest_report_header(config):
    from diofant.utilities.misc import ARCH
    from diofant.core.cache import USE_CACHE
    from diofant.core.compatibility import GROUND_TYPES, HAS_GMPY
    s = "architecture: %s\n" % ARCH
    s += "cache:        %s\n" % USE_CACHE
    version = ''
    if GROUND_TYPES == 'gmpy':
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
    from diofant import init_printing

    # hook our nice, hash-stable strprinter
    init_printing(pretty_print=False, use_unicode=False)

    # doctest restore sys.displayhook from __displayhook__,
    # see https://bugs.python.org/issue26092.
    sys.__displayhook__ = sys.displayhook
