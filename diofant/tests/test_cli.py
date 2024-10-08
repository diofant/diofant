"""Tests for the Command-Line Interface."""

import time

import pexpect
import pytest


__all__ = ()


class Console(pexpect.spawn):
    """Spawned console for testing."""

    def __init__(self, command, timeout=60):
        super().__init__(command, timeout=timeout, encoding='utf-8')

    def __del__(self):
        self.send('exit()\r\n')
        time.sleep(10)  # a delay to allow coverage finish work
        if self.isalive():
            self.terminate(force=True)


def test_bare_console():
    c = Console('python -m diofant --no-ipython')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1 + 2\r\n') == 7
    assert c.expect_exact('3\r\n>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('1/2\r\n>>> ') == 0
    assert c.send('x\r\n') == 3
    assert c.expect_exact('x\r\n>>> ') == 0
    assert c.expect_exact('>>> ') == 0
    assert c.send('ℕ = 1\r\n') == 9
    assert c.expect_exact('\r\n>>> ') == 0
    assert c.send('N = 2\r\n') == 7
    assert c.expect_exact('\r\n>>> ') == 0
    assert c.send('ℕ\r\n') == 5
    assert c.expect_exact('2\r\n>>> ') == 0


def test_bare_console_bare_division():
    c = Console('python -m diofant --no-ipython --no-wrap-division')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('0.5\r\n>>> ') == 0


def test_bare_console_with_auto():
    c = Console('python -m diofant -a --no-ipython')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1 + 2\r\n') == 7
    assert c.expect_exact('3\r\n>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('1/2\r\n>>> ') == 0
    assert c.send('x\r\n') == 3
    assert c.expect_exact('x\r\n>>> ') == 0
    assert c.send('(a/2 + \r\n1)\r\n') == 13
    assert c.expect_exact('a    \r\n─ + 1\r\n2    \r\n>>> ') == 0


def test_bare_console_with_unicode_identifiers():
    c = Console('python -m diofant --no-ipython --unicode-identifiers')

    assert c.expect_exact('>>> ') == 0
    assert c.send('ℕ = 1\r\n') == 9
    assert c.expect_exact('\r\n>>> ') == 0
    assert c.send('N = 2\r\n') == 7
    assert c.expect_exact('\r\n>>> ') == 0
    assert c.send('ℕ\r\n') == 5
    assert c.expect_exact('1\r\n>>> ') == 0


def test_bare_console_without_ipython():
    try:
        import IPython
        del IPython
        pytest.skip('IPython is available')
    except ImportError:
        pass

    c = Console('python -m diofant')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1 + 2\r\n') == 7
    assert c.expect_exact('3\r\n>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('1/2\r\n>>> ') == 0
    assert c.send('x\r\n') == 3
    assert c.expect_exact('x\r\n>>> ') == 0


def test_ipython_console():
    pytest.importorskip('IPython')

    c = Console('python -m diofant -a --unicode-identifiers '
                "--simple-prompt --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('a\r\n') == 3
    assert c.expect_exact('\r\nOut[1]: a\r\n\r\nIn [2]: ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nOut[2]: 1/2\r\n\r\nIn [3]: ') == 0
    assert c.send('ℕ = 1\r\n') == 9
    assert c.expect_exact('\r\nIn [4]: ') == 0
    assert c.send('N = 2\r\n') == 7
    assert c.expect_exact('\r\nIn [5]: ') == 0
    assert c.send('ℕ\r\n') == 5
    assert c.expect_exact('Out[5]: 1\r\n\r\nIn [6]: ') == 0


def test_ipython_console_bare_division_noauto():
    pytest.importorskip('IPython')

    c = Console('python -m diofant --simple-prompt '
                "--no-wrap-division --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nOut[1]: 0.5\r\n\r\nIn [2]: ') == 0
    assert c.send('spam\r\n') == 6
    assert c.expect(".*NameError: name 'spam' "
                    'is not defined\r\n.*\r\nIn [\\[]3[]]: ') == 0
    assert c.send('x\r\n') == 3
    assert c.expect_exact('\r\nOut[3]: x\r\n\r\nIn [4]: ') == 0
    assert c.send('ℕ = 1\r\n') == 9
    assert c.expect_exact('\r\nIn [5]: ') == 0
    assert c.send('N = 2\r\n') == 7
    assert c.expect_exact('\r\nIn [6]: ') == 0
    assert c.send('ℕ\r\n') == 5
    assert c.expect_exact('Out[6]: 2\r\n\r\nIn [7]: ') == 0


def test_ipython_console_wrap_floats():
    pytest.importorskip('IPython')

    c = Console('python -m diofant --simple-prompt '
                "--wrap-floats --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('repr(10.9)\r\n') == 12
    assert c.expect_exact("\r\nOut[1]: \"Float('10.9004', dps=3)\"\r\n\r\nIn [2]: ") == 0


def test_bare_console_wrap_floats():
    c = Console('python -m diofant --simple-prompt --no-ipython '
                "--wrap-floats --colors 'NoColor'")

    assert c.expect_exact('>>> ') == 0
    assert c.send('repr(10.9)\r\n') == 12
    assert c.expect_exact("\"Float('10.9004', dps=3)\"\r\n>>> ") == 0
    assert c.send('Float(1.) + 1\r\n') == 15
    assert c.expect_exact('2.00000000000000\r\n>>> ') == 0


def test_ipython_console_wrap_ints():
    pytest.importorskip('IPython')

    c = Console('python -m diofant --simple-prompt '
                "--wrap-ints --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('repr(10)\r\n') == 10
    assert c.expect_exact("\r\nOut[1]: \'Integer(10)\'\r\n\r\nIn [2]: ") == 0


def test_bare_console_wrap_ints():
    c = Console('python -m diofant --simple-prompt --no-ipython '
                "--wrap-ints --colors 'NoColor'")

    assert c.expect_exact('>>> ') == 0
    assert c.send('repr(10)\r\n') == 10
    assert c.expect_exact("\'Integer(10)\'\r\n>>> ") == 0


def test_diofant_version():
    c = Console('python -m diofant --version')

    assert c.expect(pexpect.EOF) == 0
    assert c.before.startswith('0.')
