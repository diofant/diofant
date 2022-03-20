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

    c = Console("python -m diofant -a --simple-prompt --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('a\r\n') == 3
    assert c.expect_exact('\r\nOut[1]: a\r\n\r\nIn [2]: ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nOut[2]: 1/2\r\n\r\nIn [3]: ') == 0


def test_ipython_console_bare_division_noauto():
    pytest.importorskip('IPython')

    c = Console('python -m diofant --simple-prompt '
                "--no-wrap-division --colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nOut[1]: 0.5\r\n\r\nIn [2]: ') == 0
    assert c.send('spam\r\n') == 6
    assert c.expect(".*NameError: name 'spam' "
                    'is not defined\r\n\r\nIn [\\[]3[]]: ') == 0
    assert c.send('x\r\n') == 3
    assert c.expect_exact('\r\nOut[3]: x\r\n\r\nIn [4]: ') == 0
