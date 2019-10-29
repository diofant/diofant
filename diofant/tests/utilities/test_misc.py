from diofant.tests.core.test_cache import clear_imports  # noqa: F401
from diofant.utilities.misc import debug


__all__ = ()


def test_debug_on(clear_imports, monkeypatch, capsys):  # noqa: F811
    monkeypatch.setenv('DIOFANT_DEBUG', 'True')

    debug('Hi')

    assert capsys.readouterr().err == 'Hi\n'


def test_debug_off(capsys):
    debug('Hi')

    assert capsys.readouterr().err == ''
