"""Tests for configuration helpers."""

from diofant.config import configure, query


def test_configure(monkeypatch):
    assert query('heu_gcd_max') == 6
    monkeypatch.setenv('DIOFANT_HEU_GCD_MAX', '7')
    configure()
    assert query('heu_gcd_max') == 7
    monkeypatch.setenv('DIOFANT_HEU_GCD_MAX', '1^2')
    configure()
    assert query('heu_gcd_max') == 6
    monkeypatch.setenv('DIOFANT_HEU_GCD_MAX', '"abcd"')
    configure()
    assert query('heu_gcd_max') == 6


def test_configure2():
    configure()
    assert query('heu_gcd_max') == 6
