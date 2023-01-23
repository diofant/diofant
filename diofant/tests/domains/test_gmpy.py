"""Tests for gmpy2-based domains."""

import pytest

from diofant import (CC, FF, CoercionFailedError, FF_gmpy, FF_python,
                     PythonRational, QQ_gmpy, QQ_python, ZZ_gmpy, ZZ_python)


__all__ = ()

gmpy = pytest.importorskip('gmpy2')


def test_convert():
    assert ZZ_gmpy.finite_field(2) == FF_gmpy(2)
    assert ZZ_python.finite_field(2) == FF_python(2)

    F3 = FF(3)
    F3_gmpy = FF_gmpy(3)
    F3_python = FF_python(3)

    assert F3.convert(gmpy.mpz(2)) == F3.dtype(2)
    assert F3.convert(gmpy.mpq(2, 1)) == F3.dtype(2)
    pytest.raises(CoercionFailedError, lambda: F3.convert(gmpy.mpq(1, 2)))

    assert ZZ_gmpy.convert(F3_python(1)) == ZZ_gmpy.dtype(1)
    assert ZZ_gmpy.convert(F3_gmpy(1)) == ZZ_gmpy.dtype(1)

    assert ZZ_gmpy.convert(PythonRational(2)) == ZZ_gmpy.dtype(2)
    pytest.raises(CoercionFailedError,
                  lambda: ZZ_gmpy.convert(PythonRational(2, 3)))

    assert QQ_python.convert(gmpy.mpz(3)) == QQ_python.dtype(3)
    assert QQ_python.convert(gmpy.mpq(2, 3)) == QQ_python.dtype(2, 3)

    assert QQ_gmpy.convert(PythonRational(2, 3)) == QQ_gmpy.dtype(2, 3)

    assert ZZ_python.convert(F3_gmpy(1)) == ZZ_python.dtype(1)
    assert ZZ_python.convert(gmpy.mpz(3)) == ZZ_python.dtype(3)
    assert ZZ_python.convert(gmpy.mpq(3, 1)) == ZZ_python.dtype(3)
    pytest.raises(CoercionFailedError,
                  lambda: ZZ_python.convert(gmpy.mpq(3, 2)))

    assert CC.convert(gmpy.mpz(3)) == CC(3)
    assert CC.convert(gmpy.mpq(1, 2)) == CC(0.5)
