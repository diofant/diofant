import pytest

from diofant.polys.polyerrors import CoercionFailed
from diofant.domains import (FF, QQ_python, QQ_gmpy, ZZ_gmpy,
                             FF_gmpy, PythonRational, ZZ_python, CC)
from diofant.external import import_module

__all__ = ()

gmpy = import_module('gmpy2')


@pytest.mark.skipif(gmpy is None, reason="no gmpy")
def test_convert():
    F3 = FF(3)
    F3g = FF_gmpy(3)

    assert F3.convert(gmpy.mpz(2)) == F3.dtype(2)
    assert F3.convert(gmpy.mpq(2, 1)) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(gmpy.mpq(1, 2)))

    assert ZZ_gmpy().convert(F3(1)) == ZZ_gmpy().dtype(1)

    assert ZZ_gmpy().convert(PythonRational(2)) == ZZ_gmpy().dtype(2)
    pytest.raises(CoercionFailed,
                  lambda: ZZ_gmpy().convert(PythonRational(2, 3)))

    assert QQ_python().convert(gmpy.mpz(3)) == QQ_python().dtype(3)
    assert QQ_python().convert(gmpy.mpq(2, 3)) == QQ_python().dtype(2, 3)

    assert QQ_gmpy().convert(PythonRational(2, 3)) == QQ_gmpy().dtype(2, 3)

    assert ZZ_python().convert(F3g(1)) == ZZ_python().dtype(1)
    assert ZZ_python().convert(gmpy.mpz(3)) == ZZ_python().dtype(3)
    assert ZZ_python().convert(gmpy.mpq(3, 1)) == ZZ_python().dtype(3)
    pytest.raises(CoercionFailed,
                  lambda: ZZ_python().convert(gmpy.mpq(3, 2)))

    assert CC.convert(gmpy.mpz(3)) == CC(3)
    assert CC.convert(gmpy.mpq(1, 2)) == CC(0.5)
