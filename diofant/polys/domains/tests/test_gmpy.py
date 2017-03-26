import pytest

from diofant.polys.polyerrors import CoercionFailed
from diofant.polys.domains import FF, QQ_python
from diofant.external import import_module

__all__ = ()

gmpy = import_module('gmpy2')


@pytest.mark.skipif(gmpy is None, reason="no gmpy")
def test_convert():
    F3 = FF(3)

    assert F3.convert(gmpy.mpz(2)) == F3.dtype(2)
    assert F3.convert(gmpy.mpq(2, 1)) == F3.dtype(2)
    pytest.raises(CoercionFailed, lambda: F3.convert(gmpy.mpq(1, 2)))

    assert QQ_python().convert(gmpy.mpz(3)) == QQ_python().dtype(3)
    assert QQ_python().convert(gmpy.mpq(2, 3)) == QQ_python().dtype(2, 3)
