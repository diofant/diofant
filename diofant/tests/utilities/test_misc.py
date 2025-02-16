import pytest

from diofant import nan, oo, zoo
from diofant.utilities import as_int
from diofant.utilities.decorator import no_attrs_in_subclass


__all__ = ()


def test_no_attrs_in_subclass():
    class A:
        x = 'test'

    A.x = no_attrs_in_subclass(A, A.x)

    class B(A):
        pass

    assert hasattr(A, 'x') is True
    assert hasattr(B, 'x') is False


def test_as_int():
    pytest.raises(ValueError, lambda: as_int(1.1))
    pytest.raises(ValueError, lambda: as_int([]))
    pytest.raises(ValueError, lambda: as_int(nan))
    pytest.raises(ValueError, lambda: as_int(oo))
    pytest.raises(ValueError, lambda: as_int(-oo))
    pytest.raises(ValueError, lambda: as_int(zoo))
