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
