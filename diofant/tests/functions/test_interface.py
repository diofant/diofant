"""This test file tests the Diofant function interface."""

from diofant import Function, cos, limit, sin, sympify, tanh
from diofant.abc import x


__all__ = ()


def test_function_series1():
    """Create our new "sin" function."""

    class MyFunction(Function):

        def fdiff(self, argindex=1):
            return cos(self.args[0])

        @classmethod
        def eval(cls, arg):
            arg = sympify(arg)
            if arg == 0:
                return sympify(0)

    # Test that the taylor series is correct
    assert MyFunction(x).series(x, 0, 10) == sin(x).series(x, 0, 10)
    assert limit(MyFunction(x)/x, x, 0) == 1


def test_function_series2():
    """Create our new "cos" function."""

    class MyFunction2(Function):

        def fdiff(self, argindex=1):
            return -sin(self.args[0])

        @classmethod
        def eval(cls, arg):
            arg = sympify(arg)
            if arg == 0:
                return sympify(1)

    # Test that the taylor series is correct
    assert MyFunction2(x).series(x, 0, 10) == cos(x).series(x, 0, 10)


def test_function_series3():
    """
    Test our easy "tanh" function.

    This test tests two things:
      * that the Function interface works as expected and it's easy to use
      * that the general algorithm for the series expansion works even when the
        derivative is defined recursively in terms of the original function,
        since tanh(x).diff(x) == 1-tanh(x)**2
    """

    class MyTanh(Function):

        def fdiff(self, argindex=1):
            return 1 - MyTanh(self.args[0])**2

        @classmethod
        def eval(cls, arg):
            arg = sympify(arg)
            if arg == 0:
                return sympify(0)

    assert tanh(x).series(x, 0, 6) == MyTanh(x).series(x, 0, 6)
