"""Implementation of :class:`PythonRationalField` class. """

from .groundtypes import PythonRational, python_factorial
from .rationalfield import RationalField


__all__ = ('PythonRationalField',)


class PythonRationalField(RationalField):
    """Rational field based on Python's rationals. """

    dtype = PythonRational
    zero = dtype(0)
    one = dtype(1)

    def __init__(self):
        pass

    @property
    def ring(self):
        """Returns ring associated with ``self``. """
        from . import PythonIntegerRing
        return PythonIntegerRing()

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(python_factorial(int(a)))
