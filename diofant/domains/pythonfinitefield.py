"""Implementation of :class:`PythonFiniteField` class. """

from .finitefield import FiniteField
from .pythonintegerring import PythonIntegerRing


__all__ = ('PythonFiniteField',)


class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_python'

    def __init__(self, mod, symmetric=True):
        return super().__init__(mod, PythonIntegerRing(), symmetric)
