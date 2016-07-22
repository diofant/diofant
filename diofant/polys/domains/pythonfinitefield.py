"""Implementation of :class:`PythonFiniteField` class. """

from diofant.polys.domains.finitefield import FiniteField
from diofant.polys.domains.pythonintegerring import PythonIntegerRing
from diofant.utilities import public


@public
class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_python'

    def __init__(self, mod, symmetric=True):
        return super(PythonFiniteField, self).__init__(mod, PythonIntegerRing(), symmetric)
