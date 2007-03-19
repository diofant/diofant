from sympy.core import Basic,exp,Symbol,Rational,I,Mul,NCSymbol

"""
This module implements Pauli algebra by subclassing NCSymbol. Only aglebraic
properties of Pauli matrices are used (we don't use the Matrix class).

See the documentation to the class Pauli for examples.

See also:
    http://en.wikipedia.org/wiki/Pauli_matrices
"""

def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

def epsilon(i,j,k):
    if (i,j,k) in [(1,2,3), (2,3,1), (3,1,2)]:
        return 1
    elif (i,j,k) in [(1,3,2), (3,2,1), (2,1,3)]:
        return -1
    else:
        return 0

class Pauli(NCSymbol):
    """
    >>> from sympy import *
    >>> Pauli(1)
    sigma1
    >>> Pauli(1)*Pauli(2)
    i*sigma3
    >>> Pauli(1)*Pauli(1)
    1
    >>> Pauli(3)**4
    1
    >>> Pauli(1)*Pauli(2)*Pauli(3)
    i

    """

    def __init__(self,i):
        if not i in [1,2,3]:
            raise "Invalid Pauli index"
        self.i=i
        NCSymbol.__init__(self, "sigma%d"%i)

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Pauli) and isinstance(y, Pauli):
            j=x.i
            k=y.i
            return delta(j,k) \
                +I*epsilon(j,k,1)*Pauli(1) \
                +I*epsilon(j,k,2)*Pauli(2) \
                +I*epsilon(j,k,3)*Pauli(3)
        return None
