"""Implementation of quotient ring elements."""

import numbers

from ..polys.polyerrors import CoercionFailedError
from .domainelement import DomainElement


class QuotientRingElement(DomainElement):
    """A class representing an element of the quotent ring."""

    @property
    def parent(self):
        return self._parent

    def __init__(self, rep):
        """Initialize self."""
        if isinstance(rep, self.__class__):
            self.rep = rep.rep % self.mod
        else:
            self.rep = self.domain.convert(rep) % self.mod

    def __reduce__(self):
        return self.parent.__call__, (self.rep,)

    def __hash__(self):
        return hash((self.rep, self.mod))

    def __int__(self):
        return int(self.rep)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.rep)

    def __add__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self.__class__(self.rep + other.rep)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self.__class__(self.rep - other.rep)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self.__class__(self.rep * other.rep)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self * other**-1

    def __rtruediv__(self, other):
        return (self**-1).__mul__(other)

    def __mod__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self.__class__(self.rep % other.rep)

    def __rmod__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return other.__mod__(self)

    def __pow__(self, exp):
        if not isinstance(exp, numbers.Integral):
            raise TypeError(f'Integer exponent expected, got {type(exp)}')
        if exp < 0:
            rep, exp = self.domain.invert(self.rep, self.mod), -exp
        else:
            rep = self.rep
        return self.__class__(rep**exp)

    def __eq__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailedError:
            return NotImplemented
        return self.rep == other.rep

    def __bool__(self):
        return bool(self.rep)
