"""Implementation of :class:`FractionField` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.polyclasses import DMF
from sympy.polys.polyerrors import GeneratorsNeeded
from sympy.polys.polyutils import (dict_from_basic, basic_from_dict,
                                   _dict_reorder)
from sympy.utilities import public


@public
class FractionField(Field, CharacteristicZero, CompositeDomain):
    """A class for representing rational function fields. """

    dtype = DMF
    is_FractionField = is_Frac = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def __init__(self, dom, *gens):
        if not gens:
            raise GeneratorsNeeded("generators not specified")

        lev = len(gens) - 1
        self.ngens = len(gens)

        self.zero = self.dtype.zero(lev, dom, ring=self)
        self.one = self.dtype.one(lev, dom, ring=self)

        self.domain = self.dom = dom
        self.symbols = self.gens = gens

    def new(self, element):
        return self.dtype(element, self.dom, len(self.gens) - 1, ring=self)

    def __str__(self):
        return str(self.dom) + '(' + ','.join(map(str, self.gens)) + ')'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.gens))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, FractionField) and \
            self.dtype == other.dtype and self.dom == other.dom and self.gens == other.gens

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        return (basic_from_dict(a.numer().to_sympy_dict(), *self.gens) /
                basic_from_dict(a.denom().to_sympy_dict(), *self.gens))

    def from_sympy(self, a):
        """Convert SymPy's expression to ``dtype``. """
        p, q = a.as_numer_denom()

        num, _ = dict_from_basic(p, gens=self.gens)
        den, _ = dict_from_basic(q, gens=self.gens)

        for k, v in num.items():
            num[k] = self.dom.from_sympy(v)

        for k, v in den.items():
            den[k] = self.dom.from_sympy(v)

        return self((num, den)).cancel()

    def from_FractionField(self, a, K0):
        """
        Convert a fraction field element to another fraction field.

        Examples
        ========

        >>> from sympy.polys.polyclasses import DMF
        >>> from sympy.polys.domains import ZZ, QQ
        >>> from sympy.abc import x

        >>> f = DMF(([ZZ(1), ZZ(2)], [ZZ(1), ZZ(1)]), ZZ)

        >>> QQx = QQ.old_frac_field(x)
        >>> ZZx = ZZ.old_frac_field(x)

        >>> QQx.from_FractionField(f, ZZx)
        (x + 2)/(x + 1)
        """
        if self.gens == K0.gens:
            if self.dom == K0.dom:
                return a
            else:
                return self((a.numer().convert(self.dom).rep,
                             a.denom().convert(self.dom).rep))
        else:  # pragma: no cover
            raise NotImplementedError
