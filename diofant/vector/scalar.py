from ..core import Expr, Integer, Symbol, sympify
from ..printing.pretty.stringpict import prettyForm


class BaseScalar(Expr):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.
    """

    is_commutative = True

    def __new__(cls, name, index, system, pretty_str, latex_str):
        from .coordsysrect import CoordSysCartesian
        index = sympify(index, strict=True)
        system = sympify(system, strict=True)
        obj = super().__new__(cls, Symbol(str(name)), index, system,
                              Symbol(str(pretty_str)), Symbol(str(latex_str)))

        if not isinstance(system, CoordSysCartesian):
            raise TypeError('system should be a CoordSysCartesian')
        if index not in range(3):
            raise ValueError('Invalid index specified.')

        # The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj.name = obj._name = str(name)
        obj._pretty_form = str(pretty_str)
        obj._latex_form = str(latex_str)
        obj._system = system

        return obj

    _diff_wrt = True

    def _eval_derivative(self, s):
        assert self != s  # == case handled in Symbol._eval_derivative
        return Integer(0)

    def _latex(self, printer=None):
        return self._latex_form

    def _pretty(self, printer=None):
        return prettyForm(self._pretty_form)

    @property
    def system(self):
        return self._system

    def __str__(self, printer=None):
        return self._name

    @property
    def free_symbols(self):
        return {self}

    __repr__ = __str__
    _diofantstr = __str__
