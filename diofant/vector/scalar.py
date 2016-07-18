from diofant.core.symbol import Symbol
from diofant.core.numbers import Integer
from diofant.printing.pretty.stringpict import prettyForm


class BaseScalar(Symbol):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, name, index, system, pretty_str, latex_str):
        name = str(name)
        pretty_str = str(pretty_str)
        latex_str = str(latex_str)
        from diofant.vector.coordsysrect import CoordSysCartesian
        obj = super(BaseScalar, cls).__new__(cls, name)
        if not isinstance(system, CoordSysCartesian):
            raise TypeError("system should be a CoordSysCartesian")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        # The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj._name = name
        obj._pretty_form = pretty_str
        obj._latex_form = latex_str
        obj._system = system

        # Change the args for the object
        obj._args = tuple([Symbol(name), Integer(index), system,
                           Symbol(pretty_str), Symbol(latex_str)])

        return obj

    def _latex(self, printer=None):
        return self._latex_form

    def _pretty(self, printer=None):
        return prettyForm(self._pretty_form)

    @property
    def system(self):
        return self._system

    def __eq__(self, other):
        # Check if the other object is a BaseScalar of same index
        # and coordinate system
        if isinstance(other, BaseScalar):
            if other._id == self._id:
                return True
        return False

    def __hash__(self):
        return self._id.__hash__()

    def __str__(self, printer=None):
        return self._name

    @property
    def free_symbols(self):
        return {self}

    __repr__ = __str__
    _diofantstr = __str__
