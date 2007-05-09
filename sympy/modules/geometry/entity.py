from sympy.core import Basic, Symbol, Number, Pow, Add, Mul


class GeometryEntity(object):
    """The base class for any geometrical entity."""

    def __init__(self, *args, **kwargs):
        pass

    @staticmethod
    def intersection(e1, e2):
        """Determine the intersection between two geometrical entities"""
        try:
            return e1.intersection(e2)
        except AttributeError, NotImplementedError:
            pass

        try:
            return e2.intersection(e1)
        except AttributeError, NotImplementedError:
            n1,n2 = type(e1).__name__, type(e2).__name__
            raise NotImplementedError("Unable to determine intersection between '%s' and '%s'" % (n1, n2))

    @staticmethod
    def _normalize_args(args):
        """Removes duplicates of points."""
        if not isinstance(args[0], GeometryEntity):
            args = args[0]
        return list(set(args))

    @staticmethod
    def _is_symnum(arg):
        """Check to see if an expression is a symbol, number, or some combination."""
        return isinstance(arg, (Symbol,Number,Mul,Add,Pow))

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __rmul__(self, a):
        return a.__mul__(self)

    def __rdiv__(self, a):
        return a.__div__(self)