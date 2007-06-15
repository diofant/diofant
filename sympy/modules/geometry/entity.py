from sympy import Basic, Symbol, Number, Pow, Add, Mul


class GeometryEntity(object):
    """The base class for any geometrical entity."""

    def __init__(self, *args, **kwargs):
        self._args = tuple(args)

    @staticmethod
    def do_intersection(e1, e2):
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
        try:
            if not isinstance(args[0], GeometryEntity):
                args = args[0]
            return list(set(args))
        except:
            return args

    def intersection(self, o):
        """
        Returns the intersection of this entity and another entity, or None if
        there is no intersection.
        """
        raise NotImplementedError()

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(self._args)

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __rmul__(self, a):
        return a.__mul__(self)

    def __rdiv__(self, a):
        return a.__div__(self)

    def __repr__(self):
        return str(self)
