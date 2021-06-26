"""Curves in 2-dimensional Euclidean space.

Contains
========
Curve
"""

from ..core import Tuple
from ..core.compatibility import is_sequence
from ..core.sympify import sympify
from .entity import GeometryEntity, GeometrySet
from .point import Point
from .util import _symbol


class Curve(GeometrySet):
    """A curve in space.

    A curve is defined by parametric functions for the coordinates, a
    parameter and the lower and upper bounds for the parameter value.

    Parameters
    ==========

    function : list of functions
    limits : 3-tuple
        Function parameter and lower and upper bounds.

    Raises
    ======

    ValueError
        When `functions` are specified incorrectly.
        When `limits` are specified incorrectly.

    See Also
    ========

    diofant.core.function.Function
    diofant.polys.polyfuncs.interpolate

    Examples
    ========

    >>> C = Curve((sin(t), cos(t)), (t, 0, 2))
    >>> C.functions
    (sin(t), cos(t))
    >>> C.limits
    (t, 0, 2)
    >>> C.parameter
    t
    >>> C = Curve((t, interpolate([1, 4, 9, 16], t)), (t, 0, 1))
    >>> C
    Curve((t, t**2), (t, 0, 1))
    >>> C.subs({t: 4})
    Point(4, 16)
    >>> C.arbitrary_point(a)
    Point(a, a**2)

    """

    def __new__(cls, function, limits):
        fun = sympify(function)
        if not is_sequence(fun) or len(fun) != 2:
            raise ValueError('Function argument should be (x(t), y(t)) '
                             f'but got {function!s}')
        if not is_sequence(limits) or len(limits) != 3:
            raise ValueError('Limit argument should be (t, tmin, tmax) '
                             f'but got {limits!s}')

        return GeometryEntity.__new__(cls, Tuple(*fun), Tuple(*limits))

    def _eval_subs(self, old, new):
        if old == self.parameter:
            return Point(*[f.subs({old: new}) for f in self.functions])

    @property
    def free_symbols(self):
        """
        Return a set of symbols other than the bound symbols used to
        parametrically define the Curve.

        Examples
        ========

        >>> Curve((t, t**2), (t, 0, 2)).free_symbols
        set()
        >>> Curve((t, t**2), (t, a, 2)).free_symbols
        {a}

        """
        free = set()
        for a in self.functions + self.limits[1:]:
            free |= a.free_symbols
        free = free.difference({self.parameter})
        return free

    @property
    def functions(self):
        """The functions specifying the curve.

        Returns
        =======

        functions : list of parameterized coordinate functions.

        See Also
        ========

        parameter

        Examples
        ========

        >>> C = Curve((t, t**2), (t, 0, 2))
        >>> C.functions
        (t, t**2)

        """
        return self.args[0]

    @property
    def parameter(self):
        """The curve function variable.

        Returns
        =======

        parameter : Diofant symbol

        See Also
        ========

        functions

        Examples
        ========

        >>> C = Curve([t, t**2], (t, 0, 2))
        >>> C.parameter
        t

        """
        return self.args[1][0]

    @property
    def limits(self):
        """The limits for the curve.

        Returns
        =======

        limits : tuple
            Contains parameter and lower and upper limits.

        See Also
        ========

        plot_interval

        Examples
        ========

        >>> C = Curve([t, t**3], (t, -2, 2))
        >>> C.limits
        (t, -2, 2)

        """
        return self.args[1]

    def rotate(self, angle=0, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        The default pt is the origin, Point(0, 0).

        Examples
        ========

        >>> Curve((x, x), (x, 0, 1)).rotate(pi/2)
        Curve((-x, x), (x, 0, 1))

        """
        from ..matrices import Matrix, rot_axis3
        pt = -Point(pt or (0, 0))
        rv = self.translate(*pt.args)
        f = list(rv.functions)
        f.append(0)
        f = Matrix(1, 3, f)
        f *= rot_axis3(angle)
        rv = self.func(f[0, :2].tolist()[0], self.limits)
        pt = -pt
        return rv.translate(*pt.args)

    def scale(self, x=1, y=1, pt=None):
        """Override GeometryEntity.scale since Curve is not made up of Points.

        Examples
        ========

        >>> Curve((x, x), (x, 0, 1)).scale(2)
        Curve((2*x, x), (x, 0, 1))

        """
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        fx, fy = self.functions
        return self.func((fx*x, fy*y), self.limits)

    def translate(self, x=0, y=0):
        """Translate the Curve by (x, y).

        Examples
        ========

        >>> Curve((x, x), (x, 0, 1)).translate(1, 2)
        Curve((x + 1, x + 2), (x, 0, 1))

        """
        fx, fy = self.functions
        return self.func((fx + x, fy + y), self.limits)

    def arbitrary_point(self, parameter='t'):
        """
        A parameterized point on the curve.

        Parameters
        ==========

        parameter : str or Symbol, optional
            Default value is 't';
            the Curve's parameter is selected with None or self.parameter
            otherwise the provided symbol is used.

        Returns
        =======

        arbitrary_point : Point

        Raises
        ======

        ValueError
            When `parameter` already appears in the functions.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> from diofant.abc import s
        >>> C = Curve([2*s, s**2], (s, 0, 2))
        >>> C.arbitrary_point()
        Point(2*t, t**2)
        >>> C.arbitrary_point(C.parameter)
        Point(2*s, s**2)
        >>> C.arbitrary_point(None)
        Point(2*s, s**2)
        >>> C.arbitrary_point(a)
        Point(2*a, a**2)

        """
        if parameter is None:
            return Point(*self.functions)

        tnew = _symbol(parameter, self.parameter)
        t = self.parameter
        if (tnew.name != t.name and
                tnew.name in (f.name for f in self.free_symbols)):
            raise ValueError(f'Symbol {tnew.name} already appears in object '
                             'and cannot be used as a parameter.')
        return Point(*[w.subs({t: tnew}) for w in self.functions])

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the curve.

        Parameters
        ==========

        parameter : str or Symbol, optional
            Default value is 't';
            otherwise the provided symbol is used.

        Returns
        =======

        plot_interval : list (plot interval)
            [parameter, lower_bound, upper_bound]

        See Also
        ========

        limits : Returns limits of the parameter interval

        Examples
        ========

        >>> from diofant.abc import s
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval()
        [t, 1, 2]
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval(s)
        [s, 1, 2]

        """
        t = _symbol(parameter, self.parameter)
        return [t] + list(self.limits[1:])
