"""Elliptical geometrical entities.

Contains
* Ellipse
* Circle
"""

import random

from ..core import Dummy, Rational, oo, pi
from ..core.logic import fuzzy_bool
from ..core.sympify import sympify
from ..functions import cos, sin, sqrt
from ..logic import false, true
from ..simplify import simplify, trigsimp
from ..solvers import solve
from ..utilities import filldedent
from ..utilities.iterables import uniq
from .entity import GeometryEntity, GeometrySet
from .exceptions import GeometryError
from .line import Line, LinearEntity
from .point import Point
from .util import _symbol, idiff


class Ellipse(GeometrySet):
    """An elliptical GeometryEntity.

    Parameters
    ==========

    center : Point, optional
        Default value is Point(0, 0)
    hradius : number or Diofant expression, optional
    vradius : number or Diofant expression, optional
    eccentricity : number or Diofant expression, optional
        Two of `hradius`, `vradius` and `eccentricity` must be supplied to
        create an Ellipse. The third is derived from the two supplied.

    Raises
    ======

    diofant.geometry.exceptions.GeometryError
        When `hradius`, `vradius` and `eccentricity` are incorrectly supplied
        as parameters.
    TypeError
        When `center` is not a Point.

    See Also
    ========

    Circle

    Notes
    -----
    Constructed from a center and two radii, the first being the horizontal
    radius (along the x-axis) and the second being the vertical radius (along
    the y-axis).

    When symbolic value for hradius and vradius are used, any calculation that
    refers to the foci or the major or minor axis will assume that the ellipse
    has its major radius on the x-axis. If this is not true then a manual
    rotation is necessary.

    Examples
    ========

    >>> e1 = Ellipse(Point(0, 0), 5, 1)
    >>> e1.hradius, e1.vradius
    (5, 1)
    >>> e2 = Ellipse(Point(3, 1), hradius=3, eccentricity=Rational(4, 5))
    >>> e2
    Ellipse(Point(3, 1), 3, 9/5)

    """

    def __new__(
        cls, center=None, hradius=None, vradius=None, eccentricity=None,
            **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center)

        if len(list(filter(None, (hradius, vradius, eccentricity)))) != 2:
            raise ValueError('Exactly two arguments of "hradius", '
                             '"vradius", and "eccentricity" must not be None."')

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(1 - eccentricity**2)
            else:
                assert vradius is None
                vradius = hradius * sqrt(1 - eccentricity**2)

        if hradius == vradius:
            return Circle(center, hradius, **kwargs)

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    @property
    def ambient_dimension(self):
        return 2

    @property
    def center(self):
        """The center of the ellipse.

        Returns
        =======

        center : number

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.center
        Point(0, 0)

        """
        return self.args[0]

    @property
    def hradius(self):
        """The horizontal radius of the ellipse.

        Returns
        =======

        hradius : number

        See Also
        ========

        vradius, major, minor

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.hradius
        3

        """
        return self.args[1]

    @property
    def vradius(self):
        """The vertical radius of the ellipse.

        Returns
        =======

        vradius : number

        See Also
        ========

        hradius, major, minor

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.vradius
        1

        """
        return self.args[2]

    @property
    def minor(self):
        """Shorter axis of the ellipse (if it can be determined) else vradius.

        Returns
        =======

        minor : number or expression

        See Also
        ========

        hradius, vradius, major

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.minor
        1

        >>> Ellipse(p1, a, b).minor
        b
        >>> Ellipse(p1, b, a).minor
        a

        >>> M = m + 1
        >>> Ellipse(p1, m, M).minor
        m

        """
        ab = self.args[1:3]
        if len(ab) == 1:
            return ab[0]
        a, b = ab
        o = a - b < 0
        if o == true:
            return a
        elif o == false:
            return b
        return self.vradius

    @property
    def major(self):
        """Longer axis of the ellipse (if it can be determined) else hradius.

        Returns
        =======

        major : number or expression

        See Also
        ========

        hradius, vradius, minor

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.major
        3

        >>> Ellipse(p1, a, b).major
        a
        >>> Ellipse(p1, b, a).major
        b

        >>> M = m + 1
        >>> Ellipse(p1, m, M).major
        m + 1

        """
        ab = self.args[1:3]
        if len(ab) == 1:
            return ab[0]
        a, b = ab
        o = b - a < 0
        if o == true:
            return a
        elif o == false:
            return b
        return self.hradius

    @property
    def area(self):
        """The area of the ellipse.

        Returns
        =======

        area : number

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.area
        3*pi

        """
        return simplify(pi * self.hradius * self.vradius)

    @property
    def circumference(self):
        """The circumference of the ellipse.

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.circumference
        12*Integral(sqrt((-8*_x**2/9 + 1)/(-_x**2 + 1)), (_x, 0, 1))

        """
        from ..integrals import Integral
        if self.eccentricity == 1:
            return 2*pi*self.hradius
        else:
            x = Dummy('x', extended_real=True)
            return 4*self.major*Integral(
                sqrt((1 - (self.eccentricity*x)**2)/(1 - x**2)), (x, 0, 1))

    @property
    def eccentricity(self):
        """The eccentricity of the ellipse.

        Returns
        =======

        eccentricity : number

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, sqrt(2))
        >>> e1.eccentricity
        sqrt(7)/3

        """
        return self.focus_distance / self.major

    @property
    def periapsis(self):
        """The periapsis of the ellipse.

        The shortest distance between the focus and the contour.

        Returns
        =======

        periapsis : number

        See Also
        ========

        apoapsis : Returns greatest distance between focus and contour

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.periapsis
        -2*sqrt(2) + 3

        """
        return self.major * (1 - self.eccentricity)

    @property
    def apoapsis(self):
        """The apoapsis of the ellipse.

        The greatest distance between the focus and the contour.

        Returns
        =======

        apoapsis : number

        See Also
        ========

        periapsis : Returns shortest distance between foci and contour

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.apoapsis
        2*sqrt(2) + 3

        """
        return self.major * (1 + self.eccentricity)

    @property
    def focus_distance(self):
        """The focal distance of the ellipse.

        The distance between the center and one focus.

        Returns
        =======

        focus_distance : number

        See Also
        ========

        foci

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.focus_distance
        2*sqrt(2)

        """
        return Point.distance(self.center, self.foci[0])

    @property
    def foci(self):
        """The foci of the ellipse.

        Notes
        -----
        The foci can only be calculated if the major/minor axes are known.

        Raises
        ======

        ValueError
            When the major and minor axis cannot be determined.

        See Also
        ========

        diofant.geometry.point.Point
        focus_distance : Returns the distance between focus and center

        Examples
        ========

        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.foci
        (Point(-2*sqrt(2), 0), Point(2*sqrt(2), 0))

        """
        c = self.center
        hr, vr = self.hradius, self.vradius
        if hr == vr:
            return c, c

        # calculate focus distance manually, since focus_distance calls this
        # routine
        fd = sqrt(self.major**2 - self.minor**2)
        if hr == self.minor:
            # foci on the y-axis
            return c + Point(0, -fd), c + Point(0, fd)
        elif hr == self.major:
            # foci on the x-axis
            return c + Point(-fd, 0), c + Point(fd, 0)
        else:
            raise NotImplementedError

    def rotate(self, angle=0, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        Note: since the general ellipse is not supported, only rotations that
        are integer multiples of pi/2 are allowed.

        Examples
        ========

        >>> Ellipse((1, 0), 2, 1).rotate(pi/2)
        Ellipse(Point(0, 1), 1, 2)
        >>> Ellipse((1, 0), 2, 1).rotate(pi)
        Ellipse(Point(-1, 0), 2, 1)

        """
        if self.hradius == self.vradius:
            return self.func(self.center.rotate(angle, pt), self.hradius)
        if (angle/pi).is_integer:
            return super().rotate(angle, pt)
        if (2*angle/pi).is_integer:
            return self.func(self.center.rotate(angle, pt), self.vradius, self.hradius)
        # XXX see https://github.com/sympy/sympy/issues/2815 for general ellipes
        raise NotImplementedError('Only rotations of pi/2 are currently supported for Ellipse.')

    def scale(self, x=1, y=1, pt=None):
        """Override GeometryEntity.scale since it is the major and minor
        axes which must be scaled and they are not GeometryEntities.

        Examples
        ========

        >>> Ellipse((0, 0), 2, 1).scale(2, 4)
        Circle(Point(0, 0), 4)
        >>> Ellipse((0, 0), 2, 1).scale(2)
        Ellipse(Point(0, 0), 4, 1)

        """
        c = self.center
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        h = self.hradius
        v = self.vradius
        return self.func(c.scale(x, y), hradius=h*x, vradius=v*y)

    def reflect(self, line):
        """Override GeometryEntity.reflect since the radius
        is not a GeometryEntity.

        Examples
        ========

        >>> Circle((0, 1), 1).reflect(Line((0, 0), (1, 1)))
        Circle(Point(1, 0), -1)
        >>> Ellipse(Point(3, 4), 1, 3).reflect(Line(Point(0, -4), Point(5, 0)))
        Traceback (most recent call last):
        ...
        NotImplementedError:
        General Ellipse is not supported but the equation of the reflected
        Ellipse is given by the zeros of: f(x, y) = (9*x/41 + 40*y/41 +
        37/41)**2 + (40*x/123 - 3*y/41 - 364/123)**2 - 1

        Notes
        =====

        Until the general ellipse (with no axis parallel to the x-axis) is
        supported a NotImplemented error is raised and the equation whose
        zeros define the rotated ellipse is given.

        """
        from .util import _uniquely_named_symbol

        if line.slope in (0, oo):
            c = self.center
            c = c.reflect(line)
            return self.func(c, -self.hradius, self.vradius)
        else:
            x, y = (_uniquely_named_symbol(name, self, line) for name in 'xy')
            expr = self.equation(x, y)
            p = Point(x, y).reflect(line)
            result = expr.subs(zip((x, y), p.args))
            raise NotImplementedError(filldedent(
                'General Ellipse is not supported but the equation '
                'of the reflected Ellipse is given by the zeros of: ' +
                f'f({x!s}, {y!s}) = {result!s}'))

    def encloses_point(self, p):
        """
        Return True if p is enclosed by (is inside of) self.

        Notes
        -----
        Being on the border of self is considered False.

        Parameters
        ==========

        p : Point

        Returns
        =======

        encloses_point : True, False or None

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> e = Ellipse((0, 0), 3, 2)
        >>> e.encloses_point((0, 0))
        True
        >>> e.encloses_point(e.arbitrary_point(t).subs({t: Rational(1, 2)}))
        False
        >>> e.encloses_point((4, 0))
        False

        """
        p = Point(p)
        if p in self:
            return False

        # if the combined distance from the foci to p (h1 + h2) is less
        # than the combined distance from the foci to the minor axis
        # (which is the same as the major axis length) then p is inside
        # the ellipse
        h1, h2 = (f.distance(p) for f in self.foci)
        test = 2*self.major - (h1 + h2)

        return fuzzy_bool(test.is_positive)

    def tangent_lines(self, p):
        """Tangent lines between `p` and the ellipse.

        If `p` is on the ellipse, returns the tangent line through point `p`.
        Otherwise, returns the tangent line(s) from `p` to the ellipse, or
        None if no tangent line is possible (e.g., `p` inside ellipse).

        Parameters
        ==========

        p : Point

        Returns
        =======

        tangent_lines : list with 1 or 2 Lines

        Raises
        ======

        NotImplementedError
            Can only find tangent lines for a point, `p`, on the ellipse.

        See Also
        ========

        diofant.geometry.point.Point, diofant.geometry.line.Line

        Examples
        ========

        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.tangent_lines(Point(3, 0))
        [Line(Point(3, 0), Point(3, -12))]

        """
        p = Point(p)
        if self.encloses_point(p):
            return []

        if p in self:
            delta = self.center - p
            rise = (self.vradius ** 2)*delta.x
            run = -(self.hradius ** 2)*delta.y
            p2 = Point(simplify(p.x + run),
                       simplify(p.y + rise))
            return [Line(p, p2)]
        else:
            # p is outside the ellipse or we can't tell. In case of the
            # latter, the solutions returned will only be valid if
            # the point is not inside the ellipse; if it is, nan will result.
            x, y = Dummy('x'), Dummy('y')
            eq = self.equation(x, y)
            dydx = idiff(eq, y, x)
            slope = Line(p, Point(x, y)).slope
            tangent_points = solve([slope - dydx, eq], [x, y])

            # handle horizontal and vertical tangent lines
            if len(tangent_points) == 1:
                assert (tangent_points[0][x] == p.x or
                        tangent_points[0][y] == p.y)
                return [Line(p, p + Point(1, 0)), Line(p, p + Point(0, 1))]

            # others
            p1 = (tangent_points[0][x], tangent_points[0][y])
            p2 = (tangent_points[1][x], tangent_points[1][y])
            return [Line(p, p1), Line(p, p2)]

    def is_tangent(self, o):
        """Is `o` tangent to the ellipse?

        Parameters
        ==========

        o : GeometryEntity
            An Ellipse, LinearEntity or Polygon

        Raises
        ======

        NotImplementedError
            When the wrong type of argument is supplied.

        Returns
        =======

        is_tangent: boolean
            True if o is tangent to the ellipse, False otherwise.

        See Also
        ========

        tangent_lines

        Examples
        ========

        >>> p0, p1, p2 = Point(0, 0), Point(3, 0), Point(3, 3)
        >>> e1 = Ellipse(p0, 3, 2)
        >>> l1 = Line(p1, p2)
        >>> e1.is_tangent(l1)
        True

        """
        from .polygon import Polygon
        inter = None
        if isinstance(o, Ellipse):
            inter = self.intersection(o)
            if isinstance(inter, Ellipse):
                return False
            return (inter is not None and len(inter) == 1
                    and isinstance(inter[0], Point))
        elif isinstance(o, LinearEntity):
            inter = self._do_line_intersection(o)
            if inter is not None and len(inter) == 1:
                return inter[0] in o
            else:
                return False
        elif isinstance(o, Polygon):
            c = 0
            for seg in o.sides:
                inter = self._do_line_intersection(seg)
                c += len([True for point in inter if point in seg])
            return c == 1
        else:
            raise NotImplementedError('Unknown argument type')

    def normal_lines(self, p, prec=None):
        """Normal lines between `p` and the ellipse.

        Parameters
        ==========

        p : Point

        Returns
        =======

        normal_lines : list with 1, 2 or 4 Lines

        Examples
        ========

        >>> e = Ellipse((0, 0), 2, 3)
        >>> c = e.center
        >>> e.normal_lines(c + Point(1, 0))
        [Line(Point(0, 0), Point(1, 0))]
        >>> e.normal_lines(c)
        [Line(Point(0, 0), Point(0, 1)), Line(Point(0, 0), Point(1, 0))]

        Off-axis points require the solution of a quartic equation. This
        often leads to very large expressions that may be of little practical
        use. An approximate solution of `prec` digits can be obtained by
        passing in the desired value:

        >>> e.normal_lines((3, 3), prec=2)
        [Line(Point(-38/47, -85/31), Point(9/47, -21/17)),
        Line(Point(19/13, -43/21), Point(32/13, -8/3))]

        Whereas the above solution has an operation count of 12, the exact
        solution has an operation count of 2020.

        """
        p = Point(p)

        # XXX change True to something like self.angle == 0 if the arbitrarily
        # rotated ellipse is introduced.
        # https://github.com/sympy/sympy/issues/2815)
        if True:  # pylint: disable=using-constant-test
            rv = []
            if p.x == self.center.x:
                rv.append(Line(self.center, slope=oo))
            if p.y == self.center.y:
                rv.append(Line(self.center, slope=0))
            if rv:
                # at these special orientations of p either 1 or 2 normals
                # exist and we are done
                return rv

        # find the 4 normal points and construct lines through them with
        # the corresponding slope
        x, y = Dummy('x', extended_real=True), Dummy('y', extended_real=True)
        eq = self.equation(x, y)
        dydx = idiff(eq, y, x)
        norm = -1/dydx
        slope = Line(p, (x, y)).slope
        seq = slope - norm
        yis = solve(seq, y)[0][y]
        xeq = eq.subs({y: yis}).as_numer_denom()[0].expand()
        if len(xeq.free_symbols) == 1:
            xsol = xeq.as_poly(x).real_roots()
            points = [Point(i, solve(eq.subs({x: i}), y)[0][y])
                      for i in xsol]
        else:
            raise NotImplementedError(
                'intersections for the general ellipse are not supported')
        slopes = [norm.subs(zip((x, y), pt.args)) for pt in points]
        if prec is not None:
            points = [pt.evalf(prec) for pt in points]
            slopes = [i if i.is_Number and i.is_finite is not False else i.evalf(prec) for i in slopes]
        return [Line(pt, slope=s) for pt, s in zip(points, slopes)]

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the ellipse.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

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

        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.arbitrary_point()
        Point(3*cos(t), 2*sin(t))

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError(filldedent(f'Symbol {t.name} already appears in object '
                                        'and cannot be used as a parameter.'))
        return Point(self.center.x + self.hradius*cos(t),
                     self.center.y + self.vradius*sin(t))

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Ellipse.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.plot_interval()
        [t, -pi, pi]

        """
        t = _symbol(parameter)
        return [t, -pi, pi]

    def random_point(self, seed=None):
        """A random point on the ellipse.

        Returns
        =======

        point : Point

        See Also
        ========

        diofant.geometry.point.Point
        arbitrary_point : Returns parameterized point on ellipse

        Examples
        ========

        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.random_point()  # gives some random point
        Point(...)
        >>> p1 = e1.random_point(seed=0)
        >>> p1.evalf(2)
        Point(2.1, 1.4)

        The random_point method assures that the point will test as being
        in the ellipse:

        >>> p1 in e1
        True

        Notes
        =====

        A random point may not appear to be on the ellipse, ie, `p in e` may
        return False. This is because the coordinates of the point will be
        floating point values, and when these values are substituted into the
        equation for the ellipse the result may not be zero because of floating
        point rounding error.

        An arbitrary_point with a random value of t substituted into it may
        not test as being on the ellipse because the expression tested that
        a point is on the ellipse doesn't simplify to zero and doesn't evaluate
        exactly to zero:

        >>> e1.arbitrary_point(t)
        Point(3*cos(t), 2*sin(t))
        >>> p2 = _.subs({t: 0.1})
        >>> p2 in e1
        False

        Note that arbitrary_point routine does not take this approach. A value
        for cos(t) and sin(t) (not t) is substituted into the arbitrary point.
        There is a small chance that this will give a point that will not
        test as being in the ellipse, so the process is repeated (up to 10
        times) until a valid point is obtained.

        """
        from ..functions import cos, sin
        t = _symbol('t')
        x, y = self.arbitrary_point(t).args
        # get a random value in [-1, 1) corresponding to cos(t)
        # and confirm that it will test as being in the ellipse
        if seed is not None:
            rng = random.Random(seed)
        else:
            rng = random
        # simplify this now or else the Float will turn s into a Float
        c = 2*Rational(rng.random()) - 1
        s = sqrt(1 - c**2)
        return Point(x.subs({cos(t): c}), y.subs({sin(t): s}))

    def equation(self, x='x', y='y'):
        """The equation of the ellipse.

        Parameters
        ==========

        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.

        Returns
        =======

        equation : diofant expression

        See Also
        ========

        arbitrary_point : Returns parameterized point on ellipse

        Examples
        ========

        >>> e1 = Ellipse(Point(1, 0), 3, 2)
        >>> e1.equation()
        y**2/4 + (x/3 - 1/3)**2 - 1

        """
        x = _symbol(x)
        y = _symbol(y)
        t1 = ((x - self.center.x) / self.hradius)**2
        t2 = ((y - self.center.y) / self.vradius)**2
        return t1 + t2 - 1

    def _do_line_intersection(self, o):
        """
        Find the intersection of a LinearEntity and the ellipse.

        All LinearEntities are treated as a line and filtered at
        the end to see that they lie in o.

        """
        hr_sq = self.hradius ** 2
        vr_sq = self.vradius ** 2
        lp = o.points

        ldir = lp[1] - lp[0]
        diff = lp[0] - self.center
        mdir = Point(ldir.x/hr_sq, ldir.y/vr_sq)
        mdiff = Point(diff.x/hr_sq, diff.y/vr_sq)

        a = ldir.dot(mdir)
        b = ldir.dot(mdiff)
        c = diff.dot(mdiff) - 1
        det = simplify(b*b - a*c)

        result = []
        if det == 0:
            t = -b / a
            result.append(lp[0] + (lp[1] - lp[0]) * t)
        # Definite and potential symbolic intersections are allowed.
        elif det.is_positive is not False:
            root = sqrt(det)
            t_a = (-b - root) / a
            t_b = (-b + root) / a
            result.append(lp[0] + (lp[1] - lp[0]) * t_a)
            result.append(lp[0] + (lp[1] - lp[0]) * t_b)

        return [r for r in result if r in o]

    def _do_ellipse_intersection(self, o):
        """The intersection of an ellipse with another ellipse or a circle.

        Private helper method for `intersection`.

        """
        x = Dummy('x', extended_real=True)
        y = Dummy('y', extended_real=True)
        seq = self.equation(x, y)
        oeq = o.equation(x, y)
        result = [(s[x], s[y]) for s in solve([seq, oeq], [x, y])]
        return [Point(*r) for r in list(uniq(result))]

    def intersection(self, o):
        """The intersection of this ellipse and another geometrical entity
        `o`.

        Parameters
        ==========

        o : GeometryEntity

        Returns
        =======

        intersection : list of GeometryEntity objects

        Notes
        -----
        Currently supports intersections with Point, Line, Segment, Ray,
        Circle and Ellipse types.

        See Also
        ========

        diofant.geometry.entity.GeometryEntity

        Examples
        ========

        >>> e = Ellipse(Point(0, 0), 5, 7)
        >>> e.intersection(Point(0, 0))
        []
        >>> e.intersection(Point(5, 0))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(0, 0), Point(0, 1)))
        [Point(0, -7), Point(0, 7)]
        >>> e.intersection(Line(Point(5, 0), Point(5, 1)))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(6, 0), Point(6, 1)))
        []
        >>> e = Ellipse(Point(-1, 0), 4, 3)
        >>> e.intersection(Ellipse(Point(1, 0), 4, 3))
        [Point(0, -3*sqrt(15)/4), Point(0, 3*sqrt(15)/4)]
        >>> e.intersection(Ellipse(Point(5, 0), 4, 3))
        [Point(2, -3*sqrt(7)/4), Point(2, 3*sqrt(7)/4)]
        >>> e.intersection(Ellipse(Point(100500, 0), 4, 3))
        []
        >>> e.intersection(Ellipse(Point(0, 0), 3, 4))
        [Point(-363/175, -48*sqrt(111)/175), Point(-363/175, 48*sqrt(111)/175), Point(3, 0)]

        >>> e.intersection(Ellipse(Point(-1, 0), 3, 4))
        [Point(-17/5, -12/5), Point(-17/5, 12/5), Point(7/5, -12/5), Point(7/5, 12/5)]

        """
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, LinearEntity):
            # LinearEntity may be a ray/segment, so check the points
            # of intersection for coincidence first
            return self._do_line_intersection(o)

        elif isinstance(o, Circle):
            return self._do_ellipse_intersection(o)

        elif isinstance(o, Ellipse):
            if o == self:
                return self
            else:
                return self._do_ellipse_intersection(o)

        return o.intersection(self)

    def evolute(self, x='x', y='y'):
        """The equation of evolute of the ellipse.

        Parameters
        ==========

        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.

        Returns
        =======

        equation : diofant expression

        Examples
        ========

        >>> e1 = Ellipse(Point(1, 0), 3, 2)
        >>> e1.evolute()
        2**(2/3)*y**(2/3) + (3*x - 3)**(2/3) - 5**(2/3)

        """
        if len(self.args) != 3:
            raise NotImplementedError('Evolute of arbitrary Ellipse is not supported.')
        x = _symbol(x)
        y = _symbol(y)
        t1 = (self.hradius*(x - self.center.x))**Rational(2, 3)
        t2 = (self.vradius*(y - self.center.y))**Rational(2, 3)
        return t1 + t2 - (self.hradius**2 - self.vradius**2)**Rational(2, 3)

    def __eq__(self, o):
        """Is the other GeometryEntity the same as this ellipse?"""
        return isinstance(o, Ellipse) and (self.center == o.center and
                                           self.hradius == o.hradius and
                                           self.vradius == o.vradius)

    __hash__ = GeometrySet.__hash__

    def __contains__(self, o):
        if isinstance(o, Point):
            x = Dummy('x', extended_real=True)
            y = Dummy('y', extended_real=True)

            res = self.equation(x, y).subs({x: o.x, y: o.y})
            return trigsimp(simplify(res)) == 0
        elif isinstance(o, Ellipse):
            return self == o
        return False


class Circle(Ellipse):
    """A circle in space.

    Constructed simply from a center and a radius, or from three
    non-collinear points.

    Parameters
    ==========

    center : Point
    radius : number or diofant expression
    points : sequence of three Points

    Raises
    ======

    diofant.geometry.exceptions.GeometryError
        When trying to construct circle from three collinear points.
        When trying to construct circle from incorrect parameters.

    See Also
    ========

    Ellipse, diofant.geometry.point.Point

    Examples
    ========

    >>> # a circle constructed from a center and radius
    >>> c1 = Circle(Point(0, 0), 5)
    >>> c1.hradius, c1.vradius, c1.radius
    (5, 5, 5)

    >>> # a circle constructed from three points
    >>> c2 = Circle(Point(0, 0), Point(1, 1), Point(1, 0))
    >>> c2.hradius, c2.vradius, c2.radius, c2.center
    (sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, Point(1/2, 1/2))

    """

    def __new__(cls, *args, **kwargs):
        c, r = None, None
        if len(args) == 3:
            args = [Point(a) for a in args]
            if Point.is_collinear(*args):
                raise GeometryError(
                    'Cannot construct a circle from three collinear points')
            from .polygon import Triangle
            t = Triangle(*args)
            c = t.circumcenter
            r = t.circumradius
        elif len(args) == 2:
            # Assume (center, radius) pair
            c = Point(args[0])
            r = sympify(args[1])

        if not (c is None or r is None):
            return GeometryEntity.__new__(cls, c, r, **kwargs)

        raise GeometryError('Circle.__new__ received unknown arguments')

    @property
    def radius(self):
        """The radius of the circle.

        Returns
        =======

        radius : number or diofant expression

        See Also
        ========

        Ellipse.major, Ellipse.minor, Ellipse.hradius, Ellipse.vradius

        Examples
        ========

        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.radius
        6

        """
        return self.args[1]

    @property
    def vradius(self):
        """
        This Ellipse property is an alias for the Circle's radius.

        Whereas hradius, major and minor can use Ellipse's conventions,
        the vradius does not exist for a circle. It is always a positive
        value in order that the Circle, like Polygons, will have an
        area that can be positive or negative as determined by the sign
        of the hradius.

        Examples
        ========

        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.vradius
        6

        """
        return abs(self.radius)

    @property
    def circumference(self):
        """The circumference of the circle.

        Returns
        =======

        circumference : number or Diofant expression

        Examples
        ========

        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.circumference
        12*pi

        """
        return 2 * pi * self.radius

    def equation(self, x='x', y='y'):
        """The equation of the circle.

        Parameters
        ==========

        x : str or Symbol, optional
            Default value is 'x'.
        y : str or Symbol, optional
            Default value is 'y'.

        Returns
        =======

        equation : Diofant expression

        Examples
        ========

        >>> c1 = Circle(Point(0, 0), 5)
        >>> c1.equation()
        x**2 + y**2 - 25

        """
        x = _symbol(x)
        y = _symbol(y)
        t1 = (x - self.center.x)**2
        t2 = (y - self.center.y)**2
        return t1 + t2 - self.major**2

    def intersection(self, o):
        """The intersection of this circle with another geometrical entity.

        Parameters
        ==========

        o : GeometryEntity

        Returns
        =======

        intersection : list of GeometryEntities

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(5, 5), Point(6, 0)
        >>> p4 = Point(5, 0)
        >>> c1 = Circle(p1, 5)
        >>> c1.intersection(p2)
        []
        >>> c1.intersection(p4)
        [Point(5, 0)]
        >>> c1.intersection(Ray(p1, p2))
        [Point(5*sqrt(2)/2, 5*sqrt(2)/2)]
        >>> c1.intersection(Line(p2, p3))
        []

        """
        if isinstance(o, Circle):
            if o.center == self.center:
                if o.radius == self.radius:
                    return o
                return []
            dx, dy = (o.center - self.center).args
            d = sqrt(simplify(dy**2 + dx**2))
            R = o.radius + self.radius
            if d > R or d < abs(self.radius - o.radius):
                return []

            a = simplify((self.radius**2 - o.radius**2 + d**2) / (2*d))

            x2 = self.center.x + (dx * a/d)
            y2 = self.center.y + (dy * a/d)

            h = sqrt(simplify(self.radius**2 - a**2))
            rx = -dy * (h/d)
            ry = dx * (h/d)

            xi_1 = simplify(x2 + rx)
            xi_2 = simplify(x2 - rx)
            yi_1 = simplify(y2 + ry)
            yi_2 = simplify(y2 - ry)

            ret = [Point(xi_1, yi_1)]
            if xi_1 != xi_2 or yi_1 != yi_2:
                ret.append(Point(xi_2, yi_2))
            return ret

        return Ellipse.intersection(self, o)

    def scale(self, x=1, y=1, pt=None):
        """Override GeometryEntity.scale since the radius
        is not a GeometryEntity.

        Examples
        ========

        >>> Circle((0, 0), 1).scale(2, 2)
        Circle(Point(0, 0), 2)
        >>> Circle((0, 0), 1).scale(2, 4)
        Ellipse(Point(0, 0), 2, 4)

        """
        c = self.center
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        c = c.scale(x, y)
        x, y = map(abs, (x, y))
        if x == y:
            return self.func(c, x*self.radius)
        h = v = self.radius
        return Ellipse(c, hradius=h*x, vradius=v*y)

    def reflect(self, line):
        """Override GeometryEntity.reflect since the radius
        is not a GeometryEntity.

        Examples
        ========

        >>> Circle((0, 1), 1).reflect(Line((0, 0), (1, 1)))
        Circle(Point(1, 0), -1)

        """
        c = self.center
        c = c.reflect(line)
        return self.func(c, -self.radius)
