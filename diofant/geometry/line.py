"""Line-like geometrical entities.

Contains
========
LinearEntity
Line
Ray
Segment
"""

from ..core import Dummy, Eq, Integer, factor_terms, oo, pi
from ..core.compatibility import is_sequence
from ..core.sympify import sympify
from ..functions import Piecewise, acos, sqrt, tan
from ..functions.elementary.trigonometric import _pi_coeff as pi_coeff
from ..logic import And, false, true
from ..simplify import simplify
from ..solvers import solve
from .entity import GeometryEntity, GeometrySet
from .exceptions import GeometryError
from .point import Point
from .util import _symbol


# TODO: this should be placed elsewhere and reused in other modules


class Undecidable(ValueError):
    """Raised when can't decide on relation."""


class LinearEntity(GeometrySet):
    """A base class for all linear entities (line, ray and segment)
    in a 2-dimensional Euclidean space.

    Notes
    =====

    This is an abstract class and is not meant to be instantiated.

    See Also
    ========

    diofant.geometry.entity.GeometryEntity

    """

    def __new__(cls, p1, p2, **kwargs):
        p1 = Point(p1)
        p2 = Point(p2)
        if p1 == p2:
            # sometimes we return a single point if we are not given two unique
            # points. This is done in the specific subclass
            raise ValueError(
                f'{cls.__name__}.__new__ requires two unique Points.')

        return GeometryEntity.__new__(cls, p1, p2, **kwargs)

    @property
    def p1(self):
        """The first defining point of a linear entity.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l = Line(p1, p2)
        >>> l.p1
        Point(0, 0)

        """
        return self.args[0]

    @property
    def p2(self):
        """The second defining point of a linear entity.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l = Line(p1, p2)
        >>> l.p2
        Point(5, 3)

        """
        return self.args[1]

    @property
    def coefficients(self):
        """The coefficients (`a`, `b`, `c`) for `ax + by + c = 0`.

        See Also
        ========

        diofant.geometry.line.Line.equation

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l = Line(p1, p2)
        >>> l.coefficients
        (-3, 5, 0)

        >>> p3 = Point(x, y)
        >>> l2 = Line(p1, p3)
        >>> l2.coefficients
        (-y, x, 0)

        """
        p1, p2 = self.points
        if p1.x == p2.x:
            return Integer(1), Integer(0), -p1.x
        elif p1.y == p2.y:
            return Integer(0), Integer(1), -p1.y
        return tuple(simplify(i)
                     for i in (self.p1.y - self.p2.y,
                               self.p2.x - self.p1.x,
                               self.p1.x*self.p2.y - self.p1.y*self.p2.x))

    @staticmethod
    def are_concurrent(*lines):
        """Is a sequence of linear entities concurrent?

        Two or more linear entities are concurrent if they all
        intersect at a single point.

        Parameters
        ==========

        lines : a sequence of linear entities.

        Returns
        =======

        True : if the set of linear entities are concurrent,
        False : otherwise.

        Notes
        =====

        Simply take the first two lines and find their intersection.
        If there is no intersection, then the first two lines were
        parallel and had no intersection so concurrency is impossible
        amongst the whole set. Otherwise, check to see if the
        intersection point of the first two lines is a member on
        the rest of the lines. If so, the lines are concurrent.

        See Also
        ========

        diofant.geometry.util.intersection

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(3, 5)
        >>> p3, p4 = Point(-2, -2), Point(0, 2)
        >>> l1, l2, l3 = Line(p1, p2), Line(p1, p3), Line(p1, p4)
        >>> Line.are_concurrent(l1, l2, l3)
        True

        >>> l4 = Line(p2, p3)
        >>> Line.are_concurrent(l2, l3, l4)
        False

        """
        # Concurrency requires intersection at a single point; One linear
        # entity cannot be concurrent.
        if len(lines) <= 1:
            return False

        # Get the intersection (if parallel)
        p = lines[0].intersection(lines[1])

        # Make sure the intersection is on every linear entity
        for line in lines[2:]:
            if p[0] not in line:
                return False
        return True

    def is_parallel(self, other):
        """Are two linear entities parallel?

        Parameters
        ==========

        self : LinearEntity
        other : LinearEntity

        Returns
        =======

        True : if self and other are parallel,
        False : otherwise.

        See Also
        ========

        coefficients

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4 = Point(3, 4), Point(6, 7)
        >>> l1, l2 = Line(p1, p2), Line(p3, p4)
        >>> Line.is_parallel(l1, l2)
        True

        >>> p5 = Point(6, 6)
        >>> l3 = Line(p3, p5)
        >>> Line.is_parallel(l1, l3)
        False

        """
        a1, b1, _ = self.coefficients
        a2, b2, _ = other.coefficients
        return bool(simplify(a1*b2 - b1*a2) == 0)

    def is_perpendicular(self, other):
        """Are two linear entities perpendicular?

        Parameters
        ==========

        self : LinearEntity
        other : LinearEntity

        Returns
        =======

        True : if self and other are perpendicular,
        False : otherwise.

        See Also
        ========

        coefficients

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(-1, 1)
        >>> l1, l2 = Line(p1, p2), Line(p1, p3)
        >>> l1.is_perpendicular(l2)
        True

        >>> p4 = Point(5, 3)
        >>> l3 = Line(p1, p4)
        >>> l1.is_perpendicular(l3)
        False

        """
        a1, b1, _ = self.coefficients
        a2, b2, _ = other.coefficients
        return bool(simplify(a1*a2 + b1*b2) == 0)

    def angle_between(self, other):
        """The angle formed between the two linear entities.

        Parameters
        ==========

        self : LinearEntity
        other : LinearEntity

        Returns
        =======

        angle : angle in radians

        Notes
        =====

        From the dot product of vectors v1 and v2 it is known that:

            ``dot(v1, v2) = |v1|*|v2|*cos(A)``

        where A is the angle formed between the two vectors. We can
        get the directional vectors of the two lines and readily
        find the angle between the two using the above formula.

        See Also
        ========

        is_perpendicular

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(0, 4), Point(2, 0)
        >>> l1, l2 = Line(p1, p2), Line(p1, p3)
        >>> l1.angle_between(l2)
        pi/2

        """
        v1 = self.p2 - self.p1
        v2 = other.p2 - other.p1
        return acos(v1.dot(v2)/(abs(v1)*abs(v2)))

    def parallel_line(self, p):
        """Create a new Line parallel to this linear entity which passes
        through the point `p`.

        Parameters
        ==========

        p : diofant.geometry.point.Point

        Returns
        =======

        line : Line

        See Also
        ========

        is_parallel

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(2, 3), Point(-2, 2)
        >>> l1 = Line(p1, p2)
        >>> l2 = l1.parallel_line(p3)
        >>> p3 in l2
        True
        >>> l1.is_parallel(l2)
        True

        """
        d = self.p1 - self.p2
        p = Point(p)
        return Line(p, p + d)

    def perpendicular_line(self, p):
        """Create a new Line perpendicular to this linear entity which passes
        through the point `p`.

        Parameters
        ==========

        p : diofant.geometry.point.Point

        Returns
        =======

        line : Line

        See Also
        ========

        is_perpendicular, perpendicular_segment

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(2, 3), Point(-2, 2)
        >>> l1 = Line(p1, p2)
        >>> l2 = l1.perpendicular_line(p3)
        >>> p3 in l2
        True
        >>> l1.is_perpendicular(l2)
        True

        """
        p = Point(p)
        d1, d2 = (self.p1 - self.p2).args
        if d2 == 0:  # If a horizontal line
            if p.y == self.p1.y:  # if p is on this linear entity
                return Line(p, p + Point(0, 1))
            else:
                p2 = Point(p.x, self.p1.y)
                return Line(p, p2)
        else:
            p2 = Point(p.x - d2, p.y + d1)
            return Line(p, p2)

    def perpendicular_segment(self, p):
        """Create a perpendicular line segment from `p` to this line.

        The enpoints of the segment are ``p`` and the closest point in
        the line containing self. (If self is not a line, the point might
        not be in self.)

        Parameters
        ==========

        p : diofant.geometry.point.Point

        Returns
        =======

        segment : Segment

        Notes
        =====

        Returns `p` itself if `p` is on this linear entity.

        See Also
        ========

        perpendicular_line

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(0, 2)
        >>> l1 = Line(p1, p2)
        >>> s1 = l1.perpendicular_segment(p3)
        >>> l1.is_perpendicular(s1)
        True
        >>> p3 in s1
        True
        >>> l1.perpendicular_segment(Point(4, 0))
        Segment(Point(2, 2), Point(4, 0))

        """
        p = Point(p)
        if p in self:
            return p
        a, b, c = self.coefficients
        if a == 0:  # horizontal
            p2 = Point(p.x, self.p1.y)
        elif b == 0:  # vertical
            p2 = Point(self.p1.x, p.y)
        else:
            # ax + by + c = 0
            y = (-c - a*p.x)/b
            m = self.slope
            d2 = 1 + m**2
            H = p.y - y
            dx = m*H/d2
            dy = m*dx
            p2 = (p.x + dx, y + dy)
        return Segment(p, p2)

    @property
    def length(self):
        """
        The length of the line.

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(3, 5)
        >>> l1 = Line(p1, p2)
        >>> l1.length
        oo

        """
        return oo

    @property
    def slope(self):
        """The slope of this linear entity, or infinity if vertical.

        Returns
        =======

        slope : Expr

        See Also
        ========

        coefficients

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(3, 5)
        >>> l1 = Line(p1, p2)
        >>> l1.slope
        5/3

        >>> p3 = Point(0, 4)
        >>> l2 = Line(p1, p3)
        >>> l2.slope
        oo

        """
        d1, d2 = (self.p1 - self.p2).args
        if d1 == 0:
            return oo
        return simplify(d2/d1)

    @property
    def points(self):
        """The two points used to define this linear entity.

        Returns
        =======

        points : tuple of Points

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 11)
        >>> l1 = Line(p1, p2)
        >>> l1.points
        (Point(0, 0), Point(5, 11))

        """
        return self.p1, self.p2

    def projection(self, o):
        """Project a point, line, ray, or segment onto this linear entity.

        Parameters
        ==========

        other : diofant.geometry.point.Point or LinearEntity

        Returns
        =======

        projection : diofant.geometry.point.Point or LinearEntity
            The return type matches the type of the parameter ``other``.

        Raises
        ======

        diofant.geometry.exceptions.GeometryError
            When method is unable to perform projection.

        Notes
        =====

        A projection involves taking the two points that define
        the linear entity and projecting those points onto a
        Line and then reforming the linear entity using these
        projections.
        A point P is projected onto a line L by finding the point
        on L that is closest to P. This point is the intersection
        of L and the line perpendicular to L that passes through P.

        See Also
        ========

        diofant.geometry.point.Point, perpendicular_line

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(Rational(1, 2), 0)
        >>> l1 = Line(p1, p2)
        >>> l1.projection(p3)
        Point(1/4, 1/4)

        >>> p4, p5 = Point(10, 0), Point(12, 1)
        >>> s1 = Segment(p4, p5)
        >>> l1.projection(s1)
        Segment(Point(5, 5), Point(13/2, 13/2))

        """
        tline = Line(self.p1, self.p2)

        def _project(p):
            """Project a point onto the line representing self."""
            if p in tline:
                return p
            l1 = tline.perpendicular_line(p)
            return tline.intersection(l1)[0]

        projected = None
        if isinstance(o, Point):
            return _project(o)
        elif isinstance(o, LinearEntity):
            n_p1 = _project(o.p1)
            n_p2 = _project(o.p2)
            if n_p1 == n_p2:
                projected = n_p1
            else:
                projected = o.__class__(n_p1, n_p2)

        # Didn't know how to project so raise an error
        if projected is None:
            n1 = self.__class__.__name__
            n2 = o.__class__.__name__
            raise GeometryError(
                f'Do not know how to project {n2} onto {n1}')

        return self.intersection(projected)[0]

    def intersection(self, o):
        """The intersection with another geometrical entity.

        Parameters
        ==========

        o : diofant.geometry.point.Point or LinearEntity

        Returns
        =======

        intersection : list of geometrical entities

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(7, 7)
        >>> l1 = Line(p1, p2)
        >>> l1.intersection(p3)
        [Point(7, 7)]

        >>> p4, p5 = Point(5, 0), Point(0, 3)
        >>> l2 = Line(p4, p5)
        >>> l1.intersection(l2)
        [Point(15/8, 15/8)]

        >>> p6, p7 = Point(0, 5), Point(2, 6)
        >>> s1 = Segment(p6, p7)
        >>> l1.intersection(s1)
        []

        """
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, LinearEntity):
            a1, b1, c1 = self.coefficients
            a2, b2, c2 = o.coefficients
            t = simplify(a1*b2 - a2*b1)
            if t.equals(0) is not False:  # assume they are parallel
                if isinstance(self, Line):
                    if o.p1 in self:
                        return [o]
                    return []
                elif isinstance(o, Line):
                    if self.p1 in o:
                        return [self]
                    return []
                elif isinstance(self, Ray):
                    if isinstance(o, Ray):
                        # case 1, rays in the same direction
                        if self.xdirection == o.xdirection and \
                                self.ydirection == o.ydirection:
                            return [self] if (self.source in o) else [o]
                        # case 2, rays in the opposite directions
                        else:
                            if o.source in self:
                                if self.source == o.source:
                                    return [self.source]
                                return [Segment(o.source, self.source)]
                            return []
                    elif isinstance(o, Segment):
                        if o.p1 in self:
                            if o.p2 in self:
                                return [o]
                            return [Segment(o.p1, self.source)]
                        elif o.p2 in self:
                            return [Segment(o.p2, self.source)]
                        return []
                    else:
                        raise NotImplementedError
                elif isinstance(self, Segment):
                    if isinstance(o, Ray):
                        return o.intersection(self)
                    elif isinstance(o, Segment):
                        # A reminder that the points of Segments are ordered
                        # in such a way that the following works. See
                        # Segment.__new__ for details on the ordering.
                        if self.p1 not in o:
                            if self.p2 not in o:
                                # Neither of the endpoints are in o so either
                                # o is contained in this segment or it isn't
                                if o in self:
                                    return [self]
                                return []
                            else:
                                # p1 not in o but p2 is. Either there is a
                                # segment as an intersection, or they only
                                # intersect at an endpoint
                                if self.p2 == o.p1:
                                    return [o.p1]
                                return [Segment(o.p1, self.p2)]
                        elif self.p2 not in o:
                            # p2 not in o but p1 is. Either there is a
                            # segment as an intersection, or they only
                            # intersect at an endpoint
                            if self.p1 == o.p2:
                                return [o.p2]
                            return [Segment(o.p2, self.p1)]

                        # Both points of self in o so the whole segment
                        # is in o
                        return [self]
                    else:
                        raise NotImplementedError
                else:
                    raise NotImplementedError

            # Not parallel, so find the point of intersection
            px = simplify((b1*c2 - c1*b2) / t)
            py = simplify((a2*c1 - a1*c2) / t)
            inter = Point(px, py)
            # we do not use a simplistic 'inter in self and inter in o'
            # because that requires an equality test that is fragile;
            # instead we employ some diagnostics to see if the intersection
            # is valid

            def inseg(self):
                def _between(a, b, c):
                    return b >= c >= a or a >= c >= b
                if _between(self.p1.x, self.p2.x, inter.x) and \
                        _between(self.p1.y, self.p2.y, inter.y):
                    return True

            def inray(self):
                if self.p1 == inter:
                    return True
                sray = Ray(self.p1, inter)
                if sray.xdirection == self.xdirection and \
                        sray.ydirection == self.ydirection:
                    return True

            prec = (Line, Ray, Segment)
            expr = self
            if prec.index(expr.func) > prec.index(o.func):
                expr, o = o, expr
            rv = [inter]
            if isinstance(expr, Line):
                if isinstance(o, Line):
                    return rv
                elif isinstance(o, Ray) and inray(o):
                    return rv
                elif isinstance(o, Segment) and inseg(o):
                    return rv
            elif isinstance(expr, Ray) and inray(expr):
                if isinstance(o, Ray) and inray(o):
                    return rv
                elif isinstance(o, Segment) and inseg(o):
                    return rv
            elif isinstance(expr, Segment) and inseg(expr):
                if isinstance(o, Segment) and inseg(o):
                    return rv
            return []

        return o.intersection(self)

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the Line.

        Parameters
        ==========

        parameter : str, optional
            The name of the parameter which will be used for the parametric
            point. The default value is 't'. When this parameter is 0, the
            first point used to define the line will be returned, and when
            it is 1 the second point will be returned.

        Returns
        =======

        point : diofant.geometry.point.Point

        Raises
        ======

        ValueError
            When ``parameter`` already appears in the Line's definition.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(1, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.arbitrary_point()
        Point(4*t + 1, 3*t)

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError(f'Symbol {t.name} already appears in object '
                             'and cannot be used as a parameter.')
        # multiply on the right so the variable gets
        # combined witht he coordinates of the point
        return self.p1 + (self.p2 - self.p1)*t

    def random_point(self):
        """A random point on a LinearEntity.

        Returns
        =======

        point : diofant.geometry.point.Point

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> p3 = l1.random_point()
        >>> # random point - don't know its coords in advance
        >>> p3
        Point(...)
        >>> # point should belong to the line
        >>> p3 in l1
        True

        """
        from random import randint

        # The lower and upper
        lower, upper = -2**32 - 1, 2**32

        if self.slope is oo:
            if isinstance(self, Ray):
                if self.ydirection is oo:
                    lower = self.p1.y
                else:
                    upper = self.p1.y
            elif isinstance(self, Segment):
                lower = self.p1.y
                upper = self.p2.y

            x = self.p1.x
            y = randint(lower, upper)
        else:
            if isinstance(self, Ray):
                if self.xdirection is oo:
                    lower = self.p1.x
                else:
                    upper = self.p1.x
            elif isinstance(self, Segment):
                lower = self.p1.x
                upper = self.p2.x

            a, b, c = self.coefficients
            x = randint(lower, upper)
            y = (-c - a*x) / b
        return Point(x, y)

    def is_similar(self, other):
        """
        Return True if self and other are contained in the same line.

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 1), Point(3, 4), Point(2, 3)
        >>> l1 = Line(p1, p2)
        >>> l2 = Line(p1, p3)
        >>> l1.is_similar(l2)
        True

        """
        def _norm(a, b, c):
            if a != 0:
                return 1, b/a, c/a
            elif b != 0:
                return a/b, 1, c/b
            else:
                raise NotImplementedError
        return _norm(*self.coefficients) == _norm(*other.coefficients)

    def __contains__(self, other):
        """Return a definitive answer or else raise an error if it cannot
        be determined that other is on the boundaries of self.

        """
        result = self.contains(other)

        if result is not None:
            return result
        else:
            raise Undecidable(
                f"can't decide whether '{self}' contains '{other}'")

    def contains(self, other):
        """Subclasses should implement this method and should return
        True if other is on the boundaries of self;
        False if not on the boundaries of self;
        None if a determination cannot be made.

        """
        raise NotImplementedError()


class Line(LinearEntity):
    """An infinite line in space.

    A line is declared with two distinct points or a point and slope
    as defined using keyword `slope`.

    Notes
    =====

    At the moment only lines in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    Parameters
    ==========

    p1 : diofant.geometry.point.Point
    pt : diofant.geometry.point.Point
    slope : Expr

    See Also
    ========

    diofant.geometry.point.Point

    Examples
    ========

    >>> L = Line(Point(2, 3), Point(3, 5))
    >>> L
    Line(Point(2, 3), Point(3, 5))
    >>> L.points
    (Point(2, 3), Point(3, 5))
    >>> L.equation()
    -2*x + y + 1
    >>> L.coefficients
    (-2, 1, 1)

    Instantiate with keyword ``slope``:

    >>> Line(Point(0, 0), slope=0)
    Line(Point(0, 0), Point(1, 0))

    Instantiate with another linear object

    >>> s = Segment((0, 0), (0, 1))
    >>> Line(s).equation()
    x

    """

    def __new__(cls, p1, pt=None, slope=None, **kwargs):
        if isinstance(p1, LinearEntity):
            p1, pt = p1.args
        else:
            p1 = Point(p1)
        if pt is not None and slope is None:
            p2 = Point(pt)
        elif slope is not None and pt is None:
            slope = sympify(slope)
            if slope.is_finite is False:
                # when infinite slope, don't change x
                dx = 0
                dy = 1
            else:
                # go over 1 up slope
                dx = 1
                dy = slope
            # XXX avoiding simplification by adding to coords directly
            p2 = Point(p1.x + dx, p1.y + dy)
        else:
            raise ValueError('A 2nd Point or keyword "slope" must be used.')

        return LinearEntity.__new__(cls, p1, p2, **kwargs)

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of line. Gives
        values that will produce a line that is +/- 5 units long (where a
        unit is the distance between the two points that define the line).

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list (plot interval)
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.plot_interval()
        [t, -5, 5]

        """
        t = _symbol(parameter)
        return [t, -5, 5]

    def equation(self, x='x', y='y'):
        """The equation of the line: ax + by + c.

        Parameters
        ==========

        x : str, optional
            The name to use for the x-axis, default value is 'x'.
        y : str, optional
            The name to use for the y-axis, default value is 'y'.

        Returns
        =======

        equation : diofant expression

        See Also
        ========

        LinearEntity.coefficients

        Examples
        ========

        >>> p1, p2 = Point(1, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.equation()
        -3*x + 4*y + 3

        """
        x, y = _symbol(x), _symbol(y)
        p1, p2 = self.points
        if p1.x == p2.x:
            return x - p1.x
        elif p1.y == p2.y:
            return y - p1.y

        a, b, c = self.coefficients
        return a*x + b*y + c

    def contains(self, other):
        """
        Return True if o is on this Line, or False otherwise.

        Examples
        ========

        >>> p1, p2 = Point(0, 1), Point(3, 4)
        >>> l = Line(p1, p2)
        >>> l.contains(p1)
        True
        >>> l.contains((0, 1))
        True
        >>> l.contains((0, 0))
        False

        """
        if is_sequence(other):
            other = Point(other)
        if isinstance(other, Point):
            other = other.func(*[simplify(i) for i in other.args])
            x, y = Dummy(), Dummy()
            eq = self.equation(x, y)
            if not eq.has(y):
                return (solve(eq, x)[0][x] - other.x).equals(0)
            if not eq.has(x):
                return (solve(eq, y)[0][y] - other.y).equals(0)
            return (solve(eq.subs({x: other.x}), y)[0][y] - other.y).equals(0)
        elif not isinstance(other, LinearEntity):
            return False
        elif isinstance(other, Line):
            return self.equal(other)
        elif not self.is_similar(other):
            return False
        else:
            return other.p1 in self and other.p2 in self

    def distance(self, o):
        """
        Finds the shortest distance between a line and a point.

        Raises
        ======

        NotImplementedError
            if o is not a Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> s = Line(p1, p2)
        >>> s.distance(Point(-1, 1))
        sqrt(2)
        >>> s.distance((-1, 2))
        3*sqrt(2)/2

        """
        if not isinstance(o, Point):
            o = Point(o)
        a, b, c = self.coefficients
        if 0 in (a, b):
            return self.perpendicular_segment(o).length
        m = self.slope
        x = o.x
        y = m*x - c/b
        return abs(factor_terms(o.y - y))/sqrt(1 + m**2)

    def equal(self, other):
        """Returns True if self and other are the same mathematical entities."""
        if not isinstance(other, Line):
            return False
        return Point.is_collinear(self.p1, other.p1, self.p2, other.p2)


class Ray(LinearEntity):
    """
    A Ray is a semi-line in the space with a source point and a direction.

    Parameters
    ==========

    p1 : diofant.geometry.point.Point
        The source of the Ray
    p2 : diofant.geometry.point.Point or Expr
        This point determines the direction in which the Ray propagates.
        If given as an angle it is interpreted in radians with the positive
        direction being ccw.

    See Also
    ========

    diofant.geometry.point.Point, Line

    Notes
    =====

    At the moment only rays in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    Examples
    ========

    >>> r = Ray(Point(2, 3), Point(3, 5))
    >>> r = Ray(Point(2, 3), Point(3, 5))
    >>> r
    Ray(Point(2, 3), Point(3, 5))
    >>> r.points
    (Point(2, 3), Point(3, 5))
    >>> r.source
    Point(2, 3)
    >>> r.xdirection
    oo
    >>> r.ydirection
    oo
    >>> r.slope
    2
    >>> Ray(Point(0, 0), angle=pi/4).slope
    1

    """

    def __new__(cls, p1, pt=None, angle=None, **kwargs):
        p1 = Point(p1)
        if pt is not None and angle is None:
            try:
                p2 = Point(pt)
            except ValueError as exc:
                from ..utilities import filldedent
                raise ValueError(filldedent("""
                    The 2nd argument was not a valid Point; if
                    it was meant to be an angle it should be
                    given with keyword "angle".""")) from exc
            if p1 == p2:
                raise ValueError('A Ray requires two distinct points.')
        elif angle is not None and pt is None:
            # we need to know if the angle is an odd multiple of pi/2
            c = pi_coeff(sympify(angle))
            p2 = None
            if c is not None:
                if c.is_Rational:
                    if c.denominator == 2:
                        if c.numerator == 1:
                            p2 = p1 + Point(0, 1)
                        else:
                            assert c.numerator == 3
                            p2 = p1 + Point(0, -1)
                    elif c.denominator == 1:
                        if c.numerator == 0:
                            p2 = p1 + Point(1, 0)
                        else:
                            assert c.numerator == 1
                            p2 = p1 + Point(-1, 0)
                if p2 is None:
                    c *= pi
            else:
                c = angle % (2*pi)
            if not p2:
                m = 2*c/pi
                left = And(1 < m, m < 3)  # is it in quadrant 2 or 3?
                x = Piecewise((-1, left), (Piecewise((0, Eq(m % 1, 0)), (1, True)), True))
                y = Piecewise((-tan(c), left), (Piecewise((1, Eq(m, 1)), (-1, Eq(m, 3)), (tan(c), True)), True))
                p2 = p1 + Point(x, y)
        else:
            raise ValueError('A 2nd point or keyword "angle" must be used.')

        return LinearEntity.__new__(cls, p1, p2, **kwargs)

    @property
    def source(self):
        """The point from which the ray emanates.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(4, 1)
        >>> r1 = Ray(p1, p2)
        >>> r1.source
        Point(0, 0)

        """
        return self.p1

    @property
    def direction(self):
        """The direction in which the ray emanates.

        See Also
        ========

        diofant.geometry.point.Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(4, 1)
        >>> r1 = Ray(p1, p2)
        >>> r1.direction
        Point(4, 1)

        """
        return self.p2 - self.p1

    @property
    def xdirection(self):
        """The x direction of the ray.

        Positive infinity if the ray points in the positive x direction,
        negative infinity if the ray points in the negative x direction,
        or 0 if the ray is vertical.

        See Also
        ========

        ydirection

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(0, -1)
        >>> r1, r2 = Ray(p1, p2), Ray(p1, p3)
        >>> r1.xdirection
        oo
        >>> r2.xdirection
        0

        """
        if self.p1.x < self.p2.x:
            return oo
        elif self.p1.x == self.p2.x:
            return Integer(0)
        else:
            return -oo

    @property
    def ydirection(self):
        """The y direction of the ray.

        Positive infinity if the ray points in the positive y direction,
        negative infinity if the ray points in the negative y direction,
        or 0 if the ray is horizontal.

        See Also
        ========

        xdirection

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(-1, -1), Point(-1, 0)
        >>> r1, r2 = Ray(p1, p2), Ray(p1, p3)
        >>> r1.ydirection
        -oo
        >>> r2.ydirection
        0

        """
        if self.p1.y < self.p2.y:
            return oo
        elif self.p1.y == self.p2.y:
            return Integer(0)
        else:
            return -oo

    def distance(self, o):
        """
        Finds the shortest distance between the ray and a point.

        Raises
        ======

        NotImplementedError
            if o is not a Point

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> s = Ray(p1, p2)
        >>> s.distance(Point(-1, -1))
        sqrt(2)
        >>> s.distance((-1, 2))
        3*sqrt(2)/2

        """
        if not isinstance(o, Point):
            o = Point(o)
        s = self.perpendicular_segment(o)
        if isinstance(s, Point):
            if self.contains(s):
                return Integer(0)
        else:
            # since arg-order is arbitrary, find the non-o point
            non_o = s.p1 if s.p1 != o else s.p2
            if self.contains(non_o):
                return Line(self).distance(o)  # = s.length but simpler
        # the following applies when neither of the above apply
        return self.source.distance(o)

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Ray. Gives
        values that will produce a ray that is 10 units long (where a unit is
        the distance between the two points that define the ray).

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

        >>> r = Ray((0, 0), angle=pi/4)
        >>> r.plot_interval()
        [t, 0, 10]

        """
        t = _symbol(parameter)
        return [t, 0, 10]

    def contains(self, other):
        """
        Is other GeometryEntity contained in this Ray?

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(4, 4)
        >>> r = Ray(p1, p2)
        >>> r.contains(p1)
        True
        >>> r.contains((1, 1))
        True
        >>> r.contains((1, 3))
        False
        >>> s = Segment((1, 1), (2, 2))
        >>> r.contains(s)
        True
        >>> s = Segment((1, 2), (2, 5))
        >>> r.contains(s)
        False
        >>> r1 = Ray((2, 2), (3, 3))
        >>> r.contains(r1)
        True
        >>> r1 = Ray((2, 2), (3, 5))
        >>> r.contains(r1)
        False

        """
        if isinstance(other, Ray):
            return (Point.is_collinear(self.p1, self.p2, other.p1, other.p2) and
                    self.xdirection == other.xdirection and
                    self.ydirection == other.ydirection)
        elif isinstance(other, Segment):
            return other.p1 in self and other.p2 in self
        elif is_sequence(other):
            other = Point(other)
        if isinstance(other, Point):
            if Point.is_collinear(self.p1, self.p2, other):
                if self.xdirection is oo:
                    rv = other.x >= self.source.x
                elif self.xdirection == -oo:
                    rv = other.x <= self.source.x
                elif self.ydirection is oo:
                    rv = other.y >= self.source.y
                else:
                    rv = other.y <= self.source.y
                if rv in (true, false):
                    return bool(rv)
                raise Undecidable(f'Cannot determine if {other} is in {self}')
            # Points are not collinear, so the rays are not parallel
            # and hence it is impossible for self to contain o
            return False

        # No other known entity can be contained in a Ray
        return False


class Segment(LinearEntity):
    """A undirected line segment in space.

    Parameters
    ==========

    p1 : diofant.geometry.point.Point
    p2 : diofant.geometry.point.Point

    See Also
    ========

    diofant.geometry.point.Point, Line

    Notes
    =====

    At the moment only segments in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    Examples
    ========

    >>> Segment((1, 0), (1, 1))  # tuples are interpreted as pts
    Segment(Point(1, 0), Point(1, 1))
    >>> s = Segment(Point(4, 3), Point(1, 1))
    >>> s
    Segment(Point(1, 1), Point(4, 3))
    >>> s.points
    (Point(1, 1), Point(4, 3))
    >>> s.slope
    2/3
    >>> s.length
    sqrt(13)
    >>> s.midpoint
    Point(5/2, 2)

    """

    def __new__(cls, p1, p2, **kwargs):
        # Reorder the two points under the following ordering:
        #   if p1.x != p2.x then p1.x < p2.x
        #   if p1.x == p2.x then p1.y < p2.y
        p1 = Point(p1)
        p2 = Point(p2)
        if p1 == p2:
            return Point(p1)
        if (p1.x - p2.x).is_positive:
            p1, p2 = p2, p1
        elif (p1.x == p2.x) and (p1.y - p2.y).is_positive:
            p1, p2 = p2, p1
        return LinearEntity.__new__(cls, p1, p2, **kwargs)

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Segment gives
        values that will produce the full segment in a plot.

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

        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> s1 = Segment(p1, p2)
        >>> s1.plot_interval()
        [t, 0, 1]

        """
        t = _symbol(parameter)
        return [t, 0, 1]

    def perpendicular_bisector(self, p=None):
        """The perpendicular bisector of this segment.

        If no point is specified or the point specified is not on the
        bisector then the bisector is returned as a Line. Otherwise a
        Segment is returned that joins the point specified and the
        intersection of the bisector and the segment.

        Parameters
        ==========

        p : diofant.geometry.point.Point

        Returns
        =======

        bisector : Line or Segment

        See Also
        ========

        LinearEntity.perpendicular_segment

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(6, 6), Point(5, 1)
        >>> s1 = Segment(p1, p2)
        >>> s1.perpendicular_bisector()
        Line(Point(3, 3), Point(9, -3))

        >>> s1.perpendicular_bisector(p3)
        Segment(Point(3, 3), Point(5, 1))

        """
        l = LinearEntity.perpendicular_line(self, self.midpoint)
        if p is None or Point(p) not in l:
            return l
        else:
            return Segment(self.midpoint, p)

    @property
    def length(self):
        """The length of the line segment.

        See Also
        ========

        diofant.geometry.point.Point.distance

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(4, 3)
        >>> s1 = Segment(p1, p2)
        >>> s1.length
        5

        """
        return Point.distance(self.p1, self.p2)

    @property
    def midpoint(self):
        """The midpoint of the line segment.

        See Also
        ========

        diofant.geometry.point.Point.midpoint

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(4, 3)
        >>> s1 = Segment(p1, p2)
        >>> s1.midpoint
        Point(2, 3/2)

        """
        return Point.midpoint(self.p1, self.p2)

    def distance(self, o):
        """
        Finds the shortest distance between a line segment and a point.

        Raises
        ======

        NotImplementedError
            if o is not a Point

        Examples
        ========

        >>> p1, p2 = Point(0, 1), Point(3, 4)
        >>> s = Segment(p1, p2)
        >>> s.distance(Point(10, 15))
        sqrt(170)
        >>> s.distance((0, 12))
        sqrt(73)

        """
        if is_sequence(o):
            o = Point(o)
        if isinstance(o, Point):
            seg_vector = self.p2 - self.p1
            pt_vector = o - self.p1
            t = seg_vector.dot(pt_vector)/self.length**2
            if t >= 1:
                distance = Point.distance(self.p2, o)
            elif t <= 0:
                distance = Point.distance(self.p1, o)
            else:
                distance = Point.distance(
                    self.p1 + Point(t*seg_vector.x, t*seg_vector.y), o)
            return distance
        else:
            raise NotImplementedError

    def contains(self, other):
        """
        Is the other GeometryEntity contained within this Segment?

        Examples
        ========

        >>> p1, p2 = Point(0, 1), Point(3, 4)
        >>> s = Segment(p1, p2)
        >>> s2 = Segment(p2, p1)
        >>> s.contains(s2)
        True

        """
        if isinstance(other, Segment):
            return other.p1 in self and other.p2 in self
        elif isinstance(other, Point):
            if Point.is_collinear(self.p1, self.p2, other):
                t = Dummy('t')
                x, y = self.arbitrary_point(t).args
                if self.p1.x != self.p2.x:
                    ti = solve(x - other.x, t)[0][t]
                else:
                    ti = solve(y - other.y, t)[0][t]
                if ti.is_number:
                    return 0 <= ti <= 1
                return

        return False
