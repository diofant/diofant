"""The definition of the base geometrical entity with attributes common to
all derived geometrical entities.

Contains
========

GeometryEntity
GeometricSet

Notes
=====

A GeometryEntity is any object that has special geometric properties.
A GeometrySet is a superclass of any GeometryEntity that can also
be viewed as a diofant.sets.Set.  In particular, points are the only
GeometryEntity not considered a Set.

Rn is a GeometrySet representing n-dimensional Euclidean space. R2 and
R3 are currently the only ambient spaces implemented.

"""

from ..core import Basic, Dummy, Tuple, oo
from ..core.compatibility import is_sequence
from ..core.sympify import sympify
from ..functions import atan, cos, sin
from ..matrices import eye
from ..sets import Set


class GeometryEntity(Basic):
    """The base class for all geometrical entities.

    This class doesn't represent any particular geometric entity, it only
    provides the implementation of some methods common to all subclasses.

    """

    def __new__(cls, *args, **kwargs):
        from .point import Point
        args = [Tuple(*a) if is_sequence(a)
                and not isinstance(a, Point) else sympify(a) for a in args]
        return Basic.__new__(cls, *args)

    def intersection(self, o):
        """
        Returns a list of all of the intersections of self with o.

        Notes
        =====

        An entity is not required to implement this method.

        See Also
        ========

        diofant.geometry.util.intersection

        """
        raise NotImplementedError()

    def rotate(self, angle, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        The default pt is the origin, Point(0, 0)

        See Also
        ========

        scale, translate

        Examples
        ========

        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t  # vertex on x axis
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.rotate(pi/2)  # vertex on y axis now
        Triangle(Point(0, 1), Point(-sqrt(3)/2, -1/2), Point(sqrt(3)/2, -1/2))

        """
        newargs = []
        for a in self.args:
            if isinstance(a, GeometryEntity):
                newargs.append(a.rotate(angle, pt))
            else:
                newargs.append(a)
        return type(self)(*newargs)

    def scale(self, x=1, y=1, pt=None):
        """Scale the object by multiplying the x,y-coordinates by x and y.

        If pt is given, the scaling is done relative to that point; the
        object is shifted by -pt, scaled, and shifted by pt.

        See Also
        ========

        rotate, translate

        Examples
        ========

        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.scale(2)
        Triangle(Point(2, 0), Point(-1, sqrt(3)/2), Point(-1, -sqrt(3)/2))
        >>> t.scale(2, 2)
        Triangle(Point(2, 0), Point(-1, sqrt(3)), Point(-1, -sqrt(3)))

        """
        from .point import Point
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        return type(self)(*[a.scale(x, y) for a in self.args])  # if this fails, override this class

    def translate(self, x=0, y=0):
        """Shift the object by adding to the x,y-coordinates the values x and y.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.translate(2)
        Triangle(Point(3, 0), Point(3/2, sqrt(3)/2), Point(3/2, -sqrt(3)/2))
        >>> t.translate(2, 2)
        Triangle(Point(3, 2), Point(3/2, sqrt(3)/2 + 2),
            Point(3/2, -sqrt(3)/2 + 2))

        """
        newargs = []
        for a in self.args:
            if isinstance(a, GeometryEntity):
                newargs.append(a.translate(x, y))
            else:
                newargs.append(a)
        return self.func(*newargs)

    def reflect(self, line):
        from .point import Point
        g = self
        l = line
        o = Point(0, 0)
        if l.slope == 0:
            y = l.args[0].y
            if not y:  # x-axis
                return g.scale(y=-1)
            reps = [(p, p.translate(y=2*(y - p.y))) for p in g.atoms(Point)]
        elif l.slope == oo:
            x = l.args[0].x
            if not x:  # y-axis
                return g.scale(x=-1)
            reps = [(p, p.translate(x=2*(x - p.x))) for p in g.atoms(Point)]
        elif hasattr(g, 'reflect'):
            a = atan(l.slope)
            c = l.coefficients
            d = -c[-1]/c[1]  # y-intercept
            # apply the transform to a single point
            x, y = Dummy(), Dummy()
            xf = Point(x, y)
            xf = xf.translate(y=-d).rotate(-a, o)
            xf = xf.scale(y=-1).rotate(a, o).translate(y=d)
            # replace every point using that transform
            reps = [(p, xf.xreplace({x: p.x, y: p.y})) for p in g.atoms(Point)]
        else:
            raise NotImplementedError(f'reflect undefined or non-Point args in {g}')
        return g.xreplace(dict(reps))

    def encloses(self, o):
        """
        Return True if o is inside (not on or outside) the boundaries of self.

        The object will be decomposed into Points and individual Entities need
        only define an encloses_point method for their class.

        See Also
        ========

        diofant.geometry.ellipse.Ellipse.encloses_point
        diofant.geometry.polygon.Polygon.encloses_point

        Examples
        ========

        >>> t1 = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t2 = Polygon(*RegularPolygon(Point(0, 0), 2, 3).vertices)
        >>> t2.encloses(t1)
        True
        >>> t1.encloses(t2)
        False

        """
        from .ellipse import Ellipse
        from .line import Line, Ray, Segment
        from .point import Point
        from .polygon import Polygon, RegularPolygon

        if isinstance(o, Point):
            return self.encloses_point(o)
        elif isinstance(o, Segment):
            return all(self.encloses_point(x) for x in o.points)
        elif isinstance(o, Ray) or isinstance(o, Line):
            return False
        elif isinstance(o, Ellipse):
            return (self.encloses_point(o.center) and
                    not self.intersection(o) and
                    self.encloses_point(Point(o.center.x + o.hradius,
                                              o.center.y)))
        elif isinstance(o, Polygon):
            if isinstance(o, RegularPolygon):
                if not self.encloses_point(o.center):
                    return False
            return all(self.encloses_point(v) for v in o.vertices)
        raise NotImplementedError()

    @property
    def ambient_dimension(self):
        """What is the dimension of the space that the object is contained in?"""
        raise NotImplementedError()

    def is_similar(self, other):
        """Is this geometrical entity similar to another geometrical entity?

        Two entities are similar if a uniform scaling (enlarging or
        shrinking) of one of the entities will allow one to obtain the other.

        Notes
        =====

        This method is not intended to be used directly but rather
        through the `are_similar` function found in util.py.
        An entity is not required to implement this method.
        If two different types of entities can be similar, it is only
        required that one of them be able to determine this.

        See Also
        ========

        scale

        """
        raise NotImplementedError()

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __str__(self):
        """String representation of a GeometryEntity."""
        from ..printing import sstr
        return type(self).__name__ + sstr(self.args)

    def __repr__(self):
        """String representation of a GeometryEntity that can be evaluated
        by diofant.

        """
        return type(self).__name__ + repr(self.args)

    def _eval_subs(self, old, new):
        from .point import Point
        if is_sequence(old) or is_sequence(new):
            old = Point(old)
            new = Point(new)
            return self._subs(old, new)


class GeometrySet(GeometryEntity, Set):
    """Parent class of all GeometryEntity that are also Sets
    (compatible with diofant.sets).

    """

    def _contains(self, other):
        """diofant.sets uses the _contains method, so include it for compatibility."""
        return self.__contains__(other)

    def _union(self, o):
        """Returns the union of self and o
        for use with diofant.sets.Set, if possible.

        """
        from ..sets import FiniteSet, Union

        # if its a FiniteSet, merge any points
        # we contain and return a union with the rest
        if o.is_FiniteSet:
            other_points = [p for p in o if not self._contains(p)]
            if len(other_points) == len(o):
                return
            return Union(self, FiniteSet(*other_points))
        if self._contains(o):
            return self

    def _intersection(self, o):
        """Returns a diofant.sets.Set of intersection objects,
        if possible.

        """
        from ..sets import FiniteSet, Union
        from .point import Point

        inter = self.intersection(o)

        # put the points in a FiniteSet
        points = FiniteSet(*[p for p in inter if isinstance(p, Point)])
        non_points = [p for p in inter if not isinstance(p, Point)]

        return Union(*(non_points + [points]))


def translate(x, y):
    """Return the matrix to translate a 2-D point by x and y."""
    rv = eye(3)
    rv[2, 0] = x
    rv[2, 1] = y
    return rv


def scale(x, y, pt=None):
    """Return the matrix to multiply a 2-D point's coordinates by x and y.

    If pt is given, the scaling is done relative to that point.

    """
    rv = eye(3)
    rv[0, 0] = x
    rv[1, 1] = y
    if pt:
        from .point import Point
        pt = Point(pt)
        tr1 = translate(*(-pt).args)
        tr2 = translate(*pt.args)
        return tr1*rv*tr2
    return rv


def rotate(th):
    """Return the matrix to rotate a 2-D point about the origin by ``angle``.

    The angle is measured in radians. To Point a point about a point other
    then the origin, translate the Point, do the rotation, and
    translate it back:

    >>> rot_about_11 = translate(-1, -1)*rotate(pi/2)*translate(1, 1)
    >>> Point(1, 1).transform(rot_about_11)
    Point(1, 1)
    >>> Point(0, 0).transform(rot_about_11)
    Point(2, 0)

    """
    s = sin(th)
    rv = eye(3)*cos(th)
    rv[0, 1] = s
    rv[1, 0] = -s
    rv[2, 2] = 1
    return rv
