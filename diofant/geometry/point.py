"""Geometrical Points.

Contains
========
Point
"""

from ..core import Add, Float, Integer, Tuple
from ..core.compatibility import iterable
from ..core.evaluate import global_evaluate
from ..core.sympify import sympify
from ..functions import im, sqrt
from ..matrices import Matrix
from ..simplify import nsimplify, simplify
from ..utilities import ordered
from .entity import GeometryEntity


class Point(GeometryEntity):
    """A point in a n-dimensional Euclidean space.

    Parameters
    ==========

    coords : sequence of n-coordinate values.

    Raises
    ======

    TypeError
        When trying to add or subtract points with different dimensions.
        When `intersection` is called with object other than a Point.

    See Also
    ========

    diofant.geometry.line.Segment : Connects two Points

    Examples
    ========

    >>> Point([1, 2])
    Point(1, 2)
    >>> Point(0, x)
    Point(0, x)

    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point(0.5, 0.25)
    Point(1/2, 1/4)
    >>> print(Point(0.5, 0.25, evaluate=False))
    Point(0.5, 0.25)

    """

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        if iterable(args[0]):
            args = args[0]

        # unpack the arguments into a friendly Tuple
        # if we were already a Point, we're doing an excess
        # iteration, but we'll worry about efficiency later
        coords = Tuple(*args)
        if any(a.is_number and im(a) for a in coords):
            raise ValueError('Imaginary coordinates not permitted.')

        # Turn any Floats into rationals and simplify
        # any expressions before we instantiate
        if evaluate:
            coords = coords.xreplace({f: simplify(nsimplify(f, rational=True))
                                      for f in coords.atoms(Float)})
        if len(coords) != 2:
            raise ValueError('Only two dimensional points currently supported')

        return GeometryEntity.__new__(cls, *coords)

    def __contains__(self, item):
        return item == self

    def is_concyclic(*points):
        """Is a sequence of points concyclic?

        Test whether or not a sequence of points are concyclic (i.e., they lie
        on a circle).

        Parameters
        ==========

        points : sequence of Points

        Returns
        =======

        is_concyclic : boolean
            True if points are concyclic, False otherwise.

        See Also
        ========

        diofant.geometry.ellipse.Circle

        Notes
        =====

        No points are not considered to be concyclic. One or two points
        are definitely concyclic and three points are conyclic iff they
        are not collinear.

        For more than three points, create a circle from the first three
        points. If the circle cannot be created (i.e., they are collinear)
        then all of the points cannot be concyclic. If the circle is created
        successfully then simply check the remaining points for containment
        in the circle.

        Examples
        ========

        >>> p1, p2 = Point(-1, 0), Point(1, 0)
        >>> p3, p4 = Point(0, 1), Point(-1, 2)
        >>> Point.is_concyclic(p1, p2, p3)
        True
        >>> Point.is_concyclic(p1, p2, p3, p4)
        False

        """
        if not all(isinstance(p, Point) for p in points):
            raise TypeError('Must pass only Point objects')

        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True

        ppoints = list(ordered(Point(p) for p in points))

        from .ellipse import Circle
        try:
            c = Circle(ppoints[0], ppoints[1], ppoints[2])
        except ValueError:
            return False
        for point in ppoints[3:]:
            if point not in c:
                return False
        return True

    def is_collinear(*args):
        """Is a sequence of points collinear?

        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Parameters
        ==========

        points : sequence of Point

        Returns
        =======

        is_collinear : boolean

        Notes
        =====

        Slope is preserved everywhere on a line, so the slope between
        any two points on the line should be the same. Take the first
        two points, p1 and p2, and create a translated point v1
        with p1 as the origin. Now for every other point we create
        a translated point, vi with p1 also as the origin. Note that
        these translations preserve slope since everything is
        consistently translated to a new origin of p1. Since slope
        is preserved then we have the following equality:

              * v1_slope = vi_slope
              * v1.y/v1.x = vi.y/vi.x (due to translation)
              * v1.y*vi.x = vi.y*v1.x
              * v1.y*vi.x - vi.y*v1.x = 0           (*)

        Hence, if we have a vi such that the equality in (*) is False
        then the points are not collinear. We do this test for every
        point in the list, and if all pass then they are collinear.

        See Also
        ========

        diofant.geometry.line.Line

        Examples
        ========

        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(x, x), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False

        """
        # Coincident points are irrelevant; use only unique points.
        uniq_args = list(set(args))
        if not all(isinstance(p, Point) for p in uniq_args):
            raise TypeError('Must pass only Point objects')

        if len(uniq_args) == 0:
            return False
        if len(uniq_args) <= 2:
            return True

        # translate our points
        points = [p - uniq_args[0] for p in uniq_args[1:]]
        for p in points[1:]:
            if not Point.is_scalar_multiple(points[0], p):
                return False
        return True

    def is_scalar_multiple(self, other):
        """Returns whether `self` and `other` are scalar multiples
        of each other.

        """
        # if the vectors self and other are linearly dependent, then they must
        # be scalar multiples of each other
        m = Matrix([self.args, other.args])
        return m.rank() < 2

    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> p = Point(0, 1)
        >>> p.length
        0

        """
        return Integer(0)

    @property
    def origin(self):
        """A point of all zeros of the same ambient dimension
        as the current point

        """
        return Point([0]*len(self))

    @property
    def is_zero(self):
        """True if every coordinate is zero, otherwise False."""
        return all(x == 0 for x in self.args)

    @property
    def ambient_dimension(self):
        """The dimension of the ambient space the point is in.
        I.e., if the point is in R^n, the ambient dimension
        will be n

        """
        return len(self)

    def distance(self, p):
        """The Euclidean distance from self to point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        distance : number or symbolic expression.

        See Also
        ========

        diofant.geometry.line.Segment.length

        Examples
        ========

        >>> p1, p2 = Point(1, 1), Point(4, 5)
        >>> p1.distance(p2)
        5

        >>> p3 = Point(x, y)
        >>> p3.distance(Point(0, 0))
        sqrt(x**2 + y**2)

        """
        p = Point(p)
        return sqrt(sum((a - b)**2
                        for a, b in zip(self.args,
                                        p.args if isinstance(p, Point) else p)))

    def midpoint(self, p):
        """The midpoint between self and point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        midpoint : Point

        See Also
        ========

        diofant.geometry.line.Segment.midpoint

        Examples
        ========

        >>> p1, p2 = Point(1, 1), Point(13, 5)
        >>> p1.midpoint(p2)
        Point(7, 3)

        """
        return Point([simplify((a + b)/2) for a, b in zip(self.args, p.args)])

    def evalf(self, dps=15, **options):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers to
        the decimal precision dps.

        Returns
        =======

        point : Point

        Examples
        ========

        >>> p1 = Point(Rational(1, 2), Rational(3, 2))
        >>> p1
        Point(1/2, 3/2)
        >>> print(p1.evalf())
        Point(0.5, 1.5)

        """
        coords = [x.evalf(dps, **options) for x in self.args]
        return Point(*coords, evaluate=False)

    n = evalf

    def intersection(self, o):
        """The intersection between this point and another point.

        Parameters
        ==========

        other : Point

        Returns
        =======

        intersection : list of Points

        Notes
        =====

        The return value will either be an empty list if there is no
        intersection, otherwise it will contain this point.

        Examples
        ========

        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(0, 0)
        >>> p1.intersection(p2)
        []
        >>> p1.intersection(p3)
        [Point(0, 0)]

        """
        if isinstance(o, Point):
            if self == o:
                return [self]
            return []

        return o.intersection(self)

    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = Point(p2)
        return Add(*[a*b for a, b in zip(self, p2)])

    def __len__(self):
        return len(self.args)

    def __iter__(self):
        return self.args.__iter__()

    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        diofant.geometry.entity.GeometryEntity.translate

        """
        if iterable(other) and len(other) == len(self):
            return Point([simplify(a + b) for a, b in zip(self, other)])
        else:
            raise ValueError(
                'Points must have the same number of dimensions')

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates.

        """
        return self + (-other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point([simplify(x*factor) for x in self.args])

    def __truediv__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point([simplify(x/divisor) for x in self.args])

    def __neg__(self):
        """Negate the point."""
        return Point([-x for x in self.args])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0]*len(self))
        return Point.distance(origin, self)

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> p = Point(0, 1)
        >>> p.x
        0

        """
        return self.args[0]

    @property
    def y(self):
        """
        Returns the Y coordinate of the Point.

        Examples
        ========

        >>> p = Point(0, 1)
        >>> p.y
        1

        """
        return self.args[1]

    @property
    def bounds(self):
        """Return a tuple (xmin, ymin, xmax, ymax) representing the bounding
        rectangle for the geometric figure.

        """
        return self.x, self.y, self.x, self.y

    def rotate(self, angle, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> t = Point(1, 0)
        >>> t.rotate(pi/2)
        Point(0, 1)
        >>> t.rotate(pi/2, (2, 0))
        Point(2, -1)

        """
        from ..functions import cos, sin

        c = cos(angle)
        s = sin(angle)

        rv = self
        if pt is not None:
            pt = Point(pt)
            rv -= pt
        x, y = rv.args
        rv = Point(c*x - s*y, s*x + c*y)
        if pt is not None:
            rv += pt
        return rv

    def scale(self, x=1, y=1, pt=None):
        """Scale the coordinates of the Point by multiplying by
        ``x`` and ``y`` after subtracting ``pt`` -- default is (0, 0) --
        and then adding ``pt`` back again (i.e. ``pt`` is the point of
        reference for the scaling).

        See Also
        ========

        rotate, translate

        Examples
        ========

        >>> t = Point(1, 1)
        >>> t.scale(2)
        Point(2, 1)
        >>> t.scale(2, 2)
        Point(2, 2)

        """
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        return Point(self.x*x, self.y*y)

    def translate(self, x=0, y=0):
        """Shift the Point by adding x and y to the coordinates of the Point.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> t = Point(0, 1)
        >>> t.translate(2)
        Point(2, 1)
        >>> t.translate(2, 2)
        Point(2, 3)
        >>> t + Point(2, 2)
        Point(2, 3)

        """
        return Point(self.x + x, self.y + y)

    def transform(self, matrix):
        """Return the point after applying the transformation described
        by the 3x3 Matrix, ``matrix``.

        See Also
        ========

        diofant.geometry.entity.GeometryEntity.rotate
        diofant.geometry.entity.GeometryEntity.scale
        diofant.geometry.entity.GeometryEntity.translate

        """
        try:
            col, row = matrix.shape
            valid_matrix = matrix.is_square and col == 3
        except AttributeError:
            # We hit this block if matrix argument is not actually a Matrix.
            valid_matrix = False
        if not valid_matrix:
            raise ValueError('The argument to the transform function must be '
                             + 'a 3x3 matrix')
        x, y = self.args
        return Point(*(Matrix(1, 3, [x, y, 1])*matrix).tolist()[0][:2])
