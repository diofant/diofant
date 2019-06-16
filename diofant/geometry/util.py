"""Utility functions for geometrical entities.

Contains
========

idiff
intersection
convex_hull
are_coplanar
are_similar
"""

from ..core import Dummy, Function, Symbol
from ..core.compatibility import is_sequence
from ..solvers import solve


def idiff(eq, y, x, n=1):
    """Return ``dy/dx`` assuming that ``eq == 0``.

    Parameters
    ==========

    y : the dependent variable or a list of dependent variables (with y first)
    x : the variable that the derivative is being taken with respect to
    n : the order of the derivative (default is 1)

    Examples
    ========

    >>> from diofant.abc import a

    >>> circ = x**2 + y**2 - 4
    >>> idiff(circ, y, x)
    -x/y
    >>> idiff(circ, y, x, 2).simplify()
    -(x**2 + y**2)/y**3

    Here, ``a`` is assumed to be independent of ``x``:

    >>> idiff(x + a + y, y, x)
    -1

    Now the x-dependence of ``a`` is made explicit by listing ``a`` after
    ``y`` in a list.

    >>> idiff(x + a + y, [y, a], x)
    -Derivative(a, x) - 1

    See Also
    ========

    diofant.core.function.Derivative: represents unevaluated derivatives
    diofant.core.function.diff: explicitly differentiates wrt symbols

    """
    if is_sequence(y):
        dep = set(y)
        y = y[0]
    elif isinstance(y, (Dummy, Symbol)):
        dep = {y}
    else:
        raise ValueError("expecting x-dependent symbol(s) but got: %s" % y)

    f = {s: Function(s.name)(x) for s in eq.free_symbols if s != x and s in dep}
    dydx = Function(y.name)(x).diff(x)
    eq = eq.subs(f)
    derivs = {}
    for i in range(n):
        yp = solve(eq.diff(x), dydx)[0][dydx].subs(derivs)
        if i == n - 1:
            return yp.subs([(v, k) for k, v in f.items()])
        derivs[dydx] = yp
        eq = dydx - yp
        dydx = dydx.diff(x)


def _symbol(s, matching_symbol=None):
    """Return s if s is a Symbol, else return either a new Symbol (extended_real=True)
    with the same name s or the matching_symbol if s is a string and it matches
    the name of the matching_symbol.

    >>> _symbol('y')
    y
    >>> _.is_extended_real
    True
    >>> _symbol(x)
    x
    >>> _.is_extended_real is None
    True
    >>> arb = Symbol('foo')
    >>> _symbol('arb', arb) # arb's name is foo so foo will not be returned
    arb
    >>> _symbol('foo', arb) # now it will
    foo

    NB: the symbol here may not be the same as a symbol with the same
    name defined elsewhere as a result of different assumptions.

    See Also
    ========

    diofant.core.symbol.Symbol

    """
    if isinstance(s, str):
        if matching_symbol and matching_symbol.name == s:
            return matching_symbol
        return Symbol(s, extended_real=True)
    elif isinstance(s, (Dummy, Symbol)):
        return s
    else:
        raise ValueError('symbol must be string for symbol name or Symbol')


def _uniquely_named_symbol(xname, *exprs):
    """Return a symbol which, when printed, will have a name unique
    from any other already in the expressions given. The name is made
    unique by prepending underscores.
    """
    prefix = '%s'
    x = prefix % xname
    syms = set().union(*[e.free_symbols for e in exprs])
    while any(x == str(s) for s in syms):
        prefix = '_' + prefix
        x = prefix % xname
    return _symbol(x)


def intersection(*entities):
    """The intersection of a collection of GeometryEntity instances.

    Parameters
    ==========

    entities : sequence of GeometryEntity

    Returns
    =======

    intersection : list of GeometryEntity

    Raises
    ======

    NotImplementedError
        When unable to calculate intersection.

    Notes
    =====

    The intersection of any geometrical entity with itself should return
    a list with one item: the entity in question.
    An intersection requires two or more entities. If only a single
    entity is given then the function will return an empty list.
    It is possible for `intersection` to miss intersections that one
    knows exists because the required quantities were not fully
    simplified internally.
    Reals should be converted to Rationals, e.g. Rational(str(real_num))
    or else failures due to floating point issues may result.

    See Also
    ========

    diofant.geometry.entity.GeometryEntity.intersection

    Examples
    ========

    >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(-1, 5)
    >>> l1, l2 = Line(p1, p2), Line(p3, p2)
    >>> c = Circle(p2, 1)
    >>> intersection(l1, p2)
    [Point2D(1, 1)]
    >>> intersection(l1, l2)
    [Point2D(1, 1)]
    >>> intersection(c, p2)
    []
    >>> intersection(c, Point(1, 0))
    [Point2D(1, 0)]
    >>> intersection(c, l2)
    [Point2D(-sqrt(5)/5 + 1, 2*sqrt(5)/5 + 1),
     Point2D(sqrt(5)/5 + 1, -2*sqrt(5)/5 + 1)]

    """
    from .entity import GeometryEntity
    from .point import Point

    if len(entities) <= 1:
        return []

    # entities may be an immutable tuple
    entities = list(entities)
    for i, e in enumerate(entities):
        if not isinstance(e, GeometryEntity):
            try:
                entities[i] = Point(e)
            except NotImplementedError:
                raise ValueError('%s is not a GeometryEntity and cannot be made into Point' % str(e))

    res = entities[0].intersection(entities[1])
    for entity in entities[2:]:
        newres = []
        for x in res:
            newres.extend(x.intersection(entity))
        res = newres
    return res


def convex_hull(*args):
    """The convex hull surrounding the Points contained in the list of entities.

    Parameters
    ==========

    args : a collection of Points, Segments and/or Polygons

    Returns
    =======

    convex_hull : Polygon

    Notes
    =====

    This can only be performed on a set of non-symbolic points.

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Graham_scan

    [2] Andrew's Monotone Chain Algorithm
    (A.M. Andrew,
    "Another Efficient Algorithm for Convex Hulls in Two Dimensions", 1979)
    http://geomalgorithms.com/a10-_hull-1.html

    See Also
    ========

    diofant.geometry.point.Point, diofant.geometry.polygon.Polygon

    Examples
    ========

    >>> points = [(1, 1), (1, 2), (3, 1), (-5, 2), (15, 4)]
    >>> convex_hull(*points)
    Polygon(Point2D(-5, 2), Point2D(1, 1), Point2D(3, 1), Point2D(15, 4))

    """
    from .entity import GeometryEntity
    from .point import Point
    from .line import Segment
    from .polygon import Polygon

    p = set()
    for e in args:
        if not isinstance(e, GeometryEntity):
            try:
                e = Point(e)
            except NotImplementedError:
                raise ValueError('%s is not a GeometryEntity and cannot be made into Point' % str(e))
        if isinstance(e, Point):
            p.add(e)
        elif isinstance(e, Segment):
            p.update(e.points)
        elif isinstance(e, Polygon):
            p.update(e.vertices)
        else:
            raise NotImplementedError(
                'Convex hull for %s not implemented.' % type(e))

    # make sure all our points are of the same dimension
    if any(len(x) != 2 for x in p):
        raise ValueError('Can only compute the convex hull in two dimensions')

    p = list(p)
    if len(p) == 1:
        return p[0]
    elif len(p) == 2:
        return Segment(p[0], p[1])

    def _orientation(p, q, r):
        """Return positive if p-q-r are clockwise, neg if ccw, zero if
        collinear.
        """
        return (q.y - p.y)*(r.x - p.x) - (q.x - p.x)*(r.y - p.y)

    # scan to find upper and lower convex hulls of a set of 2d points.
    U = []
    L = []
    p.sort(key=lambda x: x.args)
    for p_i in p:
        while len(U) > 1 and _orientation(U[-2], U[-1], p_i) <= 0:
            U.pop()
        while len(L) > 1 and _orientation(L[-2], L[-1], p_i) >= 0:
            L.pop()
        U.append(p_i)
        L.append(p_i)
    U.reverse()
    convexHull = tuple(L + U[1:-1])

    if len(convexHull) == 2:
        return Segment(convexHull[0], convexHull[1])
    return Polygon(*convexHull)


def are_coplanar(*e):
    """ Returns True if the given entities are coplanar otherwise False

    Parameters
    ==========

    e: entities to be checked for being coplanar

    Returns
    =======

    Boolean

    Examples
    ========

    >>> a = Line3D(Point3D(5, 0, 0), Point3D(1, -1, 1))
    >>> b = Line3D(Point3D(0, -2, 0), Point3D(3, 1, 1))
    >>> c = Line3D(Point3D(0, -1, 0), Point3D(5, -1, 9))
    >>> are_coplanar(a, b, c)
    False

    """
    from .entity import GeometryEntity
    from .line3d import LinearEntity3D
    from .point import Point3D, Point
    from .plane import Plane
    # XXX update tests for coverage

    e = set(e)
    # first work with a Plane if present
    for i in list(e):
        if isinstance(i, Plane):
            e.remove(i)
            return all(p.is_coplanar(i) for p in e)

    if all(isinstance(i, Point3D) for i in e):
        if len(e) < 3:
            return False

        # remove pts that are collinear with 2 pts
        a, b = e.pop(), e.pop()
        for i in list(e):
            if Point3D.are_collinear(a, b, i):
                e.remove(i)

        if not e:
            return False
        else:
            # define a plane
            p = Plane(a, b, e.pop())
            for i in e:
                if i not in p:
                    return False
            return True
    else:
        pt3d = []
        for i in e:
            if isinstance(i, Point3D):
                pt3d.append(i)
            elif isinstance(i, LinearEntity3D):
                pt3d.extend(i.args)
            elif isinstance(i, GeometryEntity):  # XXX we should have a GeometryEntity3D class so we can tell the difference between 2D and 3D -- here we just want to deal with 2D objects; if new 3D objects are encountered that we didn't hanlde above, an error should be raised
                # all 2D objects have some Point that defines them; so convert those points to 3D pts by making z=0
                for p in i.args:
                    if isinstance(p, Point):
                        pt3d.append(Point3D(*(p.args + (0,))))
        return are_coplanar(*pt3d)


def are_similar(e1, e2):
    """Are two geometrical entities similar.

    Can one geometrical entity be uniformly scaled to the other?

    Parameters
    ==========

    e1 : GeometryEntity
    e2 : GeometryEntity

    Returns
    =======

    are_similar : boolean

    Raises
    ======

    diofant.geometry.exceptions.GeometryError
        When `e1` and `e2` cannot be compared.

    Notes
    =====

    If the two objects are equal then they are similar.

    See Also
    ========

    diofant.geometry.entity.GeometryEntity.is_similar

    Examples
    ========

    >>> c1, c2 = Circle(Point(0, 0), 4), Circle(Point(1, 4), 3)
    >>> t1 = Triangle(Point(0, 0), Point(1, 0), Point(0, 1))
    >>> t2 = Triangle(Point(0, 0), Point(2, 0), Point(0, 2))
    >>> t3 = Triangle(Point(0, 0), Point(3, 0), Point(0, 1))
    >>> are_similar(t1, t2)
    True
    >>> are_similar(t1, t3)
    False

    """
    from .exceptions import GeometryError

    if e1 == e2:
        return True
    try:
        return e1.is_similar(e2)
    except AttributeError:
        try:
            return e2.is_similar(e1)
        except AttributeError:
            n1 = e1.__class__.__name__
            n2 = e2.__class__.__name__
            raise GeometryError(
                "Cannot test similarity between %s and %s" % (n1, n2))


def centroid(*args):
    """Find the centroid (center of mass) of the collection containing only Points,
    Segments or Polygons. The centroid is the weighted average of the individual centroid
    where the weights are the lengths (of segments) or areas (of polygons).
    Overlapping regions will add to the weight of that region.

    If there are no objects (or a mixture of objects) then None is returned.

    See Also
    ========

    diofant.geometry.point.Point, diofant.geometry.line.Segment,
    diofant.geometry.polygon.Polygon

    Examples
    ========

    >>> p = Polygon((0, 0), (10, 0), (10, 10))
    >>> q = p.translate(0, 20)
    >>> p.centroid, q.centroid
    (Point2D(20/3, 10/3), Point2D(20/3, 70/3))
    >>> centroid(p, q)
    Point2D(20/3, 40/3)
    >>> p, q = Segment((0, 0), (2, 0)), Segment((0, 0), (2, 2))
    >>> centroid(p, q)
    Point2D(1, -sqrt(2) + 2)
    >>> centroid(Point(0, 0), Point(2, 0))
    Point2D(1, 0)

    Stacking 3 polygons on top of each other effectively triples the
    weight of that polygon:

        >>> p = Polygon((0, 0), (1, 0), (1, 1), (0, 1))
        >>> q = Polygon((1, 0), (3, 0), (3, 1), (1, 1))
        >>> centroid(p, q)
        Point2D(3/2, 1/2)
        >>> centroid(p, p, p, q) # centroid x-coord shifts left
        Point2D(11/10, 1/2)

    Stacking the squares vertically above and below p has the same
    effect:

        >>> centroid(p, p.translate(0, 1), p.translate(0, -1), q)
        Point2D(11/10, 1/2)

    """

    from .polygon import Polygon, Segment
    from .point import Point
    if args:
        if all(isinstance(g, Point) for g in args):
            c = Point(0, 0)
            for g in args:
                c += g
            den = len(args)
        elif all(isinstance(g, Segment) for g in args):
            c = Point(0, 0)
            L = 0
            for g in args:
                l = g.length
                c += g.midpoint*l
                L += l
            den = L
        elif all(isinstance(g, Polygon) for g in args):
            c = Point(0, 0)
            A = 0
            for g in args:
                a = g.area
                c += g.centroid*a
                A += a
            den = A
        c /= den
        return c.func(*[i.simplify() for i in c.args])
