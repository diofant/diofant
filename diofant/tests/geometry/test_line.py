import itertools

import pytest

from diofant import (Float, I, Rational, Symbol, cos, oo, pi, simplify, sin,
                     sqrt, symbols, tan)
from diofant.geometry import (Circle, GeometryError, Line, Point, Ray, Segment,
                              Triangle, intersection)
from diofant.geometry.line import Undecidable
from diofant.geometry.polygon import _asa as asa


__all__ = ()

x = Symbol('x', extended_real=True)
y = Symbol('y', extended_real=True)
z = Symbol('z', extended_real=True)
t = Symbol('t', extended_real=True)
k = Symbol('k', extended_real=True)
x1 = Symbol('x1', extended_real=True)
x2 = Symbol('x2', extended_real=True)
x3 = Symbol('x3', extended_real=True)
y1 = Symbol('y1', extended_real=True)
y2 = Symbol('y2', extended_real=True)
y3 = Symbol('y3', extended_real=True)
z1 = Symbol('z1', extended_real=True)
z2 = Symbol('z2', extended_real=True)
z3 = Symbol('z3', extended_real=True)
half = Rational(1, 2)


def feq(a, b):
    """Test if two floating point values are 'equal'."""
    t = Float('1.0E-10')
    return -t < a - b < t


def test_line_geom():
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(x1, x1)
    p4 = Point(y1, y1)
    p5 = Point(x1, 1 + x1)
    p6 = Point(1, 0)
    p7 = Point(0, 1)
    p8 = Point(2, 0)
    p9 = Point(2, 1)

    l1 = Line(p1, p2)
    l2 = Line(p3, p4)
    l3 = Line(p3, p5)
    l4 = Line(p1, p6)
    l5 = Line(p1, p7)
    l6 = Line(p8, p9)
    l7 = Line(p2, p9)
    pytest.raises(ValueError, lambda: Line(Point(0, 0), Point(0, 0)))

    # Basic stuff
    assert Line((1, 1), slope=1) == Line((1, 1), (2, 2))
    assert Line((1, 1), slope=oo) == Line((1, 1), (1, 2))
    assert Line((1, 1), slope=-oo) == Line((1, 1), (1, 2))
    pytest.raises(ValueError, lambda: Line((1, 1), 1))
    assert Line(p1, p2) == Line(p1, p2)
    assert Line(p1, p2) != Line(p2, p1)
    assert l1 != l2
    assert l1 != l3
    assert l1.slope == 1
    assert l1.length == oo
    assert l3.slope == oo
    assert l4.slope == 0
    assert l4.coefficients == (0, 1, 0)
    assert l4.equation(x=x, y=y) == y
    assert l5.slope == oo
    assert l5.coefficients == (1, 0, 0)
    assert l5.equation() == x
    assert l6.equation() == x - 2
    assert l7.equation() == y - 1
    assert p1 in l1  # is p1 on the line l1?
    assert p1 not in l3
    assert Line((-x, x), (-x + 1, x - 1)).coefficients == (1, 1, 0)

    assert simplify(l1.equation()) in (x - y, y - x)
    assert simplify(l3.equation()) in (x - x1, x1 - x)

    assert Line(p1, p2).scale(2, 1) == Line(p1, p9)

    assert l2.arbitrary_point() in l2
    for ind in range(5):
        assert l3.random_point() in l3
    pytest.raises(ValueError, lambda: l3.arbitrary_point('x1'))

    assert Line(Point(0, 0), Point(1, 0)).is_similar(Line(Point(1, 0),
                                                          Point(2, 0))) is True

    assert l1.equal(l1) is True
    assert l1.equal(l2) is True
    assert l1.equal(l3) is False
    assert l1.equal(object()) is False

    # Orthogonality
    p1_1 = Point(-x1, x1)
    l1_1 = Line(p1, p1_1)
    assert l1.perpendicular_line(p1.args) == Line(Point(0, 0), Point(1, -1))
    assert l1.perpendicular_line(p1) == Line(Point(0, 0), Point(1, -1))
    assert Line.is_perpendicular(l1, l1_1)
    assert Line.is_perpendicular(l1, l2) is False
    p = l1.random_point()
    assert l1.perpendicular_segment(p) == p

    assert l4.perpendicular_line(p2) == Line(Point(1, 1), Point(1, 0))

    # Parallelity
    l2_1 = Line(p3, p5)
    assert l2.parallel_line(p1_1) == Line(Point(-x1, x1), Point(-y1, 2*x1 - y1))
    assert l2_1.parallel_line(p1.args) == Line(Point(0, 0), Point(0, -1))
    assert l2_1.parallel_line(p1) == Line(Point(0, 0), Point(0, -1))
    assert Line.is_parallel(l1, l2)
    assert Line.is_parallel(l2, l3) is False
    assert Line.is_parallel(l2, l2.parallel_line(p1_1))
    assert Line.is_parallel(l2_1, l2_1.parallel_line(p1))

    # Intersection
    assert intersection(l1, p1) == [p1]
    assert intersection(l1, p5) == []
    assert intersection(l1, l2) in [[l1], [l2]]
    assert intersection(l1, l1.parallel_line(p5)) == []

    # Concurrency
    l3_1 = Line(Point(5, x1), Point(-Rational(3, 5), x1))
    assert Line.are_concurrent(l1) is False
    assert Line.are_concurrent(l1, l3)
    assert Line.are_concurrent(l1, l3, l3_1)
    assert Line.are_concurrent(l1, l1_1, l3) is False

    # Projection
    assert l2.projection(p4) == p4
    assert l1.projection(p1_1) == p1
    assert l3.projection(p2) == Point(x1, 1)
    pytest.raises(GeometryError, lambda: Line(Point(0, 0), Point(1, 0))
                  .projection(Circle(Point(0, 0), 1)))

    # Finding angles
    l1_1 = Line(p1, Point(5, 0))
    assert feq(Line.angle_between(l1, l1_1).evalf(), pi.evalf()/4)

    # Testing Rays and Segments (very similar to Lines)
    pytest.raises(ValueError, lambda: Ray((1, 1), I))
    pytest.raises(ValueError, lambda: Ray(p1, p1))
    pytest.raises(ValueError, lambda: Ray(p1))
    assert Ray((1, 1), angle=pi/4) == Ray((1, 1), (2, 2))
    assert Ray((1, 1), angle=pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=-pi/2) == Ray((1, 1), (1, 0))
    assert Ray((1, 1), angle=-3*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=5*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=5.0*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=pi) == Ray((1, 1), (0, 1))
    assert Ray((1, 1), angle=3.0*pi) == Ray((1, 1), (0, 1))
    assert Ray((1, 1), angle=4.0*pi) == Ray((1, 1), (2, 1))
    assert Ray((1, 1), angle=0) == Ray((1, 1), (2, 1))
    assert Ray((1, 1), angle=4.05*pi) == Ray(Point(1, 1),
                                             Point(2, -sqrt(5)*sqrt(2*sqrt(5) + 10)/4 - sqrt(2*sqrt(5) + 10)/4 + 2 + sqrt(5)))
    assert Ray((1, 1), angle=4.02*pi) == Ray(Point(1, 1),
                                             Point(2, 1 + tan(4.02*pi)))
    assert Ray((1, 1), angle=5) == Ray((1, 1), (2, 1 + tan(5)))
    pytest.raises(ValueError, lambda: Ray((1, 1), 1))

    # issue sympy/sympy#7963
    r = Ray((0, 0), angle=x)
    assert r.subs({x: 3*pi/4}) == Ray((0, 0), (-1, 1))
    assert r.subs({x: 5*pi/4}) == Ray((0, 0), (-1, -1))
    assert r.subs({x: -pi/4}) == Ray((0, 0), (1, -1))
    assert r.subs({x: pi/2}) == Ray((0, 0), (0, 1))
    assert r.subs({x: -pi/2}) == Ray((0, 0), (0, -1))

    r1 = Ray(p1, Point(-1, 5))
    r2 = Ray(p1, Point(-1, 1))
    r3 = Ray(p3, p5)
    r4 = Ray(p1, p2)
    r5 = Ray(p2, p1)
    r6 = Ray(Point(0, 1), Point(1, 2))
    r7 = Ray(Point(0.5, 0.5), Point(1, 1))
    assert l1.projection(r1) == Ray(Point(0, 0), Point(2, 2))
    assert l1.projection(r2) == p1
    assert r3 != r1
    t = Symbol('t', extended_real=True)
    assert Ray((1, 1), angle=pi/4).arbitrary_point() == \
        Point(t + 1, t + 1)
    r8 = Ray(Point(0, 0), Point(0, 4))
    r9 = Ray(Point(0, 1), Point(0, -1))
    assert r8.intersection(r9) == [Segment(Point(0, 0), Point(0, 1))]

    s1 = Segment(p1, p2)
    s2 = Segment(p1, p1_1)
    assert s1.midpoint == Point(Rational(1, 2), Rational(1, 2))
    assert s2.length == sqrt( 2*(x1**2) )
    assert Segment((1, 1), (2, 3)).arbitrary_point() == Point(1 + t, 1 + 2*t)
    assert s1.perpendicular_bisector() == \
        Line(Point(1/2, 1/2), Point(3/2, -1/2))
    # intersections
    assert s1.intersection(Line(p6, p9)) == []
    s3 = Segment(Point(0.25, 0.25), Point(0.5, 0.5))
    assert s1.intersection(s3) == [s1]
    assert s3.intersection(s1) == [s3]
    assert r4.intersection(s3) == [s3]
    assert r4.intersection(Segment(Point(2, 3), Point(3, 4))) == []
    assert r4.intersection(Segment(Point(-1, -1), Point(0.5, 0.5))) == \
        [Segment(p1, Point(0.5, 0.5))]
    s3 = Segment(Point(1, 1), Point(2, 2))
    assert s1.intersection(s3) == [Point(1, 1)]
    s3 = Segment(Point(0.5, 0.5), Point(1.5, 1.5))
    assert s1.intersection(s3) == [Segment(Point(0.5, 0.5), p2)]
    assert s1.intersection(Segment(Point(4, 4), Point(5, 5))) == []
    assert s1.intersection(Segment(Point(-1, -1), p1)) == [p1]
    assert s1.intersection(Segment(Point(-1, -1), Point(0.5, 0.5))) == \
        [Segment(p1, Point(0.5, 0.5))]
    assert r4.intersection(r5) == [s1]
    assert r5.intersection(r6) == []
    assert r4.intersection(r7) == r7.intersection(r4) == [r7]

    # Segment contains
    a, b = symbols('a,b')
    s = Segment((0, a), (0, b))
    assert Point(0, (a + b)/2) in s
    s = Segment((a, 0), (b, 0))
    assert Point((a + b)/2, 0) in s

    pytest.raises(Undecidable, lambda: Point(2*a, 0) in s)

    # Testing distance from a Segment to an object
    s1 = Segment(Point(0, 0), Point(1, 1))
    s2 = Segment(Point(half, half), Point(1, 0))
    pt1 = Point(0, 0)
    pt2 = Point(Rational(3, 2), Rational(3, 2))
    assert s1.distance(pt1) == 0
    assert s1.distance((0, 0)) == 0
    assert s2.distance(pt1) == 2**half/2
    assert s2.distance(pt2) == 2**half
    # Line to point
    p1, p2 = Point(0, 0), Point(1, 1)
    s = Line(p1, p2)
    assert s.distance(Point(-1, 1)) == sqrt(2)
    assert s.distance(Point(1, -1)) == sqrt(2)
    assert s.distance(Point(2, 2)) == 0
    assert s.distance((-1, 1)) == sqrt(2)
    assert Line((0, 0), (0, 1)).distance(p1) == 0
    assert Line((0, 0), (0, 1)).distance(p2) == 1
    assert Line((0, 0), (1, 0)).distance(p1) == 0
    assert Line((0, 0), (1, 0)).distance(p2) == 1
    m = symbols('m')
    l = Line((0, 5), slope=m)
    p = Point(2, 3)
    assert l.distance(p) == 2*abs(m + 1)/sqrt(m**2 + 1)
    # Ray to point
    r = Ray(p1, p2)
    assert r.distance(Point(-1, -1)) == sqrt(2)
    assert r.distance(Point(1, 1)) == 0
    assert r.distance(Point(-1, 1)) == sqrt(2)
    assert Ray((1, 1), (2, 2)).distance(Point(1.5, 3)) == 3*sqrt(2)/4
    assert r.distance((1, 1)) == 0
    assert r.distance((-1, Rational(1, 2))) == sqrt(5)/2

    # Line contains
    p1, p2 = Point(0, 1), Point(3, 4)
    l = Line(p1, p2)
    assert l.contains(p1) is True
    assert l.contains((0, 1)) is True
    assert l.contains((0, 0)) is False
    assert l.contains(Circle(p1, 1)) is False
    assert l.contains(Ray(p1, p1 + p2)) is False

    # Ray contains
    p1, p2 = Point(0, 0), Point(4, 4)
    r = Ray(p1, p2)
    assert r.contains(p1) is True
    assert r.contains((1, 1)) is True
    assert r.contains((1, 3)) is False
    assert r.contains(object()) is False
    s = Segment((1, 1), (2, 2))
    assert r.contains(s) is True
    s = Segment((1, 2), (2, 5))
    assert r.contains(s) is False
    r1 = Ray((2, 2), (3, 3))
    assert r.contains(r1) is True
    r1 = Ray((2, 2), (3, 5))
    assert r.contains(r1) is False
    r1 = Ray(p1, angle=-pi)
    assert r1.contains(Point(1, 0)) is False
    r1 = Ray(p1, angle=-pi/2)
    assert r1.contains(Point(0, 1)) is False
    pytest.raises(Undecidable, lambda: r1.contains(Point(0, x)))

    # Special cases of projection and intersection
    r1 = Ray(Point(1, 1), Point(2, 2))
    r2 = Ray(Point(2, 2), Point(0, 0))
    r3 = Ray(Point(1, 1), Point(-1, -1))
    r4 = Ray(Point(0, 4), Point(-1, -5))
    r5 = Ray(Point(2, 2), Point(3, 3))
    assert intersection(r1, r2) == [Segment(Point(1, 1), Point(2, 2))]
    assert intersection(r1, r3) == [Point(1, 1)]
    assert r1.projection(r3) == Point(1, 1)
    assert r1.projection(r4) == Segment(Point(1, 1), Point(2, 2))

    r5 = Ray(Point(0, 0), Point(0, 1))
    r6 = Ray(Point(0, 0), Point(0, 2))
    assert r5 in r6
    assert r6 in r5

    s1 = Segment(Point(0, 0), Point(2, 2))
    s2 = Segment(Point(-1, 5), Point(-5, -10))
    s3 = Segment(Point(0, 4), Point(-2, 2))
    assert intersection(r1, s1) == [Segment(Point(1, 1), Point(2, 2))]
    assert r1.projection(s2) == Segment(Point(1, 1), Point(2, 2))
    assert s3.projection(r1) == Segment(Point(0, 4), Point(-1, 3))

    l1 = Line(Point(0, 0), Point(3, 4))
    r1 = Ray(Point(0, 0), Point(3, 4))
    s1 = Segment(Point(0, 0), Point(3, 4))
    assert intersection(l1, l1) == [l1]
    assert intersection(l1, r1) == [r1]
    assert intersection(l1, s1) == [s1]
    assert intersection(r1, l1) == [r1]
    assert intersection(s1, l1) == [s1]

    entity1 = Segment(Point(-10, 10), Point(10, 10))
    entity2 = Segment(Point(-5, -5), Point(-5, 5))
    assert intersection(entity1, entity2) == []

    r1 = Ray(p1, Point(0, 1))
    r2 = Ray(Point(0, 1), p1)
    r3 = Ray(p1, p2)
    r4 = Ray(p2, p1)
    s1 = Segment(p1, Point(0, 1))
    assert Line(r1.source, r1.random_point()).slope == r1.slope
    assert Line(r2.source, r2.random_point()).slope == r2.slope
    assert Segment(Point(0, -1), s1.random_point()).slope == s1.slope
    p_r3 = r3.random_point()
    p_r4 = r4.random_point()
    assert p_r3.x >= p1.x and p_r3.y >= p1.y
    assert p_r4.x <= p2.x and p_r4.y <= p2.y
    p10 = Point(2000, 2000)
    s1 = Segment(p1, p10)
    p_s1 = s1.random_point()
    assert p1.x <= p_s1.x and p_s1.x <= p10.x and \
        p1.y <= p_s1.y and p_s1.y <= p10.y
    s2 = Segment(p10, p1)
    assert hash(s1) == hash(s2)
    p11 = p10.scale(2, 2)
    assert s1.is_similar(Segment(p10, p11))
    assert s1.is_similar(r1) is False
    assert (r1 in s1) is False
    assert Segment(p1, p2) in s1
    assert s1.plot_interval() == [t, 0, 1]
    assert s1 in Line(p1, p10)
    assert Line(p1, p10) != Line(p10, p1)
    assert Line(p1, p10) != p1
    assert Line(p1, p10).plot_interval() == [t, -5, 5]
    assert Ray((0, 0), angle=pi/4).plot_interval() == \
        [t, 0, 10]

    p1, p2 = Point(0, 0), Point(4, 1)
    r1 = Ray(p1, p2)
    assert r1.direction == p2

    p1, p2, p3 = Point(0, 0), Point(-1, -1), Point(-1, 0)
    r1, r2 = Ray(p1, p2), Ray(p1, p3)
    assert r1.ydirection == -oo
    assert r2.ydirection == 0

    p1, p2, p3 = Point(0, 0), Point(6, 6), Point(5, 1)
    s1 = Segment(p1, p2)
    assert s1.perpendicular_bisector() == Line(Point(3, 3), Point(9, -3))
    assert s1.perpendicular_bisector(p3) == Segment(Point(3, 3), Point(5, 1))


def test_line_intersection():
    assert asa(120, 8, 52) == \
        Triangle(
            Point(0, 0),
            Point(8, 0),
            Point(-4*cos(19*pi/90)/sin(2*pi/45),
                  4*sqrt(3)*cos(19*pi/90)/sin(2*pi/45)))
    assert Line((0, 0), (1, 1)).intersection(Ray((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Line((0, 0), (1, 1)).intersection(Segment((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Ray((0, 0), (1, 1)).intersection(Ray((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Ray((0, 0), (1, 1)).intersection(Segment((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert (Ray((0, 0), angle=-pi).intersection(Segment((-1, 0), (2, 0))) ==
            [Segment((-1, 0), (0, 0))])
    assert Ray((0, 0), (10, 10)).contains(Segment((1, 1), (2, 2))) is True
    assert Segment((1, 1), (2, 2)) in Line((0, 0), (10, 10))
    x = 8*tan(13*pi/45)/(tan(13*pi/45) + sqrt(3))
    y = (-8*sqrt(3)*tan(13*pi/45)**2 + 24*tan(13*pi/45)) / \
        (-3 + tan(13*pi/45)**2)
    assert Line(Point(0, 0), Point(1, -sqrt(3))).contains(Point(x, y)) is True

    # issue sympy/sympy#2941
    def _check():
        for f, g in itertools.product(*[(Line, Ray, Segment)]*2):
            l1 = f(a, b)
            l2 = g(c, d)
            assert l1.intersection(l2) == l2.intersection(l1)
    # intersect at end point
    c, d = (-2, -2), (-2, 0)
    a, b = (0, 0), (1, 1)
    _check()
    # midline intersection
    c, d = (-2, -3), (-2, 0)
    a, b = (0, 0), (1, 1)
    _check()


def test_symbolic_intersection():
    # Issue sympy/sympy#7814.
    circle = Circle(Point(x, 0), y)
    line = Line(Point(k, z), slope=0)
    assert line.intersection(circle) == [
        Point(x - sqrt(y**2 - z**2), z), Point(x + sqrt(y**2 - z**2), z)]
