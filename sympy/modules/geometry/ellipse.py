from sympy import Basic, Rational, Symbol, sqrt, cos, sin, pi
from sympy.modules.simplify import simplify
from entity import GeometryEntity
from point import Point
from line import LinearEntity, Line

class Ellipse(GeometryEntity):
    """An ellipse in space."""

    def __init__(self, center, horizontal_radius, vertical_radius, **kwargs):
        GeometryEntity.__init__(self, **kwargs)
        if not isinstance(center, Point):
            raise TypeError("center must be be a Point")
        self._c = center
        # XXX Check for zero radii?
        self._hr = horizontal_radius
        self._vr = vertical_radius

    @property
    def horizontal_radius(self):
        """Returns the horizontal radius of the ellipse."""
        return self._hr

    @property
    def vertical_radius(self):
        """Returns the vertical radius of the ellipse."""
        return self._vr

    @property
    def center(self):
        """Returns the center of the ellipse."""
        return self._c

    @property
    def area(self):
        """Returns the area of the ellipse."""
        return simplify(pi * self._hr * self._vr)

    @property
    def circumference(self):
        """Returns the circumference of the ellipse."""
        # TODO It's fairly complicated, but we could use Ramanujan's approximation.
        raise NotImplementedError

    @property
    def foci(self):
        """Returns the foci of the ellipse."""
        pass

    def tangent_line(self, p):
        """
        If p is on the ellipse, returns the tangent line through point p.
        Otherwise, returns the tangent line(s) from p to the ellipse, or
        None if no tangent line is possible (e.g., p inside ellipse).
        """
        if p in self:
            rise = (self._hr ** Rational(2))*(self._c[0] - p[0])
            run = (self._vr ** Rational(2))*(p[1] - self._c[1])
            p2 = Point(simplify(p[0] + run),
                       simplify(p[1] + rise))
            return Line(p, p2)
        else:
            # TODO If p is not on the ellipse, attempt to create the
            #      tangent(s) from point p to the ellipse.
            return None

    def is_tangent(self, o):
        """Returns True if o is tangent to the ellipse, False otherwise."""
        inter = None
        if isinstance(o, Ellipse):
            inter = self.intersection(o)
            return (inter is not None) and isinstance(inter[0], Point)
        elif isinstance(o, LinearEntity):
            inter = self._do_line_intersection(o)
            if (inter is not None) and len(inter) == 1:
                return (inter[0] in o)
            else:
                return False
        else:
            raise NotImplementedError("Ellipse.is_tangent requires a Line or Ellipse instance")

    def arbitrary_point(self, parameter_name='t'):
        """Returns a symbolic point that is on the ellipse."""
        t = Symbol(parameter_name)
        return Point(self._c[0] + self._hr*cos(t), self._c[1] + self._vr*sin(t))

    def random_point(self):
        """Returns a random point on the ellipse."""
        from random import randint
        from sys import maxint
        t = Symbol('t')
        p = self.arbitrary_point('t')
        subs_val = randint(-maxint-1, maxint)
        return Point(p[0].subs(t, subs_val).eval(), p[1].subs(t, subs_val).eval())

    #@property
    def equation(self, xaxis_name='x', yaxis_name='y'):
        """Returns the equation of the ellipse."""
        x = Symbol(xaxis_name)
        y = Symbol(yaxis_name)
        t1 = ((x - self._c[0]) / self._hr)**Rational(2)
        t2 = ((y - self._c[1]) / self._vr)**Rational(2)
        return t1 + t2 - Rational(1)

    def _do_line_intersection(self, o):
        """
        Find the intersection of a LinearEntity and the ellipse. Makes no regards
        to what the LinearEntity is because it assumes a Line. To ensure correct
        intersection results one must invoke intersection() to remove bad results.
        """
        def dot(p1, p2):
            sum = Rational(0)
            for ind in xrange(0, len(p1)):
                sum += simplify(p1[ind] * p2[ind])
            return sum

        hr_sq = self._hr ** Rational(2)
        vr_sq = self._vr ** Rational(2)
        lp = o.points

        ldir = lp[1] - lp[0]
        diff = lp[0] - self._c

        mdir = (ldir[0] / hr_sq, ldir[1] / vr_sq)
        mdiff = (diff[0] / hr_sq, diff[1] / vr_sq)

        a = dot(ldir, mdir)
        b = dot(ldir, mdiff)
        c = dot(diff, mdiff) - Rational(1)
        det = simplify(b*b - a*c);

        result = []
        if det == 0:
            t = -b / a
            i1 = lp[0] + (lp[0] - lp[1]) * t;
            result.append(i1)
        else:
            is_good = True
            try:
                is_good = (det < 0)
            except NotImplementedError: #symbolic, allow
                is_good = True

            if is_good:
                root = sqrt(det)
                t_a = (-b - root) / a
                t_b = (-b + root) / a
                try:
                    if 0 <= t_a <= 1:
                        result.append( lp[0] + (lp[0] - lp[1]) * t_a )
                except NotImplementedError: #symbolic, allow
                    result.append( lp[0] + (lp[0] - lp[1]) * t_a )

                try:
                    if 0 <= t_a <= 1:
                        result.append( lp[0] + (lp[0] - lp[1]) * t_b )
                except NotImplementedError: #symbolic, allow
                    result.append( lp[0] + (lp[0] - lp[1]) * t_b )

        if len(result) == 0:
            return None
        return result

    def intersection(self, o):
        """
        Returns the intersection of the ellipse and another entity, or None if
        there is no intersection.
        """
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return None
        elif isinstance(o, LinearEntity):
            # LinearEntity may be a ray/segment, so check the points
            # of intersection for coincidence first
            result = self._do_line_intersection(o)
            for ind in xrange(0, len(result)):
                if p not in o:
                    del result[ind]
            return result
        elif isinstance(o, Ellipse):
            if o == self:
                return self
            else:
                # TODO
                raise NotImplementedError("Unable to find intersection with " + n)
        else:
            n = type(o).__name__
            raise NotImplementedError("Unable to find intersection with " + n)

    def __eq__(self, o):
        return ((self._c == o._c) and
                (self._hr == o._hr) and (self._vr == o._vr))

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((self._c, self._hr, self._vr))

    def __contains__(self, o):
        if isinstance(o, Point):
            x = Symbol('x')
            y = Symbol('y')
            res = self.equation('x', 'y').subs_dict({x: o[0], y: o[1]})
            return (simplify(res) == 0)
        elif isinstance(o, Ellipse):
            return (self == o)
        else:
            return False

    def __str__(self):
        args = (str(self._c), str(self._hr), str(self._vr))
        return "Ellipse(center=%s, horizontal radius=%s, vertical radius=%s)" % args

    def __repr__(self):
        return self.__str__()

class Circle(Ellipse):
    """A circle in space."""

    def __init__(self, *args, **kwargs):
        if len(args) == 2:
            # Assume (center, radius) pair
            Ellipse.__init__(self, args[0], args[1], args[1], **kwargs)
        elif len(args) == 3:
            # Assume three points
            raise NotImplementedError()

    @property
    def radius(self):
        return self._hr

    @property
    def circumference(self):
        return simplify(Rational(2) * pi * self.radius)

    def equation(self, xaxis_name='x', yaxis_name='y'):
        """Returns the equation of the circle."""
        x = Symbol(xaxis_name)
        y = Symbol(yaxis_name)
        t1 = (x - self._c[0])**Rational(2)
        t2 = (y - self._c[1])**Rational(2)
        return t1 + t2 - self._hr**Rational(2)

    def intersection(self, o):
        try:
            return Ellipse.intersection(self, o)
        except NotImplementedError:
            if isinstance(o, Circle):
                dx,dy = o._c - self._c
                d = sqrt( simplify(dy**2 + dx**2) )
                a = simplify( (self._hr**2 - o._hr**2 + d**2) / (Rational(2) * d) )

                x2 = simplify(self._c[0] + (dx * a/d))
                y2 = simplify(self._c[1] + (dy * a/d))

                h = sqrt( simplify(r0**2 - a**2) )
                rx = simplify(-dy * (h/d))
                ry = simplify(dx * (h/d))

                xi_1 = simplify(x2 + rx)
                xi_2 = simplify(x2 - rx)
                yi_1 = simplify(y2 + ry)
                yi_2 = simplify(y2 - ry)

                ret = [Point(xi_1, yi_1)]
                if xi_1 != xi_2 or yi_1 != yi_2:
                    ret.append(Point(xi_2, yi_2))
                return ret
            elif isinstance(o, Ellipse):
                a,b,r = o.horizontal_radius,o.vertical_radius,self._hr
                t1 = sqrt(simplify((r**2 - b**2)/(a**2 - b**2)))
                t2 = sqrt(simplify((a**2 - r**2)/(a**2 - b**2)))
                return [simplify(a*t1),
                        simplify(b*t1),
                        simplify(-a*t1),
                        simplify(-b*t1)]

        raise NotImplementedError("Unable to find intersection with " + n)