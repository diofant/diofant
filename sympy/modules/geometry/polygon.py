from sympy import Basic, Rational, Symbol
from sympy.modules.simplify import simplify
from entity import GeometryEntity
from point import Point
from ellipse import Circle
from line import LinearEntity, Line, Segment


class Polygon(GeometryEntity):
    """A polygon in space."""
    def __init__(self, *args, **kwargs):
        GeometryEntity.__init__(self, **kwargs)
        if not isinstance(args[0], GeometryEntity):
            self._vertices = tuple(args[0])
        else:
            self._vertices = tuple(args)

        for p in self._vertices:
            if not isinstance(p, Point):
                raise TypeError("__init__ requires Point instancess")

        if len(self._vertices) < 3:
            raise ValueError("A polygon must have at least three points")

    @property
    def area(self):
        """Returns the area of the polygon."""
        area = 0
        for ind in xrange(0, len(self._vertices)-1):
            pi = self._vertices[ind]
            pii = self._vertices[ind+1]
            area += pi[0]*pii[1]-pii[0]*pi[1]
        return simplify(area) / Rational(2)

    @property
    def perimeter(self):
        """Returns the perimeter of the polygon."""
        p = sum([side.length for side in self.sides])
        return simplify(p)

    @property
    def vertices(self):
        """Returns the vertices that define the polygon."""
        return self._vertices

    @property
    def centroid(self):
        """Returns the centroid of the polygon."""
        # TODO Double check this formula
        x = Rational(0)
        y = Rational(0)
        n = Rational(len(self._points))
        for p in self._vertices:
            x += p[0]
            y += p[1]
        x = simplify(x/n)
        y = simplify(y/n)
        return Point(x, y)

    @property
    def sides(self):
        """Returns a list of the sides."""
        res = []
        for ind in xrange(1, len(self._vertices)):
            res.append( Segment(self._vertices[ind-1], self._vertices[ind]) )
        res.append( Segment(self._vertices[-1], self._vertices[0]) )
        return res

    def intersection(self, o):
        """
        Returns the intersection of the Polygon and another entity, or None if
        there is no intersection.
        """
        res = []
        for side in self.sides:
            inter = GeometryEntity.intersection(side, o)
            if inter is not None:
                res.extend(inter)

        if len(res) == 0:
            return None
        return res

    def __len__(self):
        return len(self._vertices)

    def __getitem__(self, ind):
        return self._vertices[ind]

    def __eq__(self, o):
        if not isinstance(o, Polygon):
            return False

        # Find index of the first point that is the same
        n1,n2 = len(self._vertices), len(o._vertices)
        start_indices = []
        for ind in xrange(0, n1):
            if self._vertices[ind] == o._vertices[0]:
                start_indices.append(ind)

        if len(start_indices) == 0:
            return False

        # Check vertices
        imax = max(n1, n2)
        for start_ind in start_indices:
            i = start_ind

            # Check to see what orientation we should check
            dir = 0
            if self._vertices[(i + 1) % n1] == o._vertices[1]:
                dir = 1
            elif self._vertices[(i - 1) % n1] == o._vertices[1]:
                dir = -1

            # If either point to the left or right if the first point
            # is value (i.e., dir is nonzero) then check in that direction
            if dir != 0:
                areEqual = True
                for ind in xrange(2, imax):
                    if self._vertices[(i + dir*ind) % n1] != o._vertices[ind % n2]:
                        areEqual = False
                        break
                if areEqual: return True

        return False

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(self._vertices)

    def __contains__(self, o):
        if isinstance(o, Polygon):
            return self == o
        elif isinstance(o, Segment):
            return o in self.sides
        elif isinstance(o, Point):
            if o in self._vertices:
                return True
            for side in self.sides:
                if o in side:
                    return True
            return False
        else:
            return False

    def __str__(self):
        #what_am_i = {
        #    3: "Triangle",
        #    4: "Quadrilateral",
        #    5: "Pentagon",
        #    6: "Hexagon",
        #    7: "Heptagon",
        #    8: "Octagon",
        #    9: "Nonagon",
        #    10: "Decagon"
        #}.get(len(self._vertices), "Polygon")
        return "Polygon(%d sides)" % len(self._vertices)

    def __repr__(self):
        return str(self)

class RegularPolygon(Polygon):
    def __init__(self, *args, **kwargs):
        Polygon.__init__(self, *args, **kwargs)

class Triangle(Polygon):
    @staticmethod
    def are_similar(t1, t2):
        """
        Returns the True if triangles t1 and t2 are similar,
        False otherwise.
        """
        s1_1, s1_2, s1_3 = [side.length for side in t1.sides]
        s2 = [side.length for side in t2.sides]
        #print [s1_1, s1_2, s1_3] + s2
        def _are_similar(u1, u2, u3, v1, v2, v3):
            e1 = simplify(u1/v1)
            e2 = simplify(u2/v2)
            e3 = simplify(u3/v3)
            return bool(e1 == e2) and bool(e2 == e3)

        # There's only 6 permutations, so write them out
        return _are_similar(s1_1, s1_2, s1_3, *s2) or \
               _are_similar(s1_1, s1_3, s1_2, *s2) or \
               _are_similar(s1_2, s1_1, s1_3, *s2) or \
               _are_similar(s1_2, s1_3, s1_1, *s2) or \
               _are_similar(s1_3, s1_1, s1_2, *s2) or \
               _are_similar(s1_3, s1_2, s1_1, *s2)

    def is_equilateral(self):
        """Returns True if the triangle is equilateral, False otherwise."""
        s = self.sides
        return bool(s[0].length == s[1].length) and bool(s[1].length == s[2].length)

    def is_right(self):
        """Returns True if the triangle is right-angled, False otherwise."""
        #for angle in self.angles:
        #    if angle == pi/2: return True
        #return False
        s = self.sides
        return LinearEntity.are_perpendicular(s[0], s[1]) or \
               LinearEntity.are_perpendicular(s[1], s[2]) or \
               LinearEntity.are_perpendicular(s[0], s[2])

    @property
    def altitudes(self):
        """
        Returns the altitudes of the triangle in a dictionary where the key
        is the vertex and the value is the altitude at that point.
        """
        # XXX Is this abusing the fact that we know how self.sides
        #     constructs its side, or shall we consider the way
        #     self.sides is constructed is standard?
        s = self.sides
        v = self._vertices
        return {v[0]: s[1].perpendicular_segment(v[0]),
                v[1]: s[2].perpendicular_segment(v[1]),
                v[2]: s[0].perpendicular_segment(v[2])}

    @property
    def orthocenter(self):
        """Returns the orthocenter of the triangle."""
        a = self.altitudes
        return GeometryEntity.intersect(a[1], a[2])[0]

    @property
    def circumcenter(self):
        return self.orthocenter

    @property
    def circumradius(self):
        return Point.distance(self.circumcenter, self._vertices[0])

    @property
    def circumcircle(self):
        return Circle(self.circumcenter, self.circumradius)

    @property
    def bisectors(self):
        """
        Returns the bisectors of the triangle in a dictionary where the key
        is the vertex and the value is the bisector at that point.
        """
        s = self.sides
        v = self._vertices
        c = self.incenter
        l1 = Segment(v[0], GeometryEntity.intersection(Line(v[0], c), s[1])[0])
        l2 = Segment(v[1], GeometryEntity.intersection(Line(v[1], c), s[2])[0])
        l3 = Segment(v[2], GeometryEntity.intersection(Line(v[2], c), s[0])[0])
        return {v[0]: l1, v[1]: l2, v[2]: l3}

    @property
    def incenter(self):
        """Returns the incenter of the triangle."""
        s = self.sides
        v = self._vertices
        A,B,C = v[0],v[1],v[2]
        a,b,c = s[1].length,s[2].length,s[0].length
        x = simplify( (a*A[0] + b*B[0] + c*C[0]) / (a+b+c) )
        y = simplify( (a*A[1] + b*B[1] + c*C[1]) / (a+b+c) )
        return Point(x, y)

    @property
    def inradius(self):
        """Returns the inradius of the triangle."""
        return simplify(self.area / self.perimeter)

    @property
    def incircle(self):
        """Returns the incircle of the triangle."""
        return Circle(self.incenter, self.inradius)

    @property
    def medians(self):
        """
        Returns the medians of the triangle in a dictionary where the key
        is the vertex and the value is the median at that point.
        """
        # XXX See Triangle.altitudes for comments on the usage of self.sides 
        s = self.sides
        v = self._vertices
        return {v[0]: Segment(s[1].midpoint, v[0]),
                v[1]: Segment(s[2].midpoint, v[1]),
                v[2]: Segment(s[0].midpoint, v[2])}

    @property
    def medial(self):
        """Returns the medial triangle of the triangle."""
        s = self.sides
        return Triangle(s[0].midpoint, s[1].midpoint, s[2].midpoint)

    @property
    def excircles(self):
        """
        Returns a list of the three excircles for this triangle.
        """
        pass

    def __str__(self):
        fmt_tuple = (str(self._vertices[0]), str(self._vertices[1]), str(self._vertices[2]))
        return "Triangle(%s, %s, %s)" % fmt_tuple