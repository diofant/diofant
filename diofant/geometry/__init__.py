"""
A geometry module for the Diofant library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Notes
=====

Currently the geometry module supports 2-dimensional
and 3-dimensional Euclidean space.

"""

from .point import Point, Point2D, Point3D
from .line import Line, Ray, Segment
from .line3d import Line3D, Segment3D, Ray3D
from .plane import Plane
from .ellipse import Ellipse, Circle
from .polygon import Polygon, RegularPolygon, Triangle, rad, deg
from .util import are_similar, centroid, convex_hull, idiff, intersection
from .exceptions import GeometryError
from .curve import Curve
