"""
A geometry module for the Diofant library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Notes
=====

Currently the geometry module supports 2-dimensional
and 3-dimensional Euclidean space.
"""

from .point import Point, Point2D, Point3D  # noqa: F401
from .line import Line, Ray, Segment  # noqa: F401
from .line3d import Line3D, Segment3D, Ray3D  # noqa: F401
from .plane import Plane  # noqa: F401
from .ellipse import Ellipse, Circle  # noqa: F401
from .polygon import Polygon, RegularPolygon, Triangle, rad, deg  # noqa: F401
from .util import are_similar, centroid, convex_hull, idiff, intersection  # noqa: F401
from .exceptions import GeometryError  # noqa: F401
from .curve import Curve  # noqa: F401
