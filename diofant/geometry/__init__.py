"""
A geometry module for the Diofant library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Notes
=====

Currently the geometry module supports 2-dimensional Euclidean space.

"""

from .curve import Curve
from .ellipse import Circle, Ellipse
from .exceptions import GeometryError
from .line import Line, Ray, Segment
from .point import Point
from .polygon import Polygon, RegularPolygon, Triangle, deg, rad
from .util import are_similar, centroid, convex_hull, idiff, intersection


__all__ = ('Curve', 'Circle', 'Ellipse', 'GeometryError', 'Line', 'Ray',
           'Segment', 'Point', 'Polygon', 'RegularPolygon',
           'Triangle', 'deg', 'rad', 'are_similar', 'centroid', 'convex_hull',
           'idiff', 'intersection')
