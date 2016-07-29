"""
A geometry module for the Diofant library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Usage:
======


Notes:
======
    Currently the geometry module supports 2-dimensional
    and 3 -dimensional Euclidean space.

Examples
========

"""
from diofant.geometry.point import Point, Point2D, Point3D
from diofant.geometry.line import Line, Ray, Segment
from diofant.geometry.line3d import Line3D, Segment3D, Ray3D
from diofant.geometry.plane import Plane
from diofant.geometry.ellipse import Ellipse, Circle
from diofant.geometry.polygon import Polygon, RegularPolygon, Triangle, rad, deg
from diofant.geometry.util import (are_similar, centroid, convex_hull, idiff,
                                   intersection)
from diofant.geometry.exceptions import GeometryError
from diofant.geometry.curve import Curve
