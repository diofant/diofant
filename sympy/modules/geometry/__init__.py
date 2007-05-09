"""
A geometry module for the SymPy library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.
"""
from point import Point
from line import LinearEntity, Line, Ray, Segment
from ellipse import Ellipse, Circle
from polygon import Polygon, RegularPolygon, Triangle


def intersection(e1, e2):
    from entity import GeometryEntity
    return GeometryEntity.intersection(e1, e2)