from sympy.modules.plotting.scene import Renderable
from linesegment import LineSegment

from OpenGL.GL import *

class Simple3dAxes(Renderable):
    def __init__(self):
        self.axes = []
        self.axes.append(LineSegment((0, 0, 0), (1, 0, 0), (1.0, 0.0, 0.0), (1.0, 0.0, 0.0) ))
        self.axes.append(LineSegment((0, 0, 0), (0, 1, 0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0) ))
        self.axes.append(LineSegment((0, 0, 0), (0, 0, 1), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0) ))

    def render(self):
        for a in self.axes:
            a.render()
