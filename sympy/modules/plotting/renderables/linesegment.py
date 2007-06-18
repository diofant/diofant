from sympy.modules.plotting.scene import Renderable

from OpenGL.GL import *

class LineSegment(Renderable):
    """
    Renders a line segment.
    """

    colors = [(1.0,0.0,0.0), (0.0,1.0,0.0)]

    def __init__(self, *args):
        if len(args) != 2 and len(args) != 4:
            raise ValueError("First two args must be coordinate triplets (x,y,z).")
        self.vertices = args[0:2]
        if len(args) == 4:
            self.colors = args[2:]

    def render(self):
        glBegin(GL_LINES)

        glColor3f(*self.colors[0])
        glVertex3f(*self.vertices[0])
        glColor3f(*self.colors[1])
        glVertex3f(*self.vertices[1])

        glEnd()
