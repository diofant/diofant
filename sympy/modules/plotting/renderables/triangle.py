from sympy.modules.plotting.scene import Renderable

from OpenGL.GL import *

class Triangle(Renderable):
    """
    Renders an arbitrary triangle. Mostly for development and testing.
    """

    colors = [(1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0)]

    def __init__(self, *args):
        if len(args) != 3 and len(args) != 6:
            raise ValueError("First three args must be coordinate triplets (x,y,z).")
        self.vertices = args[0:3]
        if len(args) == 6:
            self.colors = args[3:]

    def render(self):
        glBegin(GL_TRIANGLES)

        glColor3f(*self.colors[0])
        glVertex3f(*self.vertices[0])
        glColor3f(*self.colors[1])
        glVertex3f(*self.vertices[1])
        glColor3f(*self.colors[2])
        glVertex3f(*self.vertices[2])

        glEnd()
