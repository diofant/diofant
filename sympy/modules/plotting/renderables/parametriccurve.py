from sympy import Basic, Symbol
from sympy.modules.plotting.scene import Renderable
from util import frange, rinterpolate, interpolatecolor

from OpenGL.GL import *

class ParametricCurve(Renderable):

    vertices = []
    
    def __init__(self, x, y, z, i_t):
        t, t_min, t_max, t_steps = i_t
        t_set = frange(t_min, t_max, t_steps)
        self.vertices = [ (float(x.subs(t, t_e)), float(z.subs(t, t_e)), float(y.subs(t, t_e))) for t_e in t_set]
        self.calczmax()

    def calczmax(self):
        self.z_max, self.z_min = None, None
        for i in range(len(self.vertices)):
            if self.z_max == None:
                self.z_max = self.vertices[i][1]
                self.z_min = self.vertices[i][1]
            else:
                self.z_max = max( [self.vertices[i][1], self.z_max] )
                self.z_min = min( [self.vertices[i][1], self.z_min] )

    def render(self):
        glBegin(GL_LINE_STRIP)
        for x in range(0, len(self.vertices)):
            z_r = rinterpolate(self.z_min, self.z_max, self.vertices[x][1])
            glColor3f(*interpolatecolor((0.0,0.0,1.0),(1.0,0.0,0.0),z_r))
            glVertex3f(*self.vertices[x]);
        glEnd()
