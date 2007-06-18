from sympy import Basic, Symbol
from sympy.modules.plotting.scene import Renderable
from util import frange, interpolatecolor, rinterpolate

from OpenGL.GL import *

class CartesianCurve(Renderable):

    colors = [[(1.0,0.0,0.0), (0.0,1.0,0.0)], [(0.0,0.0,1.0), (1.0,0.0,1.0)]]
    vertices = []
    
    def __init__(self, f, i_x):
        x, x_min, x_max, x_steps = i_x
        x_set = frange(x_min, x_max, x_steps)
        def eval(f, x, x_e):
            try:
                return float(f.subs(x, x_e))
            except:
                return None
        self.vertices = [ (x_e, eval(f, x, x_e), 0.0) for x_e in x_set]
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
            if self.vertices[x][1] == None:
                glEnd()
                glBegin(GL_LINE_STRIP)
                continue
            glColor3f(*interpolatecolor((0.0,0.0,1.0),(1.0,0.0,0.0),z_r))
            glVertex3f(*self.vertices[x]);
        glEnd()
