from sympy import Basic, Symbol
from sympy.modules.plotting.scene import Renderable
from util import frange, rinterpolate, interpolatecolor
from math import sin, cos # faster, used with floats, not Basic

from OpenGL.GL import *

instantiation_count = 0

class CylindricalSurface(Renderable):
    
    colorschemes = [ [ (0.4, 0.0, 1.0), (1.0, 0.0, 0.4) ],
                     
                     [ (0.0, 0.4, 1.0), (0.0, 1.0, 0.4) ] ]
    
    def __init__(self, f, i_height, i_theta):
        f = Basic.sympify(f)
        t, self.t_min, self.t_max, self.t_steps = i_theta
        h, self.h_min, self.h_max, self.h_steps = i_height
        self.h_set = frange(self.h_min, self.h_max, self.h_steps)
        self.t_set = frange(self.t_min, self.t_max, self.t_steps)
        def eval(f, t, h, t_e, h_e):
            try:
                r = float(f.subs(t, t_e).subs(h, h_e))
                return r*cos(t_e), h_e, r*sin(t_e)
            except:
                return None, None, None
        self.vertices = [ [eval(f, t,h, t_e,h_e) for h_e in self.h_set] for t_e in self.t_set ]
        global instantiation_count
        self.colorscheme = instantiation_count % len(self.colorschemes)
        instantiation_count += 1

    def applycolor(self, x, z, y):
        h_r = rinterpolate(self.h_min, self.h_max, z)

        return interpolatecolor(self.colorschemes[self.colorscheme][0],
                                self.colorschemes[self.colorscheme][1],
                                h_r)

    def render(self):
        t_len = len(self.vertices)
        h_len = len(self.vertices[0])

        for t in range(1, t_len):
            glBegin(GL_TRIANGLE_STRIP)
            for h in range(0, h_len):

                if (self.vertices[t][h][1] == None) or (self.vertices[t-1][h][1] == None):
                    glEnd()
                    glBegin(GL_TRIANGLE_STRIP)
                    continue
                
                glColor3f(*self.applycolor(*self.vertices[t][h]))
                glVertex3f(*self.vertices[t][h]);

                glColor3f(*self.applycolor(*self.vertices[t-1][h]))
                glVertex3f(*self.vertices[t-1][h])

            glEnd()
