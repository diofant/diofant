from sympy import Basic, Symbol
from sympy.modules.plotting.scene import Renderable
from util import frange, rinterpolate, interpolatecolor
from math import sin, cos # faster, used with floats, not Basic

from OpenGL.GL import *

instantiation_count = 0

class SphericalSurface(Renderable):
    
    colorschemes = [ [ (0.4, 0.0, 1.0), (1.0, 0.0, 0.4) ],
                     
                     [ (0.0, 0.4, 1.0), (0.0, 1.0, 0.4) ] ]
    
    def __init__(self, f, i_theta, i_phi):
        f = Basic.sympify(f)
        t, self.t_min, self.t_max, self.t_steps = i_theta
        p, self.p_min, self.p_max, self.p_steps = i_phi
        self.t_set = frange(self.t_min, self.t_max, self.t_steps)
        self.p_set = frange(self.p_min, self.p_max, self.p_steps)
        def eval(f, t, p, t_e, p_e):
            try:
                r = float(f.subs(t, t_e).subs(p, p_e))
                return r*sin(p_e)*cos(t_e), r*cos(p_e), r*sin(p_e)*sin(t_e)
            except:
                return None, None, None
        self.vertices = [ [eval(f, t,p, t_e,p_e) for p_e in self.p_set] for t_e in self.t_set ]
        global instantiation_count
        self.colorscheme = instantiation_count % len(self.colorschemes)
        instantiation_count += 1

    def applycolor(self, x, z, y):
        p_r = rinterpolate(self.p_min, self.p_max, z)

        return interpolatecolor(self.colorschemes[self.colorscheme][0],
                                self.colorschemes[self.colorscheme][1],
                                p_r)

    def render(self):
        t_len = len(self.vertices)
        p_len = len(self.vertices[0])

        for t in range(1, t_len):
            glBegin(GL_TRIANGLE_STRIP)
            for p in range(0, p_len):

                if (self.vertices[t][p][1] == None) or (self.vertices[t-1][p][1] == None):
                    glEnd()
                    glBegin(GL_TRIANGLE_STRIP)
                    continue
                
                glColor3f(*self.applycolor(*self.vertices[t][p]))
                glVertex3f(*self.vertices[t][p]);

                glColor3f(*self.applycolor(*self.vertices[t-1][p]))
                glVertex3f(*self.vertices[t-1][p])

            glEnd()
