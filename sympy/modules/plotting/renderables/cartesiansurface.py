from sympy import Basic, Symbol
from sympy.modules.plotting.scene import Renderable
from util import frange, rinterpolate, interpolatecolor

from OpenGL.GL import *

instantiation_count = 0

class CartesianSurface(Renderable):
    
    colorschemes = [ [ (0.4, 0.0, 1.0), (1.0, 0.0, 0.4) ],
                     
                     [ (0.0, 0.4, 1.0), (0.0, 1.0, 0.4) ] ]
    
    def __init__(self, f, i_x, i_y):
        x, self.x_min, self.x_max, self.x_steps = i_x
        y, self.y_min, self.y_max, self.y_steps = i_y
        self.y_set = frange(self.y_min, self.y_max, self.y_steps)
        self.x_set = frange(self.x_min, self.x_max, self.x_steps)
        self.vertices = [ [(x_e, float(f.subs(x, x_e).subs(y, y_e)), y_e) for y_e in self.y_set] for x_e in self.x_set ]
        self.calczrange()
        global instantiation_count
        self.colorscheme = instantiation_count % len(self.colorschemes)
        instantiation_count += 1

    def calczrange(self):
        self.z_max, self.z_min = None, None
        for j in range(len(self.y_set)):
            for i in range(len(self.x_set)):
                if self.z_max == None:
                    self.z_max = self.vertices[i][j][1]
                    self.z_min = self.vertices[i][j][1]
                else:
                    self.z_max = max( [self.vertices[i][j][1], self.z_max] )
                    self.z_min = min( [self.vertices[i][j][1], self.z_min] )

    def applycolor(self, x, z, y):
        """
        The y and z axis are swapped compared to the function itself.
        I eventually want to change that.
        """
        x_r = rinterpolate(self.x_min, self.x_max, x)
        y_r = rinterpolate(self.y_min, self.y_max, y)
        z_r = rinterpolate(self.z_min, self.z_max, z)

        return interpolatecolor(self.colorschemes[self.colorscheme][0],
                                self.colorschemes[self.colorscheme][1],
                                z_r)

    def render(self):
        x_len = len(self.vertices)
        y_len = len(self.vertices[0])

        for x in range(1, x_len):
            glBegin(GL_TRIANGLE_STRIP)
            for y in range(0, y_len):

                glColor3f(*self.applycolor(*self.vertices[x][y]))
                glVertex3f(*self.vertices[x][y]);

                glColor3f(*self.applycolor(*self.vertices[x-1][y]))
                glVertex3f(*self.vertices[x-1][y])

            glEnd()
