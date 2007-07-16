from pyglet.gl import *
from plot_object import PlotObject
from plot_function import PlotFunction

class BoundingBox(PlotObject):
    
    x_min, x_max = 0.0, 0.0
    y_min, y_max = 0.0, 0.0
    z_min, z_max = 0.0, 0.0
    color = (0.7,0.75,0.8)    
    
    
    def __init__(self):
        pass

    def consider_function(self, f):
        self.x_min = min([self.x_min, f.x_min])
        self.x_max = max([self.x_max, f.x_max])
        self.y_min = min([self.y_min, f.y_min])
        self.y_max = max([self.y_max, f.y_max])
        self.z_min = min([self.z_min, f.z_min])
        self.z_max = max([self.z_max, f.z_max])

    def render(self):
        if self.x_min == self.x_max and self.y_min == self.y_max and self.z_min == self.z_max:
            return

        glColor3f(*self.color)

        glBegin(GL_LINE_STRIP)
        glVertex3f(self.x_min, self.y_max, self.z_max)
        glVertex3f(self.x_min, self.y_max, self.z_min)
        glVertex3f(self.x_max, self.y_max, self.z_min)
        glVertex3f(self.x_max, self.y_max, self.z_max)
        glVertex3f(self.x_min, self.y_max, self.z_max)
        glEnd()

        glBegin(GL_LINE_STRIP)
        glVertex3f(self.x_min, self.y_min, self.z_max)
        glVertex3f(self.x_min, self.y_min, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_max)
        glVertex3f(self.x_min, self.y_min, self.z_max)
        glEnd()
        
        glBegin(GL_LINES)
        glVertex3f(self.x_min, self.y_max, self.z_max)
        glVertex3f(self.x_min, self.y_min, self.z_max)
        
        glVertex3f(self.x_min, self.y_max, self.z_min)
        glVertex3f(self.x_min, self.y_min, self.z_min)
        
        glVertex3f(self.x_max, self.y_max, self.z_min)
        glVertex3f(self.x_max, self.y_min, self.z_min)
        
        glVertex3f(self.x_max, self.y_max, self.z_max)
        glVertex3f(self.x_max, self.y_min, self.z_max)
        glEnd()
        