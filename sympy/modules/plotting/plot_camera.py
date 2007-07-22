from pyglet.gl import *
from plot_rotation import get_spherical_rotatation

class PlotCamera(object):

    min_dist = 1.0
    max_dist = 500.0

    min_ortho_dist = 100.0
    max_ortho_dist = 10000.0

    _dist = 10.0
    _rot = None
    
    _default_ortho_dist = 1000.0

    def __init__(self, window, ortho = False):
        self.window = window
        self.ortho = ortho
        if self.ortho:
            self._dist = self._default_ortho_dist
        self.init_rot_matrix()

    def get_matrix(self):
        m = (c_float*16)()
        glGetFloatv(GL_MODELVIEW_MATRIX, m)
        return m

    def init_rot_matrix(self):
        glPushMatrix()
        glLoadIdentity()
        self._rot = self.get_matrix()
        glPopMatrix()

    def mult_rot_matrix(self, rot):
        glPushMatrix()
        glLoadIdentity()
        glMultMatrixf(rot)
        glMultMatrixf(self._rot)
        self._rot = self.get_matrix()
        glPopMatrix()

    def setup_projection(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.ortho:
            # yep, pseudo ortho (don't tell anyone)
            gluPerspective(0.3, float(self.window.width) / float(self.window.height), self.min_ortho_dist-0.5, self.max_ortho_dist+0.5)
        else:
            gluPerspective(30.0, float(self.window.width) / float(self.window.height), self.min_dist-0.5, self.max_dist+0.5)
        glMatrixMode(GL_MODELVIEW)

    def apply_transformation(self):
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -self._dist)
        if self._rot != None:
            glMultMatrixf(self._rot)

    def spherical_rotate(self, p1, p2, sensitivity=1.0):
        mat = get_spherical_rotatation(p1, p2, self.window.width, self.window.height, sensitivity)
        if mat != None: self.mult_rot_matrix(mat)

    def zoom_relative(self, clicks, sensitivity):
        
        if self.ortho:
            dist_d = clicks * sensitivity * 50.0
            min_dist = self.min_ortho_dist
            max_dist = self.max_ortho_dist
        else:
            dist_d = clicks * sensitivity
            min_dist = self.min_dist
            max_dist = self.max_dist
        
        new_dist = (self._dist - dist_d)
        if (clicks < 0 and new_dist < max_dist) or new_dist > min_dist:
            self._dist = new_dist
