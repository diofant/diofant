from pyglet.gl import *
from plot_rotation import get_spherical_rotatation

class PlotCamera(object):

    minimum_dist = 0.0

    _dist = 5.0
    _rot = None

    def __init__(self, window):
        self.window = window
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

    def apply(self):
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -self._dist)
        if self._rot != None:
            glMultMatrixf(self._rot)

    def spherical_rotate(self, p1, p2, sensitivity=1.0):
        mat = get_spherical_rotatation(p1, p2, self.window.width, self.window.height, sensitivity)
        if mat != None: self.mult_rot_matrix(mat)

    def zoom_relative(self, clicks, sensitivity):
        if clicks < 0 or (self._dist - clicks*sensitivity) >= self.minimum_dist:
            self._dist -= clicks * sensitivity
