from pyglet.gl import *
from window_manager import ManagedWindow

from plot_camera import PlotCamera
from plot_controller import PlotController

class PlotWindow(ManagedWindow):

    _calculating = False

    def __init__(self, plot,
                 title="SymPy Plot",
                 wireframe=False,
                 antialiasing=True,
                 ortho=False,
                 **kwargs):
        self.plot = plot
        self.title = title
        self.wireframe = wireframe
        self.antialiasing = antialiasing
        self.ortho = ortho

        kwargs['caption'] = title
        super(PlotWindow, self).__init__(**kwargs)

    def setup(self):
        self.camera = PlotCamera(self.window)

        self.controller = PlotController(self)
        self.window.push_handlers(self.controller)

        glClearColor(1.0,1.0,1.0,0.0)
        glClearDepth(1.0)

        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)

        glShadeModel(GL_SMOOTH)

        if self.wireframe:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        else:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

        if self.antialiasing:
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
            glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)

        self.init_projection()

    def init_projection(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.ortho:
            raise NotImplementedError("Orthographic projection not implemented.")
        else:
            gluPerspective(60.0, float(self.window.width)/float(self.window.height), 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def update(self, dt):
        self.update_caption()
        return self.controller.update(dt)

    def draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        self.plot.lock_begin()
        self.camera.apply()
        for r in self.plot._plotobjects:
            glPushMatrix()
            r.render()
            glPopMatrix()
        for r in self.plot._functions.itervalues():
            glPushMatrix()
            r.render()
            glPopMatrix()
        self.plot.lock_end()

    def update_caption(self):
        if not self._calculating and self.plot._calculations_in_progress > 0:
            self.window.set_caption(self.title + " (calculating...)")
            self._calculating = True

        elif self._calculating and self.plot._calculations_in_progress == 0:
            self.window.set_caption(self.title)
            self._calculating = False
