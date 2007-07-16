from threading import Thread, Lock

from pyglet.gl import *
from pyglet import clock
from pyglet.window import Window

from plot_camera import PlotCamera
from plot_controller import PlotController

class PlotWindow(object):

    def __init__(self, plot, title="SymPy Plot", width=400, height=300, fullscreen=False, wireframe=False):
        self.plot = plot
        self.title = title
        self.width = width
        self.height = height
        self.fullscreen = fullscreen
        self.wireframe = wireframe

        self._calculating = False
        self.notify_close = False

        self.thread = Thread(target=self.eventloop)
        self.thread.start()
        #self.thread.run()

    def close(self):
        self.notify_close = True

    def render(self, dt):
        self.plot.lock_begin()

        if not self._calculating and self.plot._calculations_in_progress > 0:
            self._window.set_caption(self.title + " (calculating...)")
            self._calculating = True
        elif self._calculating and self.plot._calculations_in_progress == 0:
            self._window.set_caption(self.title)
            self._calculating = False

        self._window.switch_to()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        self.controller.update(dt)
        self.camera.apply()

        for r in self.plot._plotobjects:
            if not r.visible:
                continue
            glPushMatrix()
            r.render()
            glPopMatrix()

        for r in self.plot._functions.itervalues():
            if not r.visible:
                continue
            glPushMatrix()
            r.render()
            glPopMatrix()

        self.plot.lock_end()

    def initprojection(self, ortho=False):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if ortho:
            raise NotImplementedError("Orthographic projection not implemented.")
        else:
            gluPerspective(60.0, float(self.width)/float(self.height), 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def on_resize(self, width, height):
        self.width, self.height = width, height
        self.initprojection()

    def initgl(self):
        self._window = Window(self.width, self.height, caption=self.title,
                              fullscreen=self.fullscreen)#, resizable=True)

        self.camera = PlotCamera(self)
        self.controller = PlotController(self)
        
        self._window.push_handlers(self.on_resize)
        self._window.push_handlers(self.controller)

        glClearColor(1.0,1.0,1.0,0.0)
        glClearDepth(1.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        if self.wireframe: glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        else: glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        glShadeModel(GL_SMOOTH)

        self.initprojection()

    def eventloop(self):
        self.initgl()
        clock.set_fps_limit(40)
        while not self._window.has_exit:
            if self.notify_close:
                self._window.close()
                self.notify_close = False
            else:
                dt = clock.tick()
                self._window.dispatch_events()
                self.render(dt)
                self._window.flip()
        return
