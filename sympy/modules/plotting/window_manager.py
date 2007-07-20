import sys
sys.path.append('.')

from threading import Event, Thread, Lock

from pyglet.gl import *
from pyglet import clock
from pyglet.window import Window

from time import sleep

class Singleton(object):
    def __new__(cls, *p, **k):
        if not '_the_instance' in cls.__dict__:
            cls._the_instance = object.__new__(cls)
        return cls._the_instance

class WindowManager(Singleton):

    # maximum framerate
    _framerate = 60

    def append(self, window):
        """
        Add a ManagedWindow.
        """
        self._lock.acquire()
        try:
            if window not in self._windows:
                self._windows.append(window)
                if self._kt == None or not self._kt.isAlive():
                    self.do_thread()
        finally:
            self._lock.release()

    # window list
    _windows = []
    # window list lock
    _lock = Lock()
    # kernel thread
    _kt = None

    def do_thread(self):
        self._kt = WindowManager.KernelThread()
        self._kt._windows = self._windows
        self._kt._lock = self._lock
        self._kt._framerate = self._framerate
        self._kt.start()

    class KernelThread(Thread):
        _windows = None
        _lock = None
        _framerate = None

        def update_windows(self, dt):
            self._lock.acquire()

            to_be_removed = []
            try:
                for w in self._windows:
                    if not w._do_frame(dt):
                        to_be_removed.append(w)
            finally:
                for w in to_be_removed:
                    if w in self._windows:
                        self._windows.remove(w)
                r = len(self._windows) > 0

                self._lock.release()
                return r

        def run(self):
            assert self._windows != None
            assert self._lock != None
            assert self._framerate != None

            clock.set_fps_limit(self._framerate)
            while self.update_windows( clock.tick() ):
                pass

class ManagedWindow(object):
    window = None

    default_window_args = dict(width=400,
                               height=300,
                               vsync=False,
                               resizable=False)

    def __init__(self, **window_args):
        """
        You'd be better off overriding setup instead
        in most cases.
        """
        self.window_args = dict(self.default_window_args, **window_args)
        WindowManager().append(self)

    def get_closed(self):
        return self.window == None or self.window.has_exit

    def close(self):
        """
        Mark for removal and closure.
        """
        if not self.get_closed():
            self.window.has_exit = True
            self.window.set_visible(False)

    def _do_frame(self, dt):
        """
        Called during the WindowManager's event loop.
        Do not override or modify; instead, override
        setup, update, and draw.
        """
        if self.window == None:
            self.window = Window(**self.window_args)
            self.setup()

        self.window.dispatch_events()
        if self.get_closed():
            return False
        
        self.window.switch_to()

        if not self.update(dt):
            self.close()
        else:
            self.draw()
            self.window.flip()
        return True

    def setup(self):
        """
        Override this function instead of __init__.
        """
        glClearColor(1, 1, 1, 1)

    def update(self, dt):
        """
        For logic to be performed directly before
        rendering.
        """
        pass

    def draw(self):
        """
        Rendering code.
        """
        glClear(GL_COLOR_BUFFER_BIT)


