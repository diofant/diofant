from pyglet.gl import *
from pyglet.window import Window

from threading import Thread, Lock, Event
from time import clock

global gl_lock
gl_lock = Lock()

class ManagedWindow(Window):

    fps_limit = 60
    win_args = dict(width=400,
                    height=300,
                    vsync=False,
                    resizable=False)

    def __init__(self, **win_args):
        self.win_args = dict(self.win_args, **win_args)
        self.close_later = False
        Thread(target=self.__event_loop__).start()

    def __event_loop__(self, **win_args):
        gl_lock.acquire()
        try:
            super(ManagedWindow, self).__init__(**self.win_args)
            self.switch_to()
            self.setup()
        except Exception, e:
            print "Error: %s" % str(e)
        finally:
            gl_lock.release()

        # can't use pyglet's clock because
        # we need one per event loop
        frame_duration = 1.0/self.fps_limit

        then = clock()
        while not self.has_exit and self.context != None:
            dt = 0.0
            while dt < frame_duration:
                now = clock()
                dt = now-then
            then = now

            gl_lock.acquire()
            try:
                self.switch_to()
                self.dispatch_events()

                if self.context != None:
                    self.clear()
                    self.update(dt)
                    if self.context != None:
                        self.draw()
                        self.flip()
                    if self.close_later:
                        super(ManagedWindow, self).close()

            except Exception, e:
                print "Event Loop Error: %s" % str(e)
            finally:
                gl_lock.release()

    def close(self):
        self.close_later = True

    def setup(self):
        pass

    def update(self, dt):
        pass

    def draw(self):
        pass

