import sys
import os
import threading
from time import sleep

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

class Renderable(object):
    """
    Base class for objects which can be rendered by a Scene.
    See sympy\modules\plotting\renderables\ for examples.
    """
    def __init__(self):
        raise Exception("Renderable base class can not be instantiated.")

    def render(self):
        raise Exception("Renderable child class must override render.")

class Scene(object):
    """
    Provides base OpenGL functionality. It handles the GLUT window
    and the rendering loop. The rendering loop operates by drawing
    each element in the renderables list. To add objects to the
    scene, use the renderable object list accessor functions below.
    """
    
    def __init__(self, title="Untitled Scene", width=800, height=600, fullscreen=False, dist=5.0):
        """
        Initializes the plot. The window and OpenGL initialization
        is not done until show() is called.
        """
        self.title = title
        self.width = width
        self.height = height
        self.fullscreen = fullscreen
        
        self._window = None
        self._renderables = []

        self.dist = dist
        self.yrot = 0.0
        self.xrot = 0.0


    """
    Renderable object list accessors
    """
    def clear(self): self._renderables = []
    def append(self, renderable):
        assert isinstance(renderable, Renderable)
        self._renderables.append(renderable)
    def __len__(self): return len(self._renderables)
    def __iter__(self): return iter(self._renderables)
    def __getitem__(self, i): return self._renderables[i]
    def __setitem__(self, i, v):
        assert isinstance(v, Renderable)
        self._renderables[i] = v
    def __delitem__(self, i): del self._renderables[i] 

    def glutloop(self):
        """
        Create the GLUT window and start the rendering loop.
        """
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        glutInitWindowSize(self.width, self.height)
        glutInitWindowPosition(0, 0)
        self._window = glutCreateWindow(self.title)
        glutDisplayFunc(self._onrender)
        if self.fullscreen:
            glutFullScreen()
        glutIdleFunc(self._onrender)
        glutReshapeFunc(self._onresize)
        glutKeyboardFunc(self._onkeypress)
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glClearDepth(1.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glShadeModel(GL_SMOOTH)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(self.width)/float(self.height), 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)
        glutMainLoop()
        
    def show(self):
        #while glutGetWindow() != 0: pass
        #threading.Thread(target=self.glutloop).start()
        self.glutloop()
    
    def _onrender(self):
        """
        GLUT rendering loop callback
        """
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        
        glTranslate(0.0, 0.0, -self.dist)
        #gluLookAt(0.0,0.0, self.dist, 0.0,0.0,0.0, 0.0,1.0,0.0);
        glRotate(self.xrot, 1.0, 0.0, 0.0);
        glRotate(self.yrot, 0.0, 1.0, 0.0);

        for r in self._renderables:
            glPushMatrix()
            r.render()
            glPopMatrix()

        glutSwapBuffers()    
    
    def _onkeypress(self, *args):
        """
        GLUT key press callback
        """
        ESCAPE = '\033'
        if args[0] == ESCAPE:
            #glutDestroyWindow(self._window) #doesn't work yet
            sys.exit()
        # zoom in
        elif args[0] == 'r':
            if self.dist > 1.0:
                self.dist -= 0.2
        # zoom out
        elif args[0] == 'f':
            self.dist += 0.2
        # yaw
        elif args[0] == 'a':
            self.yrot -= 3.0
        elif args[0] == 'd':
            self.yrot += 3.0
        # pitch
        elif args[0] == 'w':
            self.xrot -= 3.0
        elif args[0] == 's':
            self.xrot += 3.0
    
    def _onresize(self, width, height):
        """
        GLUT resize callback
        """
        self.width = width
        self.height = height
        if self.height == 0:
            self.height = 1
        glViewport(0, 0, self.width, self.height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(self.width)/float(self.height), 0.1, 500.0)
        glMatrixMode(GL_MODELVIEW)
