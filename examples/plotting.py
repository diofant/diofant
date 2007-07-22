"""
Plotting Example

Note: In Python < 2.5, you will need the ctypes library
to use plotting. It is included with Python 2.5 and later.
"""

import sys
sys.path.append("..")

from time import sleep
from sympy import Symbol, sin, cos
from sympy.modules.plotting import Plot

x, y = Symbol('x'), Symbol('y')

if __name__ == "__main__":

    p = Plot(width = 500, height = 450,
             grid = 'xy',
             bounding_box = False,
             wireframe = False,
             #ortho = True)
             ortho = False)

    #p[1] = x, [x,-3,3,2]
    #p[2] = 1/x, [x,-3,3]
    #p[3] =  x**2, [x,-3,3]
    #p[4] = -x**2, [x,-3,3]
    #p[5] =  x**2 + y**2, [x,-1,1], [y,-1,1]
    #p[6] = -x**2 - y**2, [x,-1,1], [y,-1,1]
    #p[7] = sin(x), cos(x), [x, 0, 6.282], 'mode=parametric;visible=false'
    #p[8] = sin(x)/2, x/5.0, cos(x)/2, [x, -3.14, 3.14, 20], 'mode=parametric'
    #p[9] = x*y**3-y*x**3, [x,-1,1], [y,-1,1]
    p[10] = -1 + x**2 - y**2, [x,-1,1], [y,-1,1]
    p[11] =  1 - x**2 + y**2, [x,-1,1], [y,-1,1]
    #p[12] = x*y, [x,-1,1], [y,-1,1]
    #p[13] =  x**2 + y**2, [x,-1,1,4], [y,-1,1,4]
