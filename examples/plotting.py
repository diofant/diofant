import sys
sys.path.append("..")

from time import sleep
from sympy import Symbol, sin, cos
from sympy.modules.plotting import Plot

x, y = Symbol('x'), Symbol('y')

if __name__ == "__main__":

    p = Plot(bbox = True)

    #p[1] = x, [x,-3,3,2]
    #p[2] = 1/x, [x,-3,3]
    #p[3] =  x**2, [x,-3,3]
    #p[4] = -x**2, [x,-3,3]
    #p[5] =  x**2 + y**2, [x,-1,1], [y,-1,1]
    #p[6] = -x**2 - y**2, [x,-1,1], [y,-1,1]
    #p[7] = sin(x), cos(x), [x, 0, 6.282], 'type=parametric;visible=false'
    #p[8] = sin(x), x/5.0, cos(x), [x, -3.14*4, 3.14*4], 'type=parametric'
    #p[9] = x*y**3-y*x**3, [x,-1,1], [y,-1,1]
    p[10] = x**2 - y**2, [x,-1,1], [y,-1,1]

    print p
