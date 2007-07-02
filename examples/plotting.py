"""
SymPy Plotting Example
======================

The plotting module is Brian Jorgensen's project for
Google Summer of Code 2007. Please send questions and comments
to the sympy mailing list, sympy@googlegroups.com.

Keyboard Controls
=================
yaw - a,d
pitch - w,s
zoom - r,f

Tutorial
========

The following example creates a Plot containing the graphs
of x**2 and -x**2 for values of x from -5 to 5, using 10
sampling intervals:

>>> from sympy import Basic, Symbol, sqrt
>>> x = Symbol('x')
>>> f1 = x**2
>>> i1 = [x,-5,5,10]
>>> p = Plot(f1, i1, show=False) # show=True by default

You can also graph 2d and 3d graphs on the same Plot:

>>> f2 = x**2-y**2
>>> p.append(f2, [x, -5, 5, 10], [y, -5, 5, 10])

append()  and __init__() take the same argument format;
*args is parsed by function parse_args. See parse_args
for the grammar. You can specify multiple sets of
intervals at once:

>>> f3 = sin(x)+sin(y)
>>> p.append(x**2, [x,0,5,25], f3, [x,-5,5,10], [y,-5,5,10])

Then...

p.show()

...opens a new Plot window containing all of the
preceding graphs.

note: no >>> on the last line because we don't want to
doctest it. doctest +SKIP isn't available in Python 2.4,
which is needed to run PyOGL on Windows.

To use a renderable besides CartesianCurve (the default for 2d
functions) or CartesianSurface (3d), you need to bypass the
Plot interface to directly access the underlying Scene object:

p = Plot( 0.1+x**2+y**2, -0.1-x**2-y**2,
          [x, -0.75, 0.75, 8], [y, -0.75, 0.75, 8],
          show=False )
p._scene.append(ParametricCurve(sin(x), cos(x), x/5.0, [x, -10, 10, 80]))
p.show() 
"""

if __name__ == "__main__":
    print """
SymPy Plotting Example
======================

The plotting module is Brian Jorgensen's project for
Google Summer of Code 2007. Please send questions and comments
to the sympy mailing list, sympy@googlegroups.com.

Keyboard Controls
=================
yaw - a,d
pitch - w,s
zoom - r,f

Documenation and Tutorial
=========================
Please see source of plotting.py.

""" 
    import sys
    sys.path.append("..")
    
    from sympy import Symbol, sin, cos, pi
    from sympy.modules.plotting import Plot
    from sympy.modules.plotting.renderables import ParametricCurve

    x, y = Symbol('x'), Symbol('y')

    print "Now building model and opening plot window..."

    Plot( x*y**3-y*x**3, x**2+y**2+1, [x, -1, 1, 12], [y, -1, 1, 12], -x, x, [x, -1, 1, 1] )
    #Plot( x**2-y**2, [x, -1, 1, 10], [y, -1, 1, 10] )
