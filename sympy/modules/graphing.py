"""
This module supports function graphing. Currently 2D graphing
is supported, and 3D is planned under Brian Jorgensen's Summer
of Code project. 2D graphing requires matplotlib.

Example
=======

>>> from sympy import Symbol
>>> x = Symbol('x')
>>> y = x**2+x

#>>> plot(y, [x, -10.0, 10.0], show=False)
"""

from sympy import Symbol, Basic

try:
    import pylab
except ImportError:
    raise ImportError("To use this module you will need matplotlib (on debian, python-matplotlib)")

def plot(f, var=None, plot_points=100, axis=True, show=True, grid=True, title=None, xlabel=None, ylabel=None):
    """
    Tries to be similar to mathematica syntax:
    http://documents.wolfram.com/mathematica/functions/Plot
    
    'f' can be a single function or a list of functions.
    
    'var' should be in the format (x, x_min, x_max) where x is a
    sympy Symbol representing the independent variable in f.

    plot([sqrt(x), log(x)], [x, -10.0, 10.0]) #doctest: +SKIP
    """

    try:
        len(f)
    except (TypeError):
        f = [f]

    if len(f[0].atoms(type=Symbol)) == 0:
        raise ValueError("First argument must be an expression of one variable.")

    if len(f[0].atoms(type=Symbol)) > 1:
        # plot3d
        raise ValueError("First argument must be an expression of one variable. 3d graphing not yet supported.")

    #plot2d
    v, v_min, v_max = None, None, None

    if var is None:
        v, v_min, v_max = f.atoms(type=Symbol)[0], -10.0, 10.0
    elif len(var) == 3:
        v, v_min, v_max = var[0], float(var[1]), float(var[2])
    else:
        raise ValueError("Second argument must be in the form (x, x_min, x_max)")
    
    plot_points = float(plot_points)

    def plot_f(f, v, v_min, v_max, plot_points):
        delta = (v_max - v_min) / plot_points
        x_a = pylab.arange(v_min, v_max, delta)
        y_a = []
        for x in x_a:
            try:
                y_i = float( f.subs(v, Basic.sympify(x)) )
            except (OverflowError, ValueError):
                y_i = None # f(x) is undefined or otherwise unplottable
            y_a.append(y_i)

        pylab.plot(x_a, y_a)

    for fx in f:
        plot_f(fx, v, v_min, v_max, plot_points)

    if title == None:
        title = ", ".join([str(fx) for fx in f])
    pylab.title(title)
    
    if xlabel != None: pylab.xlabel(xlabel)
    if ylabel != None: pylab.ylabel(ylabel)
    
    pylab.grid(grid)
    pylab.draw()

    if show:
        pylab.show()

    return
    return True

