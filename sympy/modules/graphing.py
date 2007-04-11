try:
    import matplotlib.numerix as nx
    import pylab as p
except ImportError: 
    raise ImportError( "To use this module you will need matplotlib (on debian, python-matplotlib)")

from sympy.core.basic import Basic
from sympy.core.symbol import Symbol

def plot(f, var=None, ppoints=80, axis=True, show=True):
    """ works only for one variable, lacks basic options, etc...
   
    try to follow mathematica syntax for plotting: http://documents.wolfram.com/mathematica/functions/Plot
 
       example: plot(x**2, (x, 0, 10))
    """
    if var is None:
        var = (f.atoms(type=Symbol)[0], -5, 5)
    elif len(var) != 3:
        raise ValueError("second vargument must be of the form (x, x_min, x_max)")
    t1 = nx.arange(float(var[1]), float(var[2]), (float(var[2])-float(var[1]))/ppoints)
    f_t1 = []
    for num in t1:
        f_t1.append(float(f.subs(var[0], Basic.sympify(num))))
    p.plot(t1, f_t1)
    p.title(str(f))
    p.grid(True)
    if show: 
        p.show()
    #TODO: return something
