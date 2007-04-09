import matplotlib.numerix as nx
import pylab as p
from sympy.core.basic import Basic

def plot(f, var, iter=80):
    """ works only for one variable, lacks basic options, etc...
   
    try to follow mathematica syntax for plotting: http://documents.wolfram.com/mathematica/functions/Plot
 
       example: plot(x**2, (x, 0, 10))
    """
    if len(var) != 3:
        raise ValueError("second vargument must be of the form (x, x_min, x_max)")
    t1 = nx.arange(float(var[1]), float(var[2]), (float(var[2])-float(var[1]))/iter)
    f_t1 = []
    for num in t1:
        f_t1.append(float(f.subs(var[0], Basic.sympify(num))))
    p.plot(t1, f_t1)
    p.title(str(f))
    p.show()
