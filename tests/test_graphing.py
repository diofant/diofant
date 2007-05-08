import sys
sys.path.append(".")

from sympy import Symbol, Rational, sqrt, log, sin
import sympy.modules.graphing as g

x = Symbol('x')

def xtest_plot2d_one_function():
    y = x**3+x
    # plot doesn't return anything yet, but raises an exception on failure
    g.plot(y, [x, -10.0, 10.0], show=False)
   
def xtest_plot2d_list_of_functions():
    y1 = x**Rational(2)+x
    y2 = sin(x)
    # plot doesn't return anything yet, but raises an exception on failure
    g.plot([y1, y2, x], [x, -10.0, 10.0], show=False)

def xtest_plot2d_functions_outside_domain():
    y1 = sqrt(x)
    y2 = log(x)
    # Evaluating these functions for x <= 0 raises a ValueError exception.
    # Plotting should catch this exception and skip the sample point.
    g.plot([y1, y2], [x, -10.0, 10.0], show=False)
