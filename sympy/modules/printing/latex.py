""" Some utility functions to convert sympy's objects to latex code"""

from sympy.core.basic import Basic

def latex(x):
    """
    Usage
    =====
        Returns the latex code representing the current object. 
    
    Notes
    =====
        @param x: a sympy object. It can be any expression as long 
            as it inherits from basic
        @return: a string with TeX code (something like $\int \left( 1 +yx\right)dx$)
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> from sympy.modules.printing.latex import latex
        >>> print latex( integrate(x*y-2, x, evaluate=False))
        $\int -2+x y\,dx$
        
    """
    x = Basic.sympify(x)
    return "$" + x.__latex__() + "$"

def print_latex(x):
    print latex(x)
    return
#TODO: print to dvi, pdf, etc.