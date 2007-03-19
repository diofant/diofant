""" Some utility functions to convert sympy's objects to latex code"""

from sympy.modules.mathml import mml2latex

def print_latex(x):
    """Returns latex code representing the current object. 
    @param x: a sympy object. It can be any expression as long 
        as it inherits from basic
    @return: a string with TeX code (something like $\int \left( 1 +yx\right)dx$)
    """
    return mml2latex(x.mathml)

#TODO: print to dvi, pdf, etc.