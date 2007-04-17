""" Some utility functions to convert sympy's objects to latex code"""

from sympy.modules.mathml import apply_xsl, c2p

def print_latex(x):
    """
    Usage
    =====
        Returns latex code representing the current object. 
    
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> from sympy.modules.printing import print_latex
        >>> print print_latex( integrate(x*y-2, x, evaluate=False))
        $\\int \\left(-2+xy\\right)dx$
        
    Parameters
    ==========
        @param x: a sympy object. It can be any expression as long 
            as it inherits from basic
        @return: a string with TeX code (something like $\int \left( 1 +yx\right)dx$)
    """
    return apply_xsl(c2p(x.mathml), 'mathml/data/mmltex.xsl')

#TODO: print to dvi, pdf, etc.