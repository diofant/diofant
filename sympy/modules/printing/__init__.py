"""Module for printing. Requires libxml2 and libxslt in order to do
transformations on the MathML, which is the underlying representation
that sympy uses

Examples
========
    >>> from sympy import *
    >>> from sympy.modules.printing.gtk import print_gtk
    >>> from sympy.modules.printing.latex import latex
    
    >>> x = Symbol('x')
    >>> f = integrate(exp(x), x, evaluate=False)

    #>>> print_gtk( f ) 
    >>> print latex( f )
    $\int {e}^{x}\,dx$

"""

__all__ = ["mathml", "gtk", "latex", "mathml", "pretty", "pygame_"]