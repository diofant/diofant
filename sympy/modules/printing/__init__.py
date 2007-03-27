"""Module for printing. Requires libxml2 and libxslt in order to do
transformations on the MathML, which is the underlying representation
that sympy uses

Examples
========
    >>> from sympy import *
    >>> from sympy.modules.printing import *
    
    >>> x = Symbol('x')
    >>> f = integrate(exp(x), x, evaluate=False)
    >>> print_gtk( f ) #doctest:SKIP
    
    >>> print_latex( f )
    '$\\int {e}^{x}dx$'caca

"""

from gtk import print_gtk
from latex import print_latex
from pygame_ import print_pygame
