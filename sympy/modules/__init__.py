"""
Description
===========

Modules of sympy. This is the place for the high-level functions 
of sympy, like limit, integrate and trigonometric functions.


Structure
=========

G{packagetree sympy.modules}
"""

from sympy.modules.trigonometric import sin, cos, tan, arctan
from sympy.modules.limits import limit, limitinf, Limit
from sympy.modules.integrals import integrate
from sympy.modules.polynomials import ispoly
from sympy.modules.paulialgebra import Pauli
from sympy.modules.matrices import Matrix, zero, one, gamma, sigma
from sympy.modules.matrices import minkowski_tensor
# printing modules are not imported by default because they depend on libxslt
