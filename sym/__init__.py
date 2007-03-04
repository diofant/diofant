"""Computer Algebra System (CAS) in Python

Sympy is a symbolic manipulation package, written in pure Python. 
Its aim is to become a full featured CAS in Python, while keeping 
the code as simple as possible in order to be comprehensible and 
easily extensible.

To help you get started, here is a simple example in the python
interpreter:

>>> import sym
>>> x = sym.Symbol('x')
>>> e = 1/sym.cos(x)
>>> print e.series(x,10)
1+1/2*x^2+5/24*x^4+61/720*x^6+277/8064*x^8+50521/3628800*x^10

for full documentation, see the docs at our web page:
http://code.google.com/p/sympy/wiki/Documentation
"""

from core import *
from modules import *
