===========
SymPy 0.6.1
===========

22 Jul 2008

* almost all functions and constants can be converted to Sage
* univariate factorization algorithm was fixed
* ``.evalf()`` method fixed, ``pi.evalf(106)`` calculates 1 000 000 digits of pi
* ``@threaded`` decorator
* more robust solvers, polynomials and simplification
* better simplify, that makes a solver more robust
* optional compiling of functions to machine code
* msolve: solving of nonlinear equation systems using Newton's method
* ``((x+y+z)**50).expand()`` is now 3 times faster
* caching was removed from the Order class: 1.5x speedups in series tests
