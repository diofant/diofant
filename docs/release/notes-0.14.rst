============
Diofant 0.14
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Removed support for CPython 3.9, see :pull:`1192`.
* Removed ``to_mpi()`` method of :class:`~diofant.sets.sets.Interval`, see :pull:`1194`.

Minor changes
=============

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/8?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`22487`: [integrals] Wrong result for Integral((cos(x**2)-cos(x))/x**2, (x, -oo, oo))
* :sympyissue:`22493`: Series expansion introduces new variables
* :sympyissue:`22558`: Error in ODE-Solver-Documentation
* :sympyissue:`22837`: Solve simplest algebraic equations with dummy parameter
* :sympyissue:`22836`: Series: Possible improvements for Order of expressions involving factorials
* :sympyissue:`22788`: RecursionError for unevluated expression in latex
* :sympyissue:`22863`: Hangs: integrate((3*x**3-x**2+2*x-4)/sqrt(x**2-3*x+2), (x, 0, 1))
* :sympyissue:`22862`: Problem with separable differential equation
* :sympyissue:`22893`: 'limit' in combination with 'positive=True' gives wrong result
* :sympyissue:`22878`: RecursionError in trigsimp
* :sympyissue:`22982`: limit((log(E + 1/x) - 1)**(1 - sqrt(E + 1/x)), x, oo) returns 0 instead of oo
* :sympyissue:`22986`: limit(acosh(1 + 1/x)*sqrt(x), x, oo) is evaluated incorrectly.
