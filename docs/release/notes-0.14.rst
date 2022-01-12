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
