============
Diofant 0.15
============

Not Released Yet

New features
============

* New configuration option (``MAX_INTEGER_NBITS``) to control the maximal size of evaluated integers, see :pull:`1327`.
* Added :func:`~diofant.polys.polytools.eliminate` to eliminate symbols from the equations, see :pull:`1331`.

Major changes
=============

Compatibility breaks
====================

* Removed ``itermonomials()`` and ``topological_sort()`` functions, see :pull:`1321` and :pull:`1322`.
* Removed ``Float.num`` property, use :func:`mpmath.mpmathify`, see :pull:`1323`.

Minor changes
=============

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/9?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`25142`: incorrect simplification of a complex relational
* :sympyissue:`19813`: logcombine hangs
* :sympyissue:`22450`: Rational raised to the big power hangs
* :sympyissue:`25165`: Series expansion not working
* :sympyissue:`25197`: Simple exponential integral error in an otherwise case
* :sympyissue:`23399`: Simplifying equation with function seemingly gets stuck
* :sympyissue:`20427`: Result from clear_denoms() prints like zero poly but behaves wierdly (due to unstripped DMP)
* :sympyissue:`2720` eliminate()
* :sympyissue:`16951`: integrate(sqrt(2*m*(E - x)), x)
* :sympyissue:`25341`: CoercionFailed on eq: 2*sqrt(x)/(x + 1)**2 - 1/(sqrt(x)*(x + 1)) - 1/(4*x**(3/2)))/(x + 1) = 0
* :sympyissue:`20327`: Finite Field coercion fails from Rational type
* :sympyissue:`25406`: Resultant of Polynomials Returns Wrong Output
* :sympyissue:`25451`: Incorrect simplification when mixing basic logical operators and equality
* :sympyissue:`25496`: Privileging expr.__class__ over expr.func for reconstruction
* :sympyissue:`25521`: integrate raises HeuristicGCDFailed
