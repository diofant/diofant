============
Diofant 0.15
============

Not Released Yet

New features
============

* New configuration option (``MAX_INTEGER_NBITS``) to control the maximal size of evaluated integers, see :pull:`1327`.

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
