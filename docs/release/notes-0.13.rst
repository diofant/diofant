============
Diofant 0.13
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Removed ``n()`` method from :class:`~diofant.core.evalf.EvalfMixin`, see :pull:`1114`.
* Former submodule ``diofant.polys.polyconfig`` now is :mod:`diofant.config`, see :pull:`1115`.

Minor changes
=============

Developer changes
=================

* Turn on type checking for the whole codebase, see :pull:`1114`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/7?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`20861`: reduce_inequalities() gives impossible answer
