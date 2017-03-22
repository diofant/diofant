============
SymPy 0.5.12
============

27 Jan 2008

* SymPy works with NumPy out of the box.
* ``RootOf`` implemented.
* Lambda support works now.
* Heuristic Risch method improved.
* ``cancel()`` function implemented.
* ``sqrt(x)`` is now equivalent to ``x**(1/2)``.
* ``Derivative`` is now unevaluated.
* ``list2numpy()`` implemented.
* Series expansion of hyperbolic functions fixed.
* ``sympify('lambda x: 2*x')`` works, plus other fixes.
* Simple maxima parser implemented.
* ``sin(x)[0]`` idiom changed to ``sin(x).args[0]``
* ``sin(x).series(x, 5)`` idiom changed to ``sin(x).series(x, 0, 5)``
* Caching refactored.
* Integration of trigonometry expressions improved
* Pretty-printing for list and tuples implemented.
* Python printing implemented.
* 2D plots now don't rotate in 3D, but translate instead.
