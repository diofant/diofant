============
SymPy 0.5.13
============

6 Mar 2008

* SymPy is now 2x faster in average compared to the previous release

  - first patches with 25% speedup
  - ``Basic.cos`` et. al. removed, use ``C.cos`` instead
  - ``sympy.core`` now uses direct imports
  - ``sympifyit`` decorator
  - speedup Integers creation and arithmetic
  - speedup unary operations for singleton numbers
  - remove silly slowdowns from fast-path of mul and div
  - significant speedup was achieved by reusing dummy variables
  - ``is_dummy`` is not an assumption anymore
  - Symbols & Wilds are cached
  - ``((2+3*I)**1000).expand()`` is now at least 100x faster
  - ``.expand()`` was made faster for cases where an expression is already expanded
  - rational powers of integers are now computed more efficiently
  - unknown assumptions are now cached as well as known assumptions

* ``integrate()`` can handle most of the basic integrals now
* interactive experience with isympy was improved through adding support for , () and {} to pretty-printer, and switching to it as the default ipython printer
* new ``trim()`` function to map all non-atomic expressions, ie. functions, derivatives and more complex objects, to symbols and remove common factors from numerator and denominator. also cancel() was improved
* ``.expand()`` for noncommutative symbols fixed
* bug in ``(x+y+sin(x)).as_independent()`` fixed
* ``.subs_dict()`` improved
* support for plotting geometry objects added
* bug in ``.tangent_line()`` of ellipse fixed
* new atan2 function and assotiated fixes for ``.arg()`` and expanding rational powers
* new ``.coeff()`` method for returning coefficient of a poly
* pretty-printer now uses unicode by default
* recognition of geometric sums were generalized
* ``.is_positive`` and ``.is_negative`` now fallback to ``evalf()`` when appropriate
* as the result ``oo*(pi-1)`` now correctly simplifies to oo
* support for objects which provide ``__int__`` method was added
* we finally started SymPy User's Guide
* ``BasicMeths`` merged into ``Basic``
* cache subsystem was cleaned up -- now it supports only immutable objects
