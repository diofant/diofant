===========
SymPy 0.7.6
===========

21 Nov 2014

New features
============

* New module ``calculus.finite_diff`` for generating finite difference formulae approximating derivatives of arbitrary order on arbitrarily spaced grids.

* New module ``physics.optics`` for symbolic computations related to optics.

* ``geometry`` module now supports 3D geometry.

* Support for series expansions at a point other then ``0`` or ``oo``.  See :sympypull:`2427`.

* Rules for the intersection of integer ImageSets were added. See :sympypull:`7587`. We can now do things like ``{2⋅m | m ∊ ℤ} ∩ {3⋅n | n ∊ ℤ} = {6⋅t | t ∊ ℤ}`` and ``{2⋅m | m ∊ ℤ} ∩ {2⋅n + 1 | n ∊ ℤ} = ∅``.

* ``dsolve`` module now supports system of ODEs including linear system of ODEs of 1st order for 2 and 3 equations and of 2nd order for 2 equations. It also supports homogeneous linear system of ``n`` equations.

* New support for continued fractions, including iterators for partial quotients and convergents, and reducing a continued fraction to a Rational or a quadratic irrational.

* Support for Egyptian fraction expansions, using several different algorithms.

* Addition of generalized linearization methods to ``physics.mechanics``.

* Use an LRU cache by default instead of an unbounded one. See :sympypull:`7464`.  Control cache size by the environment variable ``SYMPY_CACHE_SIZE`` (default is 500).  ``SYMPY_CACHE_SIZE=None`` restores the unbounded cache.

* Added ``fastcache`` as an optional dependency.  Requires v0.4 or newer.  Control via ``SYMPY_CACHE_SIZE``.  May result in significant speedup. See :sympypull:`7737`.

* New experimental module ``physics.unitsystems`` for computation with dimensions, units and quantities gathered into systems. This opens the way to dimensional analysis and better quantity calculus.  The old module ``physics.units`` will stay available until the new one reaches a mature state. See :sympypull:`2628`.

* New ``Complement`` class to represent relative complements of two sets. See :sympypull:`7462`.

* New trigonometric functions (asec, acsc), many enhancements for other trigonometric functions (see :sympypull:`7500`).

* New ``Contains`` class to represent the relation "is an element of" (see :sympypull:`7989`).

* The code generation tools (code printers, ``codegen``, ``autowrap``, and ``ufuncify``) have been updated to support a wider variety of constructs, and do so in a more robust way. Major changes include added support for matrices as inputs/outputs, and more robust handling of conditional (``Piecewise``) statements.

* ``ufuncify`` now uses a backend that generates actual ``numpy.ufuncs`` by default through the use of the ``numpy`` C api. This allows broadcasting on *all* arguments. The previous ``cython`` and ``f2py`` backends are still accessible through the use of the ``backend`` kwarg.

* ``CodeGen`` now generates code for Octave and Matlab from SymPy expressions.  This is supported by a new CodePrinter with interface ``octave_code``.  For example ``octave_code(Matrix([[x**2, sin(pi*x*y), ceiling(x)]]))`` gives the string ``[x.^2 sin(pi*x.*y) ceil(x)]``.

* New general 3D vector package at ``sympy.vector``.  This package provides a 3D vector object with the Del, gradient, divergence, curl, and operators. It supports arbitrary rotations of Cartesian coordinate systems and arbitrary locations of points.

Compatibility breaks
====================

* All usage of inequalities (``>``, ``>=``, ``<``, ``<=``) on SymPy objects will now return SymPy's ``S.true`` or ``S.false`` singletons instead of Python's ``True`` or ``False`` singletons.  Code that checks for, e.g., ``(a < b) is True`` should be changed to ``(a < b) == True`` or ``(a < b) == S.true``.  Use of ``is`` is not recommended here.

* The ``subset()`` method in ``sympy.core.sets`` is marked as being deprecated and will be removed in a future release (:sympyissue:`7460`). Instead, the ``is_subset()`` method should be used.

* Previously, if you compute the series expansion at a point other than 0, the result was shifted to 0.  Now SymPy returns the usual series expansion, see :sympypull:`2427`.

* In ``physics.mechanics``, ``KanesMethod.linearize`` has a new interface. Old code should be changed to use this instead. See docstring for information.

* ``physics.gaussopt`` has been moved to ``physics.optics.gaussopt``. You can still import it from the previous location but it may result in a deprecation warning.

* This is the last release with the bundled `mpmath library <http://mpmath.org/>`_. In the next release you will have to install this library from the official site.

* Previously ``lambdify`` would convert ``Matrix`` to ``numpy.matrix`` by default. This behavior is being deprecated, and will be completely phased out with the release of 0.7.7. To use the new behavior now set the modules kwarg to ``[{'ImmutableMatrix': numpy.array}, 'numpy']``. If lambdify will be used frequently it is recommended to wrap it with a ``partial`` as so: ``lambdify = functools.partial(lambdify, modules=[{'ImmutableMatrix': numpy.array}, 'numpy'])``. For more information see :sympyissue:`7853` and the ``lambdify`` docstring.

* ``Set.complement`` doesn't exists as an attribute anymore. Now we have a method ``Set.complement(<universal_set>)`` which complements the given universal set.

* Removed is_finite assumption (see :sympypull:`7891`).  Use instead a combination of ``is_bounded and is_nonzero`` assumptions.

* is_bounded and is_unbounded assumptions were renamed to is_finite and is_infinite (see :sympypull:`7947`).

* Removed is_infinitesimal assumption (see :sympypull:`7995`).

* Removed is_real property for Sets, use ``Set.is_subset(Reals)`` instead (see :sympypull:`7996`).

* For generic symbol ``x`` (SymPy's symbols are not bounded by default), inequalities with ``oo`` are no longer evaluated as they were before, e.g. ``x < oo`` no longer evaluates to True).  See :sympypull:`7861`.

* ``CodeGen`` has been refactored to make it easier to add other languages.  The main high-level tool is still ``utilities.codegen.codegen``.  But if you previously used the ``Routine`` class directly, note its ``__init__`` behaviour has changed; the new ``utilities.codegen.make_routine`` is recommended instead and by default retains the previous C/Fortran behaviour.  If needed, you can still instantiate ``Routine`` directly; it only does minimal sanity checking on its inputs.  See :sympypull:`8082`.

* ``FiniteSet([1, 2, 3, 4])`` syntax not supported anymore, use ``FiniteSet(1, 2, 3, 4)`` instead.  See :sympypull:`7622`.

Minor changes
=============

* Updated the parsing module to allow sympification of lambda statements to their SymPy equivalent.
* Lambdify can now use ``numexpr`` by specifying ``modules='numexpr'``.
* Use ``with evaluate(False)`` context manager to control automatic evaluation.  E.g. ``with evaluate(False): x + x`` is actually ``x + x``, not ``2*x``.
* IndexedBase and Indexed are changed to be commutative by default.
* ``sympy.core.sets`` moved to ``sympy.sets``.
* Changes in ``sympy.sets``:

  - Infinite ``Range`` is now allowed. See :sympypull:`7741`.
  - ``is_subset()``: The ``is_subset()`` method deprecates the ``subset()`` method.  ``self.is_subset(other)`` checks if ``self`` is a subset of ``other``. This is different from ``self.subset(other)``, which checked if ``other`` is a subset of ``self``.
  - ``is_superset()``: A new method ``is_superset()`` method is now available.  ``self.is_superset(other)`` checks if ``self`` is a superset of ``other``.
  - ``is_proper_subset`` and ``is_proper_superset``: Two new methods allow checking if one set is the proper subset and proper superset of another respectively. For e.g. ``self.is_proper_subset(other)`` and ``self.is_proper_superset(other)`` checks if ``self`` is the proper subset of ``other`` and if ``self`` is the proper superset of ``other`` respectively.
  - ``is_disjoint()``: A new method for checking if two sets are disjoint.
  - ``powerset()``: A new method ``powerset()`` has been added to find the power set of a set.
  - The cardinality of a ``ProductSet`` can be found using the ``len()`` function.

* Changes in ``sympy.plot.plot_implicit``:

  - The ``plot_implicit`` function now also allows explicitly specifying the symbols to plot on the X and Y axes. If not specified, the symbols will be assigned in the order they are sorted.
  - The ``plot_implicit`` function also allows axes labels for the plot to be specified.

* rules for simplification of ImageSet were added :sympypull:`7625`.  As a result ``{x | x ∊ ℤ}`` now simplifies to ``ℤ`` and ``{sin(n) | n ∊ {tan(m) | m ∊ ℤ}}`` automatically simplifies to ``{sin(tan(m)) | m ∊ ℤ}``.
* coth(0) now returns Complex Infinity.  See :sympypull:`7634`.
* dioptre is added to ``physics.units``.  See :sympypull:`7782`.
* ``replace`` now respects commutativity, see :sympypull:`7752`.
* The CCodePrinter gracefully handles Symbols which have string representations that match C reserved words, see :sympypull:`8199`.
* ``limit`` function now returns an unevaluated ``Limit`` instance if it can't compute given limit, see :sympypull:`8213`.
