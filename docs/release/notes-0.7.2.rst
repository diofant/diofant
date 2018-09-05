===========
SymPy 0.7.2
===========

16 Oct 2012

Major Changes
=============

* Python 3 support

  - SymPy now supports Python 3. The officially supported versions are 3.2 and
    3.3, but 3.1 should also work in a pinch. The Python 3-compatible tarballs
    will be provided separately, but it is also possible to download Python 2 code
    and convert it manually, via the bin/use2to3 utility. See the README for more.

* PyPy support

  - All SymPy tests pass in recent nightlies of PyPy, and so it should have full
    support as of the next version after 1.9.

* Combinatorics

  - A new module called Combinatorics was added which is the result of a
    successful GSoC project. It attempts to replicate the functionality of
    Combinatorica and currently has full featured support for Permutations,
    Subsets, Gray codes and Prufer codes.

  - In another GSoC project, facilities from computational group theory were added
    to the combinatorics module, mainly following the book "Handbook of
    computational group theory". Currently only permutation groups are
    supported. The main functionalities are: basic properties (orbits,
    stabilizers, random elements...), the Schreier-Sims algorithm (three
    implementations, in increasing speed: with Jerrum's filter, incremental, and
    randomized (Monte Carlo)), backtrack searching for subgroups with certain
    properties.

* Definite Integration

  - A new module called meijerint was added, which is also the result of a
    successful GSoC project. It implements a heuristic algorithm for (mainly)
    definite integration, similar to the one used in Mathematica. The code is
    automatically called by the standard integrate() function. This new algorithm
    allows computation of important integral transforms in many interesting cases,
    so helper functions for Laplace, Fourier and Mellin transforms were added as
    well.

* Random Variables

  - A new module called stats was added. This introduces a RandomSymbol type which
    can be used to model uncertainty in expressions.

* Matrix Expressions

  - A new matrix submodule named expressions was added. This introduces a
    MatrixSymbol type which can be used to describe a matrix without explicitly
    stating its entries. A new family of expression types were also added:
    Transpose, Inverse, Trace, and BlockMatrix. ImmutableMatrix was added so that
    explicitly defined matrices could interact with other SymPy expressions.

* Sets

  - A number of new sets were added including atomic sets like FiniteSet, Reals,
    Naturals, Integers, UniversalSet as well as compound sets like ProductSet and
    TransformationSet. Using these building blocks it is possible to build up a
    great variety of interesting sets.

* Classical Mechanics

  - A physics submodule named machanics was added which assists in formation of
    equations of motion for constrained multi-body systems. It is the result of 3
    GSoC projects. Some nontrivial systems can be solved, and examples are
    provided.

* Quantum Mechanics

  - Density operator module has been added. The operator can be initialized with
    generic Kets or Qubits. The Density operator can also work with TensorProducts
    as arguments. Global methods are also added that compute entropy and fidelity
    of states. Trace and partial-trace operations can also be performed on these
    density operators.

  - To enable partial trace operations a Tr module has been added to the core
    library.  While the functionality should remain same, this module is likely to
    be relocated to an alternate folder in the future. One can currently also use
    sympy.core.Tr to work on general trace operations, but this module is what is
    needed to work on trace and partial-trace operations on any
    sympy.physics.quantum objects.

  - The Density operators, Tr and Partial trace functionality was implemented as
    part of student participation in GSoC 2012.

  - Expanded angular momentum to include coupled-basis states and product-basis
    states. Operators can also be treated as acting on the coupled basis (default
    behavior) or on one component of the tensor product states. The methods for
    coupling and uncoupling these states can work on an arbitrary number of
    states.  Representing, rewriting and applying states and operators between
    bases has been improved.

* Commutative Algebra

  - A new module ``agca`` was started which seeks to support computations in
    commutative algebra (and eventually algebraic geometry) in the style of
    Macaulay2 and Singular. Currently there is support for computing Gröbner
    bases of modules over a (generalized) polynomial ring over a field. Based on
    this, there are algorithms for various standard problems in commutative
    algebra, e.g., computing intersections of submodules, equality tests in
    quotient rings, etc...

* Plotting Module

  - A new plotting module has been added which uses Matplotlib as its
    back-end. The plotting module has functions to plot the following:

    * 2D line plots
    * 2D parametric plots.
    * 2D implicit and region plots.
    * 3D surface plots.
    * 3D parametric surface plots.
    * 3D parametric line plots.

* Differential Geometry

  - Thanks to a GSoC project the beginning of a new module covering the theory of
    differential geometry was started. It can be imported with
    ``sympy.diffgeom``. It is based on "Functional Differential Geometry" by Sussman
    and Wisdom. Currently implemented are scalar, vector and form fields over
    manifolds as well as covariant and other derivatives.

Compatibility breaks
====================

- The KroneckerDelta class was moved from ``sympy/physics/quantum/kronecker.py`` to
  ``sympy/functions/special/tensor_functions.py``.

- Merged the KroneckerDelta class in ``sympy/physics/secondquant.py`` with the
  class above.

- The Dij class in ``sympy/functions/special/tensor_functions.py`` was replaced
  with KroneckerDelta.

- The errors raised for invalid ``float`` calls on SymPy objects were changed in
  order to emulate more closely the errors raised by the standard library. The
  ``__float__`` and ``__complex__`` methods of ``Expr`` are concerned with that
  change.

- The ``solve()`` function returns empty lists instead of ``None`` objects if no
  solutions were found. Idiomatic code of the form ``sol = solve(...); if
  sol:...`` will not be affected by this change.

- Piecewise no longer accepts a Set or Interval as a condition. One should
  explicitly specify a variable using ``Set().contains(x)`` to obtain a valid
  conditional.

- The statistics module has been deprecated in favor of the new stats module.

- ``sympy/galgebra/GA.py``:

  * ``set_main()`` is no longer needed
  * ``make_symbols()`` is deprecated (use ``sympy.symbols()`` instead)
  * the symbols used in this package are no longer broadcast to the main program

- The classes for Infinity, NegativeInfinity, and NaN no longer subclass from
  Rational.  Creating a Rational with 0 in the denominator will still return
  one of these classes, however.

Minor changes
=============

- A new module ``gaussopt`` was added supporting the most basic constructions
  from Gaussian optics (ray tracing matrices, geometric rays and Gaussian
  beams).

- New classes were added to represent the following special functions:
  classical and generalized exponential integrals (Ei, expint), trigonometric
  (Si, Ci) and hyperbolic integrals (Shi, Chi), the polylogarithm (polylog)
  and the Lerch transcendent (lerchphi). In addition to providing all the
  standard sympy functionality (differentiation, numerical evaluation,
  rewriting ...), they are supported by both the new meijerint module and the
  existing hypergeometric function simplification module.

- An ImmutableMatrix class was created. It has the same interface and
  functionality of the old Matrix but is immutable and inherits from Basic.

- A new function in ``geometry.util`` named ``centroid`` was added which will
  calculate the centroid of a collection of geometric entities. And the
  polygon module now allows triangles to be instantiated from combinations of
  side lengths and angles (using keywords sss, asa, sas) and defines utility
  functions to convert between degrees and radians.

- In ``ntheory.modular`` there is a function (``solve_congruence``) to solve
  congruences such as "What number is 2 mod 3, 3 mod 5 and 2 mod 7?"

- A utility function named ``find_unit`` has been added to physcis.units that
  allows one to find units that match a given pattern or contain a given unit.

- There have been some additions and modifications to Expr's methods:

  - Although the problem of proving that two expressions are equal is in general
    a difficult one (since whatever algorithm is used, there will always be an
    expression that will slip through the algorithm) the new method of Expr
    named ``equals`` will do its best to answer whether A equals B: A.equals(B)
    might given True, False or None.

  - coeff now supports a third argument ``n`` (which comes 2nd now, instead of
    ``right``). This ``n`` is used to indicate the exponent on x which one seeks:
    ``(x**2 + 3*x + 4).coeff(x, 1)`` -> ``3``.  This makes it possible to extract the
    constant term from a polynomial: ``(x**2 + 3*x + 4).coeff(x, 0)`` -> ``4``.

  - The method ``round`` has been added to round a SymPy expression to a given a
    number of decimal places (to the left or right of the decimal point).

- divmod is now supported for all SymPy numbers.

- In the simplify module, the algorithms for denesting of radicals
  (sqrtdenest) and simplifying gamma functions (in combsimp) has been
  significantly improved.

- The mathematica-similar ``TableForm`` function has been added to the
  printing.tableform module so one can easily generate tables with headings.

- The expand API has been updated.  ``expand()`` now officially supports
  arbitrary ``_eval_expand_hint()`` methods on custom
  objects. ``_eval_expand_hint()`` methods are now only responsible for
  expanding the top-level expression.  All ``deep=True`` related logic happens
  in ``expand()`` itself.  See the docstring of ``expand()``
  for more information and an example.

- Two options were added to ``isympy`` to aid in interactive usage.  ``isympy -a``
  automatically creates symbols, so that typing something like ``a`` will give
  ``Symbol('a')``, even if you never typed ``a = Symbol('a')`` or ``var('a')``.
  ``isympy -i`` automatically wraps integer literals with Integer, so that ``1/2``
  will give ``Rational(1, 2)`` instead of ``0.5``.  ``isympy -I`` is the same as
  ``isympy -a -i``.  ``isympy -I`` makes isympy act much more like a traditional
  interactive computer algebra system.  These both require IPython.

- The official documentation at https://docs.sympy.org/ now includes an
  extension that automatically hooks the documentation examples in to
  `SymPy Live <https://live.sympy.org>`_.

In addition to the more noticeable changes listed above, there have been
numerous smaller additions, improvements and bug fixes in the commits in
this release. See the git log for a full list of all changes. The command
``git log sympy-0.7.1..sympy-0.7.2`` will show all commits made between this
release and the last. You can also see the issues closed since the last
release `here <https://github.com/sympy/sympy/issues?utf8=%E2%9C%93&q=is%3Aissue%20closed%3A%222011-07-29%20..%202012-10-16%22>`_.
