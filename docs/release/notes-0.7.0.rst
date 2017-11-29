===========
SymPy 0.7.0
===========

28 Jun 2011

Major Changes
=============

* Polys

  - New internal representations of dense and sparse polynomials (see :commit:`6aecdb7`, :commit:`31c9aa4`)
  - Implemented algorithms for real and complex root isolation and counting (see :commit:`3acac67`, :commit:`4b75dae`, :commit:`fa1206e`, :commit:`103b928`, :commit:`45c9b22`, :commit:`8870c8b`, :commit:`b348b30`)
  - Improved Gröbner bases algorithm (see :commit:`ff65e9f`, :commit:`891e4de`, :commit:`310a585`)
  - Field isomorphism algorithm (see :commit:`b097b01`, :commit:`08482bf`)
  - Implemented efficient orthogonal polynomials (see :commit:`b8fbd59`)
  - Added configuration framework for polys (see :commit:`33d8cdb`, :commit:`7eb81c9`)
  - Function for computing minimal polynomials (see :commit:`88bf187`, :commit:`f800f95`)
  - Function for generating Viete's formulas (see :commit:`1027408`)
  - ``roots()`` supports more classes of polynomials (e.g. cyclotomic) (see :commit:`d8c8768`, :commit:`75c8d2d`)
  - Added a function for recognizing cyclotomic polynomials (see :commit:`b9c2a9a`)
  - Added a function for computing Horner form of polynomials (see :commit:`8d235c7`)
  - Added a function for computing symmetric reductions of polynomials (see :commit:`6d560f3`)
  - Added generators of Swinnerton-Dyer, cyclotomic, symmetric, random and interpolating polynomials (see :commit:`dad03dd`, :commit:`6ccf20c`, :commit:`dc728d6`, :commit:`2f17684`, :commit:`3004db8`)
  - Added a function computing isolation intervals of algebraic numbers (see :commit:`37a58f1`)
  - Polynomial division (``div()``, ``rem()``, ``quo()``) now defaults to a field (see :commit:`a72d188`)
  - Added wrappers for numerical root finding algorithms (see :commit:`0d98945`, :commit:`f638fcf`)
  - Added symbolic capabilities to ``factor()``, ``sqf()`` and related functions (see :commit:`d521c7f`, :commit:`548120b`, :commit:`f6f74e6`, :commit:`b1c49cd`, :commit:`3527b64`)
  - ``together()`` was significantly improved (see :commit:`dc327fe`)
  - Added support for iterable containers to ``gcd()`` and ``lcm()`` (see :commit:`e920870`)
  - Added a function for constructing domains from coefficient containers (see :commit:`a8f20e6`)
  - Implemented greatest factorial factorization (see :commit:`d4dbbb5`)
  - Added partial fraction decomposition algorithm based on undetermined coefficient approach (see :commit:`9769d49`, :commit:`496f08f`)
  - ``RootOf`` and ``RootSum`` were significantly improved (see :commit:`f3e432`, :commit:`4c88be6`, :commit:`41502d7`)
  - Added support for gmpy (GNU Multiple Precision Arithmetic Library) (see :commit:`38e1683`)
  - Allow to compile ``sympy.polys`` with Cython (see :commit:`afb3886`)
  - Improved configuration of variables in ``Poly`` (see :commit:`22c4061`)
  - Added documentation based on Wester's examples (see :commit:`1c23792`)
  - Irreducibility testing over finite fields (see :commit:`17e8f1f`)
  - Allow symmetric and non-symmetric representations over finite fields (see :commit:`60fbff4`)
  - More consistent factorization forms from ``factor()`` and ``sqf()`` (see :commit:`5df77f5`)
  - Added support for automatic recognition algebraic extensions (see :commit:`7de602c`)
  - Implemented Collins' modular algorithm for computing resultants (see :commit:`950969b`)
  - Implemented Berlekamp's algorithm for factorization over finite fields (see :commit:`70353e9`)
  - Implemented Trager's algorithm for factorization over algebraic number fields (see :commit:`bd0be06`)
  - Improved Wang's algorithm for efficient factorization of multivariate polynomials (see :commit:`425e225`)

* Quantum

  - Symbolic, abstract dirac notation in ``sympy.physics.quantum``. This includes operators,
    states (bras and kets), commutators, anticommutators, dagger, inner products, outer
    products, tensor products and Hilbert spaces
  - Symbolic quantum computing framework that is based on the general capabilities in
    ``sympy.physics.quantum``. This includes qubits (``sympy.physics.quantum.qubit``), gates
    (``sympy.physics.quantum.gate``), Grover's algorithm (``sympy.physics.quantum.grover``),
    the quantum Fourier transform (``sympy.physics.quantum.qft``), Shor's algorithm
    (``sympy.physics.quantum.shor``) and circuit plotting (``sympy.physics.quantum.circuitplot``)
  - Second quantization framework that inclues creation/anihilation operators for
    both Fermions and Bosons and Wick's theorem for Fermions (``sympy.physics.secondquant``).
  - Symbolic quantum angular momentum (spin) algebra (``sympy.physics.quantum.spin``)
  - Hydrogen wave functions (Schroedinger) and energies (both Schroedinger and Dirac)
  - Wave functions and energies for 1D harmonic oscillator
  - Wave functions and energies for 3D spherically symmetric harmonic oscillator
  - Wigner and Clebsch Gordan coefficients

* Everything else

  - Implement symarray, providing numpy nd-arrays of symbols.
  - update mpmath to 0.16
  - Add a tensor module (see `this report <https://code.google.com/archive/p/sympy/wikis/CodeGenerationReport.wiki>`_)
  - A lot of stuff was being imported with ``from sympy import *`` that shouldn't have been (like ``sys``).  This has been fixed.

* Assumptions:

  - Refine
  - Added predicates (see :commit:`7c0b857`, :commit:`53f0e1a`, :commit:`d1dd6a3`)
  - Added query handlers for algebraic numbers (see :commit:`f3bee7a`)
  - Implement a SAT solver (see `this <https://code.google.com/archive/p/sympy/wikis/SuperchargingAssumptionsReport.wiki>`_, :commit:`2d96329`, :commit:`acfbe75`, etc)

* Concrete

  - Finalized implementation of Gosper's algorithm (see :commit:`0f187e5`, :commit:`5888024`)
  - Removed redundant ``Sum2`` and related classes (see :commit:`ef1f6a7`)

* Core:

  - Split ``Atom`` into ``Atom`` and ``AtomicExpr`` (see :commit:`965aa91`)
  - Various ``sympify()`` improvements
  - Added functionality for action verbs (many functions can be called both as global functions and as methods e.g. ``a.simplify() == simplify(a)``)
  - Improve handling of rational strings (see :commit:`053a045`, :sympyissue:`4877`)
  - Major changes to factoring of integers (see :commit:`273f450`, :sympyissue:`5102`)
  - Optimized ``.has()`` (see :commit:`c83c9b0`, :sympyissue:`5079`, :commit:`d86d08f`)
  - Improvements to power (see :commit:`c8661ef`, :sympyissue:`5062`)
  - Added range and lexicographic syntax to ``symbols()`` and ``var()`` (see :commit:`f6452a8`, :commit:`9aeb220`, :commit:`957745a`)
  - Added ``modulus`` argument to ``expand()`` (see :commit:`1ea5be8`)
  - Allow to convert ``Interval`` to relational form (see :commit:`4c269fe`)
  - SymPy won't manipulate minus sign of expressions any more (see :commit:`6a26941`, :commit:`9c6bf0f`, :commit:`e9f4a0a`)
  - ``Real`` and ``.is_Real`` were renamed to ``Float`` and ``.is_Float``.  ``Real`` and ``.is_Real`` still remain as deprecated shortcuts to ``Float`` and ``is_Float`` for backwards compatibility. (see :commit:`abe1c49`)
  - Methods coeff and as_coefficient are now non-commutative aware. (see :commit:`a4ea170`)

* Geometry:

  - Various improvements to Ellipse
  - Updated documentation to numpy standard
  - Polygon and Line improvements
  - Allow all geometry objects to accept a tuple as ``Point`` args

* Integrals:

  - Various improvements (see e.g. :sympyissue:`4871`, :sympyissue:`5098`, :sympyissue:`5091`, :sympyissue:`5086`)

* isympy

  - Fixed the ``-p`` switch (see :commit:`e8cb04a`)
  - Caching can be disabled using ``-C`` switch (see :commit:`0d8d748`)
  - Ground types can be set using ``-t`` switch (see :commit:`75734f8`)
  - Printing ordering can be set using ``-o`` switch (see :commit:`fcc6b13`, :commit:`4ec9dc5`)

* Logic

  - implies object adheres to negative normal form
  - Create new boolean class, ``logic.boolalg.Boolean``
  - Added XOR operator (^) support
  - Added If-then-else (ITE) support
  - Added the dpll algorithm

* Functions:

  - Added Piecewise, B-splines
  - Spherical Bessel function of the second kind implemented
  - Add series expansions of multivariate functions (see :commit:`d4d351d`)

* Matrices:

  - Add elementwise product (Hadamard product)
  - Extended QR factorization for general full ranked mxn matrices
  - Remove deprecated functions ``zero()``, ``zeronm()``, ``one()`` (see :commit:`5da0884`)
  - Added cholesky and LDL factorizations, and respective solves.
  - Added functions for efficient triangular and diagonal solves.
  - ``SMatrix`` was renamed to ``SparseMatrix`` (see :commit:`acd1685`)

* Printing:

  - Implemented pretty printing of binomials (see :commit:`58c1dad`)
  - Implemented pretty printing of Sum() (see :commit:`84f2c22`, :commit:`95b4321`)
  - ``sympy.printing`` now supports ordering of terms and factors (see :commit:`859bb33`)
  - Lexicographic order is now the default. Now finally things will print as ``x**2 + x + 1`` instead of ``1 + x + x**2``, however series still print using reversed ordering, e.g. ``x - x**3/6 + O(x**5)``. You can get the old order (and other orderings) by setting the ``-o`` option to isympy (see :commit:`08b4932`, :commit:`a30c5a3`)

* Series:

  - Implement a function to calculate residues, ``residue()``
  - Implement nseries and lseries to handle ``x0 != 0``, series should be more robust now (see :commit:`2c99999`, :sympyissue:`5221` - :sympyissue:`5223`)
  - Improvements to Gruntz algorithm

* Simplify:

  - Added ``use()`` (see :commit:`147c142`)
  - ``ratsimp()`` now uses ``cancel()`` and ``reduced()`` (see :commit:`108fb41`)
  - Implemented EPath (see :commit:`696139d`, :commit:`bf90689`)
  - a new keyword ``rational`` was added to nsimplify which will replace Floats with Rational approximations. (see :commit:`053a045`)

* Solvers:

  - ODE improvements (see :commit:`d12a2aa`, :commit:`3542041`; :commit:`73fb9ac`)
  - Added support for solving inequalities (see :commit:`328eaba`, :commit:`8455147`, :commit:`f8fcaa7`)

* Utilities:

  - Improve cartes, for generating the Cartesian product (see :commit:`b1b10ed`)
  - Added a function computing topological sort of graphs (see :commit:`b2ce27b`)
  - Allow to setup a customized printer in ``lambdify()`` (see :commit:`c1ad905`)
  - ``flatten()`` was significantly improved (see :commit:`31ed8d7`)
  - Major improvements to the Fortran code generator (see `report <https://code.google.com/archive/p/sympy/wikis/CodeGenerationReport.wiki>`_, :commit:`3383aa3`, :commit:`7ab2da2`, etc)

Compatibility breaks
====================

* This will be the last release of SymPy to support Python 2.4.  Dropping support for Python 2.4 will let us move forward with things like supporting Python 3, and will let us use things that were introduced in Python 2.5, like with-statement context managers.
* no longer support creating matrices without brackets (see :sympyissue:`4029`)
* Renamed ``sum()`` to ``summation()`` (see :commit:`3e763a8`, :sympyissue:`4475`, :sympyissue:`4826`). This was changed so that it no longer overrides the built-in ``sum()``. The unevaluated summation is still called ``Sum()``.
* Renamed ``abs()`` to ``Abs()`` (see :commit:`64a12a4`, :sympyissue:`4826`).  This was also changed so that it no longer overrides the built-in ``abs()``. Note that because of ``__abs__`` magic, you can still do ``abs(expr)`` with the built-in ``abs()``, and it will return ``Abs(expr)``.
* Renamed ``max_()`` and ``min_()`` to now ``Max()`` and ``Min()`` (see :commit:`99a271e`, :sympyissue:`5252`)
* Changed behaviour of ``symbols()``. ``symbols('xyz')`` gives now a single symbol (``'xyz'``), not three (``'x'``, ``'y'`` and ``'z'``) (see :commit:`f6452a8`). Use ``symbols('x,y,z')`` or ``symbols('x y z')`` to get three symbols. The ``each_char`` option will still work but is being deprecated.
* Split class ``Basic`` into new classes ``Expr``, ``Boolean`` (see :commit:`a0ab479`, :commit:`635d89c`). Classes that are designed to be part of standard symbolic expressions (like ``x**2*sin(x)``) should subclass from ``Expr``. More generic objects that do not work in symbolic expressions but still want the basic SymPy structure like ``.args`` and basic methods like ``.subs()`` should only subclass from ``Basic``.
* ``as_basic()`` method was renamed to ``as_expr()`` to reflect changes in the core (see :commit:`e61819d`, :commit:`80dfe91`)
* Methods ``as_coeff_terms`` and ``as_coeff_factors`` were renamed to ``as_coeff_mul`` and ``as_coeff_add``, respectively.
* Removed the ``trim()`` function.  The function is redundant with the new polys.  Use the ``cancel()`` function instead.
* The ``assume_pos_real`` option to ``logcombine()`` was renamed to ``force`` to be consistant with similar ``force`` options to other functions.

In addition to the more noticeable changes listed above, there have been numerous other smaller additions, improvements and bug fixes in the ~2000 commits in this release. See the git log for a full list of all changes.  The command ``git log sympy-0.6.7..sympy-0.7.0`` will show all commits made between this release and the last. You can also see the issues closed since the last release `here <https://github.com/sympy/sympy/issues?utf8=%E2%9C%93&q=is%3Aissue%20closed%3A%222010-03-17%20..%202011-06-28%22>`_.
