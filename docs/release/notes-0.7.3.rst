===========
SymPy 0.7.3
===========

13 Jul 2013

Major changes
=============

* Integration

  - This release includes Risch integration algorithm from
    `Aaron Meurer's 2010 Google Summer of Code project <https://github.com/sympy/sympy/wiki/GSoC-2010-Report-Aaron-Meurer:-Risch-Integration>`_.
    This makes ``integrate`` much more powerful and much faster for the supported
    functions.  The algorithm is called automatically from ``integrate()``.  For
    now, only transcendental elementary functions containing ``exp`` or ``log`` are
    supported.  To access the algorithm directly, use ``integrate(expr, x,
    risch=True)``.  The algorithm has the ability to prove that integrals are
    nonelementary.  To determine if a function is nonelementary, integrate using
    ``risch=True``.  If the resulting ``Integral`` class is an instance of
    ``NonElementaryIntegral``, then it is not elementary (otherwise, that part of
    the algorithm has just not been implemented yet).

* ODE

  - Built basic infrastructure of the PDE module (:sympypull:`1970`)

* Theano Interaction

  - SymPy expressions can now be translated into `Theano <http://deeplearning.net/software/theano/>`_ expressions for numeric evaluation.  This includes most standard scalar operations (e.g. ``sin``, ``exp``, ``gamma``, but not ``beta`` or ``MeijerG``) and matrices.  This system generally outperforms ``lambdify`` and ``autowrap`` but does require Theano to be installed.

* Matrix Expressions

  - Matrix expressions now support inference using the new assumptions system.  New predicates include ``invertible``, ``symmetric``, ``positive_definite``, ``orthogonal``, ....

  - New operators include ``Adjoint``, ``HadamardProduct``, ``Determinant``, ``MatrixSlice``, ``DFT``.  Also, preliminary support exists for factorizations like ``SVD`` and ``LU``.

* Context manager for New Assumptions

  - Added the ``with assuming(*facts)`` context manager for new assumptions.  See `blogpost <https://web.archive.org/web/20181111092915/https://matthewrocklin.com/blog/work/2013/02/05/Assuming>`_.

Compatibility breaks
====================

- This is the last version of SymPy to support Python 2.5.

- The IPython extension, i.e., ``%load_ext sympy.interactive.ipythonprinting``
  is deprecated.  Use ``from sympy import init_printing; init_printing()``
  instead. See :sympyissue:`7013`.

- The ``viewer='file'`` option to ``preview`` without a file name is
  deprecated. Use ``filename='name'`` in addition to ``viewer='file'``. See
  :sympyissue:`7018`.

- The deprecated syntax ``Symbol('x', dummy=True)``, which had been deprecated
  since 0.7.0, has been removed. Use ``Dummy('x')`` or ``symbols('x', cls=Dummy)``
  instead. See :sympyissue:`6477`.

- The deprecated ``Expr`` methods ``as_coeff_terms`` and ``as_coeff_factors``, which
  have been deprecated in favor of ``as_coeff_mul`` and ``as_coeff_add``,
  respectively (see also ``as_coeff_Mul`` and ``as_coeff_Add``), were removed.
  The methods had been deprecated since SymPy 0.7.0.  See
  :sympyissue:`6476`.

- The spherical harmonics have been completely rewritten. See :sympypull:`1510`.

Minor changes
=============

* Solvers

  - Added enhancements and improved the methods of solving exact differential equation.  See :sympypull:`1955` and :sympypull:`1823`.
  - Support for differential equations with linear coefficients and those that can be reduced to separable and linear form.  See :sympypull:`1940`, :sympypull:`1864` and :sympypull:`1883`.
  - Support for first order linear general PDE's with constant coefficients (:sympypull:`2109`).
  - Return all found independent solutions for underdetermined systems.
  - Handle recursive problems for which ``y(0) = 0``.
  - Handle matrix equations.

* Integration

  - ``integrate`` will split out integrals into Piecewise expressions when
    conditions must hold for the answer to be true. For example,
    ``integrate(x**n, x)`` now gives ``Piecewise((log(x), Eq(n, -1), (x**(n +
    1)/(n + 1), True))`` (previously it just gave ``x**(n + 1)/(n + 1)``).
  - Calculate Gauss-Legendre and Gauss-Laguerre points and weights (:sympypull:`1497`).
  - Various new error and inverse error functions (:sympypull:`1703`).
  - Use in heurisch for more symmetric and nicer results.
  - Gruntz for expintegrals and all new erf*.
  - Li, li logarithmic integrals (:sympypull:`1708`).
  - Integration of li/Li by heurisch (:sympypull:`1712`).
  - elliptic integrals, complete and incomplete.
  - Integration of complete elliptic integrals by meijerg.
  - Integration of Piecewise with symbolic conditions.
  - Fixed many wrong results of DiracDelta integrals.

* Logic

  - Addition of SOPform and POSform functions to sympy.logic to generate boolean expressions from truth tables.
  - Addition of simplify_logic function and enabling ``simplify()`` to reduce logic expressions to their simplest forms.
  - Addition of bool_equals function to check equality of boolean expressions and return a mapping of variables from one expr to other that leads to the equality.
  - Addition of disjunctive normal form methods - to_dnf, is_dnf

* Others

  - gmpy version 2 is now supported
  - Added ``is_algebraic_expr()`` method (:sympypull:`2176`).
  - Many improvements to the handling of noncommutative symbols:

    - Better support in simplification functions, e.g. ``factor``, ``trigsimp``
    - Better integration with ``Order()``
    - Better pattern matching

  - Improved pattern matching including matching the identity.
  - normalizes Jacobi polynomials
  - Quadrature rules for orthogonal polynomials in arbitrary precision
    (hermite, laguerre, legendre, gen_legendre, jacobi)
  - summation of harmonic numbers
  - Many improvements of the polygamma functions
  - evaluation at special arguments
  - Connections to harmonic numbers
  - structured full partial fraction decomposition (mainly interesting for developers)
  - besselsimp improvements
  - Karr summation convention
  - New spherical harmonics
  - improved minimal_polynomial using composition of algebraic numbers (:sympypull:`2038`).
  - faster integer polynomial factorization (:sympypull:`2148`).
  - Euler-Descartes method for quartic equations (:sympypull:`1947`)
  - algebraic operations on tensors (:sympypull:`1700`).
  - tensor canonicalization (:sympypull:`1644`).
  - Handle the simplification of summations and products over a KroneckerDelta.
  - Implemented LaTeX printing of DiracDelta, Heaviside, KroneckerDelta and LeviCivita, also many Matrix expressions.
  - Improved LaTeX printing of fractions, Mul in general.
  - IPython integration and printing issues have been ironed out.
  - Stats now supports discrete distributions (e.g. ``Poisson``) by relying on ``Summation`` objects
  - Added DOT printing for visualization of expression trees
  - Added information about solvability and nilpotency of named groups.
