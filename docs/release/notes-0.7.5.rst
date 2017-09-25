===========
SymPy 0.7.5
===========

22 Feb 2014

Major changes
=============

* The version of mpmath included in SymPy has been updated to 0.18.

* New routines for efficiently compute the *dispersion* of a polynomial
  or a pair thereof.

* Fancy indexing of matrices is now provided, e.g. ``A[:, [1, 2, 5]]`` selects all rows and only 3 columns.

* Enumeration of multiset partitions is now based on an implementation of
  Algorithm 7.1.2.5M from Knuth's *The Art of Computer Programming*.  The new
  version is much faster, and includes fast methods for enumerating only those
  partitions with a restricted range of sizes, and counting multiset
  partitions.  (See the new file ``sympy.utilities.enumerative.py``.)

* ``distance`` methods were added to Line and Ray to compute the shortest distance to them from a point.

* The ``normal_lines`` method was added to Ellipse to compute the lines from a point that strike the Ellipse at a normal angle.

* ``inv_quick`` and ``det_quick`` were added as functions in solvers.py to facilitate fast solution of small symbolic matrices; their use in solve has reduced greatly the time needed to solve such systems.

* ``solve_univariate_inequality`` has been added to sympy.solvers.inequalities.py.

* ``as_set`` attribute for Relationals and Booleans has been added.

* Several classes and functions strictly associated with vector calculus were moved from ``sympy.physics.mechanics`` to a new package ``sympy.physics.vector``. (See :sympypull:`2732`, :sympypull:`2862` and :sympypull:`2894`).

* New implementation of the Airy functions ``Ai`` and ``Bi`` and their derivatives
  ``Ai'`` and ``Bi'`` (called ``airyai``, ``airybi``, ``airyaiprime`` and ``airybiprime``,
  respectively). Most of the usual features of SymPy special function are
  present.  Notable exceptions are Gruntz limit computation helpers and
  meijerg special functions integration code.

* Euler-Lagrange equations (function ``euler_equations``) in a new package ``sympy.calculus`` (:sympypull:`2431`).

Compatibility breaks
====================

* the ``submatrix`` method of matrices was removed; access the functionality by
  providing slices or list of rows/columns to matrix directly,
  e.g. ``A[:, [1, 2]]``.

* ``Matrix([])`` and ``Matrix([[]])`` now both return a 0x0 matrix

* ``terms_gcd`` no longer removes a -1.0 from expressions

* ``extract_multiplicatively`` will not remove a negative Number from a positive one, so
  ``(4*x*y).extract_multiplicatively(-2*x)`` will return None.

* the shape of the result from ``M.cross(B)`` now has the same shape as matrix M.

* The factorial of negative numbers is now ``zoo`` instead of ``0``. This is
  consistent with the definition ``factorial(n) = gamma(n + 1)``.

* ``1/0`` returns ``zoo``, not ``oo`` (:sympypull:`2813`).

* ``zoo.is_number is True`` (:sympypull:`2823`).

* ``oo < I`` raises ``TypeError``, just as for finite numbers (:sympypull:`2734`).

* ``1**oo == nan`` instead of ``1``, better documentation for ``Pow`` class (:sympypull:`2606`).

Minor changes
=============

* Some improvements to the gamma function.

* ``generate_bell`` now generates correct permutations for any number of elements.

* It is no longer necessary to provide ``nargs`` to objects subclassed from
  Function unless an ``eval`` class method is not defined. (If ``eval`` is defined,
  the number of arguments will be inferred from its signature.)

* geometric ``Point`` creation will be faster since simplification is done only on Floats

* Some improvements to the intersection method of the Ellipse.

* solutions from solve of equations involving multiple log terms are more robust

* ``idiff`` can now return higher order derivatives

* Added ``to_matrix()`` method to ``sympy.physics.vector.Vector`` and ``sympy.physics.dyadic.Dyadic``. (:sympypull:`2686`).

* Printing improvements for ``sympy.physics.vector`` objects and mechanics printing. (See :sympypull:`2687`, :sympypull:`2728`, :sympypull:`2772`, :sympypull:`2862` and :sympypull:`2894`).

* Functions with LaTeX symbols now print correct LaTeX. (:sympypull:`2772`).

* ``init_printing`` has several new options, including a flag ``print_builtin`` to prevent SymPy
  printing of basic Python types (:sympypull:`2683`), and flags to let you supply custom
  printers (:sympypull:`2894`).

* improvements in evaluation of ``imageset`` for Intervals (:sympypull:`2723`).

* Set properties to determine boundary and interior (:sympypull:`2744`).

* input to a function created by lambdify no longer needs to be flattened.
