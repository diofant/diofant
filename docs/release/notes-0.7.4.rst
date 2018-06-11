===========
SymPy 0.7.4
===========

9 Dec 2013

Major changes
=============

* Python 3

  - SymPy now uses a single code-base for Python 2 and Python 3.

* Geometric Algebra

  - The internal representation of a multivector has been changes to more fully
    use the inherent capabilities of SymPy. A multivector is now represented by a
    linear combination of real commutative SymPy expressions and a collection of
    non-commutative SymPy symbols.  Each non-commutative symbol represents a base
    in the geometric algebra of an N-dimensional vector space. The total number of
    non-commutative bases is ``2**N - 1`` (``N`` of which are a basis for the vector
    space) which when including scalars give a dimension for the geometric algebra
    of ``2**N``.  The different products of geometric algebra are implemented as
    functions that take pairs of bases symbols and return a multivector for each
    pair of bases.

  - The LaTeX printing module for multivectors has been rewritten to simply extend
    the existing sympy LaTeX printing module and the sympy LaTeX module is now
    used to print the bases coefficients in the multivector representation instead
    of writing an entire LaTeX printing module from scratch.

  - The main change in the geometric algebra module from the viewpoint of the user
    is the interface for the gradient operator and the implementation of vector
    manifolds:

    - The gradient operator is now implemented as a special vector (the user can
      name it ``grad`` if they wish) so the if ``F`` is a multivector field all the
      operations of ``grad`` on ``F`` can be written ``grad*F``, ``F*grad``, ``grad^F``,
      ``F^grad``, ``grad|F``, ``F|grad``, ``grad<F``, ``F<grad``, ``grad>F``, and ``F>grad`` where
      ``**``, ``^``, ``|``, ``<``, and ``>`` are the geometric product, outer product, inner
      product, left contraction, and right contraction, respectively.

    - The vector manifold is defined as a parametric vector field in an embedding
      vector space. For example a surface in a 3-dimensional space would be a vector
      field as a function of two parameters. Then multivector fields can be defined
      on the manifold. The operations available to be performed on these fields are
      directional derivative, gradient, and projection. The weak point of the
      current manifold representation is that all fields on the manifold are
      represented in terms of the bases of the embedding vector space.

* Classical Cryptography, implements:

  - Affine ciphers
  - Vigenere ciphers
  - Bifid ciphers
  - Hill ciphers
  - RSA and "kid RSA"
  - linear feedback shift registers.

* Common Subexpression Elimination (CSE).  Major changes have been done in
  cse internals resulting in a big speedup for larger expressions.  Some
  changes reflect on the user side:

  - Adds and Muls are now recursively matched (``[w*x, w*x*y, w*x*y*z]`` ǹow turns
    into ``[(x0, w*x), (x1, x0*y)], [x0, x1, x1*z]``)
  - CSE is now not performed on the non-commutative parts of multiplications
    (it avoids some bugs).
  - Pre and post optimizations are not performed by default anymore. The
    ``optimizations`` parameter still exists and ``optimizations='basic'`` can be used
    to apply previous default optimizations. These optimizations could really slow
    down cse on larger expressions and are no guarantee of better results.
  - An ``order`` parameter has been introduced to control whether Adds and Muls
    terms are ordered independently of hashing implementation. The default
    ``order='canonical'`` will independently order the terms. ``order='none'`` will
    not do any ordering (hashes order is used) and will represent a major
    performance improvement for really huge expressions.
  - In general, the output of cse will be slightly different from the previous
    implementation.

* Diophantine Equation Module.
  This is a new addition to SymPy as a result of a GSoC project. With the current
  release, following five types of equations are supported.

  - Linear Diophantine equation, `a_{1}x_{1} + a_{2}x_{2} + . . . + a_{n}x_{n} = b`
  - General binary quadratic equation, `ax^2 + bxy + cy^2 + dx + ey + f = 0`
  - Homogeneous ternary quadratic equation, `ax^2 + by^2 + cz^2 + dxy + eyz + fzx = 0`
  - Extended Pythagorean equation, `a_{1}x_{1}^2 + a_{2}x_{2}^2 + . . . + a_{n}x_{n}^2 = a_{n+1}x_{n+1}^2`
  - General sum of squares, `x_{1}^2 + x_{2}^2 + . . . + x_{n}^2 = k`

* Unification of Sum, Product, and Integral classes

  - A new superclass has been introduced to unify the treatments of indexed
    expressions, such as Sum, Product, and Integral.  This enforced common
    behavior accross the objects, and provides more robust support for a number of
    operations.  For example, Sums and Integrals can now be factored or expanded.
    ``S.subs()`` can be used to substitute for expressions inside a
    Sum/Integral/Product that are independent of the index variables, including
    unknown functions, for instance, ``Integral(f(x), (x, 1, 3)).subs(f(x), x**2)``,
    while ``Sum.change_index()`` or ``Integral.transform`` are now used for other
    changes of summation or integration variables.  Support for finite and
    infinite sequence products has also been restored.

  - In addition there were a number of fixes to the evaluation of nested sums
    and sums involving Kronecker delta functions, see :sympyissue:`7023`
    and :sympyissue:`7086`.

* Series

  - The ``Order`` object used to represent the growth of a function in series expansions
    as a variable tend to zero can now also represent growth as a variable tend to infinity.
    This also fixed a number of issues with limits.  See :sympyissue:`3333` and :sympyissue:`5769`.

  - Division by ``Order`` is disallowed, see :sympyissue:`4855`.

  - Addition of ``Order`` object is now commutative, see :sympyissue:`4279`.

* Physics

  - Initial work on gamma matrices, depending on the tensor module.

* Logic

  - New objects ``true`` and ``false`` which are ``Basic`` versions of the
    Python builtins ``True`` and ``False``.

* Other

  - Arbitrary comparisons between expressions (like ``x < y``) no longer have a
    boolean truth value. This means code like ``if x < y`` or ``sorted(exprs)`` will
    raise ``TypeError`` if ``x < y`` is symbolic.  A typical fix of the former is ``if
    (x < y) is True`` (assuming the ``if`` block should be skipped if ``x < y`` is
    symbolic), and of the latter is ``sorted(exprs, key=default_sort_key)``, which
    will order the expressions in an arbitrary, but consistent way, even across
    platforms and Python versions.  See :sympyissue:`5931`.

  - Arbitrary comparisons between complex *numbers* (for example,
    ``I > 1``) now raise ``TypeError`` as well (see :sympypull:`2510`).

  - ``minimal_polynomial`` now works with algebraic functions, like ``minimal_polynomial(sqrt(x) + sqrt(x + 1), y)``.

  - ``exp`` can now act on any matrix, even those which are not diagonalizable. It
    is also more comfortable to call it, ``exp(m)`` instead of just ``m.exp()``, as
    was required previously.

  - ``sympify`` now has an option ``evaluate=False`` that will not automatically simplify expressions like ``x+x``.

  - Deep processing of ``cancel`` and ``simplify`` functions.
    ``simplify`` is now recursive through the expression tree.
    See e.g. :sympyissue:`7022`.

  - Improved the modularity of the codebase for potential subclasses,
    see :sympyissue:`6751`.

  - The SymPy cheatsheet was cleaned up.

Compatibility breaks
====================

- Removed deprecated Real class and is_Real property of Basic, see :sympyissue:`4820`.
- Removed deprecated 'each_char' option for ``symbols()``, see :sympyissue:`5018`.
- The ``viewer="StringIO"`` option to ``preview()`` has been deprecated.  Use
  ``viewer="BytesIO"`` instead. See :sympyissue:`7083`.
- ``TransformationSet`` has been renamed to ``ImageSet``.  Added public facing ``imageset`` function.
