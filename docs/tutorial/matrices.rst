==========
 Matrices
==========

..
    >>> init_printing(pretty_print=True, use_unicode=True)

To make a matrix in Diofant, use the
:class:`~diofant.matrices.Matrix` object.  A matrix is
constructed by providing a list of row vectors that make up the
matrix.

    >>> Matrix([[1, -1], [3, 4], [0, 2]])
    ⎡1  -1⎤
    ⎢     ⎥
    ⎢3  4 ⎥
    ⎢     ⎥
    ⎣0  2 ⎦

A list of elements is considered to be a column vector.

    >>> Matrix([1, 2, 3])
    ⎡1⎤
    ⎢ ⎥
    ⎢2⎥
    ⎢ ⎥
    ⎣3⎦

One important thing to note about Diofant matrices is that, unlike
every other object in Diofant, they are mutable.  This means that they
can be modified in place, as we will see below.  Use
:class:`~diofant.matrices.immutable.ImmutableMatrix` in places that
require immutability, such as inside other Diofant expressions or as
keys to dictionaries.

Indexing
========

Diofant matrices support subscription of matrix elements with pair of
integers or :class:`slice` instances.  In last case, new
:class:`~diofant.matrices.Matrix` instances will be returned.

    >>> M = Matrix([[1, 2, 3], [4, 5, 6]])
    >>> M[0, 1]
    2
    >>> M[1, 1]
    5
    >>> M[:, 1]
    ⎡2⎤
    ⎢ ⎥
    ⎣5⎦
    >>> M[1, :-1]
    [4  5]
    >>> M[0, :]
    [1  2  3]
    >>> M[:, -1]
    ⎡3⎤
    ⎢ ⎥
    ⎣6⎦

It's possible to modify matrix elements.

    >>> M[0, 0] = 0
    >>> M
    ⎡0  2  3⎤
    ⎢       ⎥
    ⎣4  5  6⎦
    >>> M[1, 1:] = Matrix([[0, 0]])
    >>> M
    ⎡0  2  3⎤
    ⎢       ⎥
    ⎣4  0  0⎦

Reshape and Rearrange
=====================

To get the shape of a matrix use
:attr:`~diofant.matrices.matrices.MatrixBase.shape` property

    >>> M = Matrix([[1, 2, 3], [-2, 0, 4]])
    >>> M
    ⎡1   2  3⎤
    ⎢        ⎥
    ⎣-2  0  4⎦
    >>> M.shape
    (2, 3)

To delete a row or column, use :keyword:`del`

    >>> del M[:, 0]
    >>> M
    ⎡2  3⎤
    ⎢    ⎥
    ⎣0  4⎦
    >>> del M[1, :]
    >>> M
    [2  3]

To insert rows or columns, use methods
:meth:`~diofant.matrices.matrices.MatrixBase.row_insert` or
:meth:`~diofant.matrices.matrices.MatrixBase.col_insert`.

    >>> M
    [2  3]
    >>> M = M.row_insert(1, Matrix([[0, 4]]))
    >>> M
    ⎡2  3⎤
    ⎢    ⎥
    ⎣0  4⎦
    >>> M = M.col_insert(0, Matrix([1, -2]))
    >>> M
    ⎡1   2  3⎤
    ⎢        ⎥
    ⎣-2  0  4⎦

.. note::

   You can see, that these methods will modify the Matrix **in
   place**.  In general, as a rule, such methods will return ``None``.

To swap two given rows or columns, use methods
:meth:`~diofant.matrices.dense.MutableDenseMatrix.row_swap` or
:meth:`~diofant.matrices.dense.MutableDenseMatrix.col_swap`.

    >>> M.row_swap(0, 1)
    >>> M
    ⎡-2  0  4⎤
    ⎢        ⎥
    ⎣1   2  3⎦
    >>> M.col_swap(1, 2)
    >>> M
    ⎡-2  4  0⎤
    ⎢        ⎥
    ⎣1   3  2⎦

To take the transpose of a Matrix, use
:attr:`~diofant.matrices.matrices.MatrixBase.T` property.

    >>> M.T
    ⎡-2  1⎤
    ⎢     ⎥
    ⎢4   3⎥
    ⎢     ⎥
    ⎣0   2⎦

Algebraic Operations
====================

Simple operations like addition and multiplication are done just by
using ``+``, ``*``, and ``**``.  To find the inverse of a matrix, just
raise it to the ``-1`` power.

    >>> M, N = Matrix([[1, 3], [-2, 3]]), Matrix([[0, 3], [0, 7]])
    >>> M + N
    ⎡1   6 ⎤
    ⎢      ⎥
    ⎣-2  10⎦
    >>> M*N
    ⎡0  24⎤
    ⎢     ⎥
    ⎣0  15⎦
    >>> 3*M
    ⎡3   9⎤
    ⎢     ⎥
    ⎣-6  9⎦
    >>> M**2
    ⎡-5  12⎤
    ⎢      ⎥
    ⎣-8  3 ⎦
    >>> M**-1
    ⎡1/3  -1/3⎤
    ⎢         ⎥
    ⎣2/9  1/9 ⎦
    >>> N**-1
    Traceback (most recent call last):
    ...
    ValueError: Matrix det == 0; not invertible.

Special Matrices
=================

Several constructors exist for creating common matrices.  To create an
identity matrix, use :func:`~diofant.matrices.dense.eye` function.

    >>> eye(3)
    ⎡1  0  0⎤
    ⎢       ⎥
    ⎢0  1  0⎥
    ⎢       ⎥
    ⎣0  0  1⎦
    >>> eye(4)
    ⎡1  0  0  0⎤
    ⎢          ⎥
    ⎢0  1  0  0⎥
    ⎢          ⎥
    ⎢0  0  1  0⎥
    ⎢          ⎥
    ⎣0  0  0  1⎦

To create a matrix of all zeros, use
:func:`~diofant.matrices.dense.zeros` function.

    >>> zeros(2, 3)
    ⎡0  0  0⎤
    ⎢       ⎥
    ⎣0  0  0⎦

Similarly, function :func:`~diofant.matrices.dense.ones` creates a
matrix of ones.

    >>> ones(3, 2)
    ⎡1  1⎤
    ⎢    ⎥
    ⎢1  1⎥
    ⎢    ⎥
    ⎣1  1⎦

To create diagonal matrices, use function
:func:`~diofant.matrices.dense.diag`.  Its arguments can be either
numbers or matrices.  A number is interpreted as a `1\times 1`
matrix. The matrices are stacked diagonally.

    >>> diag(1, 2, 3)
    ⎡1  0  0⎤
    ⎢       ⎥
    ⎢0  2  0⎥
    ⎢       ⎥
    ⎣0  0  3⎦
    >>> diag(-1, ones(2, 2), Matrix([5, 7, 5]))
    ⎡-1  0  0  0⎤
    ⎢           ⎥
    ⎢0   1  1  0⎥
    ⎢           ⎥
    ⎢0   1  1  0⎥
    ⎢           ⎥
    ⎢0   0  0  5⎥
    ⎢           ⎥
    ⎢0   0  0  7⎥
    ⎢           ⎥
    ⎣0   0  0  5⎦

Advanced Methods
================

To compute the determinant of a matrix, use
:meth:`~diofant.matrices.matrices.MatrixBase.det` method.

    >>> Matrix([[1, 0, 1], [2, -1, 3], [4, 3, 2]])
    ⎡1  0   1⎤
    ⎢        ⎥
    ⎢2  -1  3⎥
    ⎢        ⎥
    ⎣4  3   2⎦
    >>> _.det()
    -1

To put a matrix into reduced row echelon form, use method
:meth:`~diofant.matrices.matrices.MatrixBase.rref`.  It returns a
tuple of two elements.  The first is the reduced row echelon form, and
the second is a list of indices of the pivot columns.

    >>> Matrix([[1, 0, 1, 3], [2, 3, 4, 7], [-1, -3, -3, -4]])
    ⎡1   0   1   3 ⎤
    ⎢              ⎥
    ⎢2   3   4   7 ⎥
    ⎢              ⎥
    ⎣-1  -3  -3  -4⎦
    >>> _.rref()
    ⎛⎡1  0   1    3 ⎤        ⎞
    ⎜⎢              ⎥        ⎟
    ⎜⎢0  1  2/3  1/3⎥, [0, 1]⎟
    ⎜⎢              ⎥        ⎟
    ⎝⎣0  0   0    0 ⎦        ⎠

To find the nullspace of a matrix, use method
:meth:`~diofant.matrices.matrices.MatrixBase.nullspace`.  It returns a
list of column vectors that span the nullspace of the matrix.

    >>> Matrix([[1, 2, 3, 0, 0], [4, 10, 0, 0, 1]])
    ⎡1  2   3  0  0⎤
    ⎢              ⎥
    ⎣4  10  0  0  1⎦
    >>> _.nullspace()
    ⎡⎡-15⎤  ⎡0⎤  ⎡ 1  ⎤⎤
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 6 ⎥  ⎢0⎥  ⎢-1/2⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 1 ⎥, ⎢0⎥, ⎢ 0  ⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎢⎢ 0 ⎥  ⎢1⎥  ⎢ 0  ⎥⎥
    ⎢⎢   ⎥  ⎢ ⎥  ⎢    ⎥⎥
    ⎣⎣ 0 ⎦  ⎣0⎦  ⎣ 1  ⎦⎦

To find the eigenvalues of a matrix, use method
:meth:`~diofant.matrices.matrices.MatrixBase.eigenvals`.  It returns a
dictionary of roots including its multiplicity (similar to the output
of :func:`~diofant.polys.polyroots.roots` function).

    >>> M = Matrix([[3, -2, 4, -2], [5, +3, -3, -2],
    ...             [5, -2, 2, -2], [5, -2, -3, +3]])
    >>> M
    ⎡3  -2  4   -2⎤
    ⎢             ⎥
    ⎢5  3   -3  -2⎥
    ⎢             ⎥
    ⎢5  -2  2   -2⎥
    ⎢             ⎥
    ⎣5  -2  -3  3 ⎦
    >>> M.eigenvals()
    {-2: 1, 3: 1, 5: 2}

This means that ``M`` has eigenvalues -2, 3, and 5, and that the
eigenvalues -2 and 3 have algebraic multiplicity 1 and that the
eigenvalue 5 has algebraic multiplicity 2.

Matrices can have symbolic elements.

    >>> Matrix([[x, x + y], [y, x]])
    ⎡x  x + y⎤
    ⎢        ⎥
    ⎣y    x  ⎦
    >>> _.eigenvals()
    ⎧      ___________           ___________   ⎫
    ⎨x - ╲╱ y⋅(x + y) : 1, x + ╲╱ y⋅(x + y) : 1⎬
    ⎩                                          ⎭

To find the eigenvectors of a matrix, use method
:meth:`~diofant.matrices.matrices.MatrixBase.eigenvects`.

    >>> M.eigenvects()
    ⎡⎛       ⎡⎡0⎤⎤⎞  ⎛      ⎡⎡1⎤⎤⎞  ⎛      ⎡⎡1⎤  ⎡0 ⎤⎤⎞⎤
    ⎢⎜       ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥  ⎢  ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥  ⎢-1⎥⎥⎟⎥
    ⎢⎜-2, 1, ⎢⎢ ⎥⎥⎟, ⎜3, 1, ⎢⎢ ⎥⎥⎟, ⎜5, 2, ⎢⎢ ⎥, ⎢  ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥⎥⎟  ⎜      ⎢⎢1⎥  ⎢0 ⎥⎥⎟⎥
    ⎢⎜       ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥⎥⎟  ⎜      ⎢⎢ ⎥  ⎢  ⎥⎥⎟⎥
    ⎣⎝       ⎣⎣1⎦⎦⎠  ⎝      ⎣⎣1⎦⎦⎠  ⎝      ⎣⎣0⎦  ⎣1 ⎦⎦⎠⎦

This shows us that, for example, the eigenvalue 5 also has geometric
multiplicity 2, because it has two eigenvectors.  Because the
algebraic and geometric multiplicities are the same for all the
eigenvalues, ``M`` is diagonalizable.

To diagonalize a matrix, use method
:meth:`~diofant.matrices.matrices.MatrixBase.diagonalize`.  It returns
a tuple `(P, D)`, where `D` is diagonal and `M = PDP^{-1}`.

    >>> M.diagonalize()
    ⎛⎡0  1  1  0 ⎤  ⎡-2  0  0  0⎤⎞
    ⎜⎢           ⎥  ⎢           ⎥⎟
    ⎜⎢1  1  1  -1⎥  ⎢0   3  0  0⎥⎟
    ⎜⎢           ⎥, ⎢           ⎥⎟
    ⎜⎢1  1  1  0 ⎥  ⎢0   0  5  0⎥⎟
    ⎜⎢           ⎥  ⎢           ⎥⎟
    ⎝⎣1  1  0  1 ⎦  ⎣0   0  0  5⎦⎠
    >>> _[0]*_[1]*_[0]**-1 == M
    True

If all you want is the characteristic polynomial, use method
:meth:`~diofant.matrices.matrices.MatrixBase.charpoly`.  This is more
efficient than :meth:`~diofant.matrices.matrices.MatrixBase.eigenvals`
method, because sometimes symbolic roots can be expensive to
calculate.

    >>> M.charpoly(x)
    PurePoly(x**4 - 11*x**3 + 29*x**2 + 35*x - 150, x, domain='ZZ')
    >>> factor(_)
           2
    (x - 5) ⋅(x - 3)⋅(x + 2)

To compute Jordan canonical form `J` for matrix `M` and its similarity
transformation `P` (i.e. such that `J = P M P^{-1}`), use method
:meth:`~diofant.matrices.matrices.MatrixBase.jordan_form`.

    >>> Matrix([[-2, 4], [1, 3]]).jordan_form()
    ⎛⎡      ____              ⎤  ⎡    -4            -4      ⎤⎞
    ⎜⎢1   ╲╱ 41               ⎥  ⎢────────────  ────────────⎥⎟
    ⎜⎢─ + ──────       0      ⎥  ⎢    ____              ____⎥⎟
    ⎜⎢2     2                 ⎥  ⎢  ╲╱ 41    5    5   ╲╱ 41 ⎥⎟
    ⎜⎢                        ⎥, ⎢- ────── - ─  - ─ + ──────⎥⎟
    ⎜⎢                ____    ⎥  ⎢    2      2    2     2   ⎥⎟
    ⎜⎢              ╲╱ 41    1⎥  ⎢                          ⎥⎟
    ⎜⎢    0       - ────── + ─⎥  ⎣     1             1      ⎦⎟
    ⎝⎣                2      2⎦                              ⎠
