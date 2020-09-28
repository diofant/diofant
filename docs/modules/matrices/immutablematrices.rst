Immutable Matrices
==================

The standard Matrix class in Diofant is mutable. This is important for
performance reasons but means that standard matrices cannot interact well with
the rest of Diofant. This is because the :class:`~diofant.core.basic.Basic` object, from which most
Diofant classes inherit, is immutable.

The mission of the :class:`~diofant.matrices.immutable.ImmutableMatrix` class is to bridge the tension
between performance/mutability and safety/immutability. Immutable matrices can
do almost everything that normal matrices can do but they inherit from
:class:`~diofant.core.basic.Basic` and can thus interact more naturally with the rest of Diofant.
:class:`~diofant.matrices.immutable.ImmutableMatrix` also inherits from :class:`~diofant.matrices.expressions.MatrixExpr`, allowing it to
interact freely with Diofant's Matrix Expression module.

You can turn any Matrix-like object into an :class:`~diofant.matrices.immutable.ImmutableMatrix` by calling
the constructor

    >>> M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> M[1, 1] = 0
    >>> IM = ImmutableMatrix(M)
    >>> IM
    Matrix([
    [1, 2, 3],
    [4, 0, 6],
    [7, 8, 9]])
    >>> IM[1, 1] = 5
    Traceback (most recent call last):
    ...
    TypeError: Can not set values in Immutable Matrix. Use Matrix instead.


ImmutableMatrix Class Reference
-------------------------------

.. module:: diofant.matrices.immutable

.. autoclass:: ImmutableMatrix
   :members:
