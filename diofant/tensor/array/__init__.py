r"""
N-dim array module.

Four classes are provided to handle N-dim arrays, given by the combinations
dense/sparse (i.e. whether to store all elements or only the non-zero ones in
memory) and mutable/immutable (immutable classes are Diofant objects, but cannot
change after they have been created).

Examples
========

The following examples show the usage of ``Array``. This is an abbreviation for
``ImmutableDenseNDimArray``, that is an immutable and dense N-dim array, the
other classes are analogous. For mutable classes it is also possible to change
element values after the object has been constructed.

Array construction can detect the shape of nested lists and tuples:

>>> a1 = Array([[1, 2], [3, 4], [5, 6]])
>>> a1
[[1, 2], [3, 4], [5, 6]]
>>> a1.shape
(3, 2)
>>> a1.rank()
2
>>> a2 = Array([[[x, y], [z, x*z]], [[1, x*y], [1/x, x/y]]])
>>> a2
[[[x, y], [z, x*z]], [[1, x*y], [1/x, x/y]]]
>>> a2.shape
(2, 2, 2)
>>> a2.rank()
3

Otherwise one could pass a 1-dim array followed by a shape tuple:

>>> m1 = Array(range(12), (3, 4))
>>> m1
[[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]
>>> m2 = Array(range(12), (3, 2, 2))
>>> m2
[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]]
>>> m2[1, 1, 1]
7
>>> m2.reshape(4, 3)
[[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]

Slice support:

>>> m2[:, 1, 1]
[3, 7, 11]

Elementwise derivative:

>>> m3 = Array([x**3, x*y, z])
>>> m3.diff(x)
[3*x**2, y, 0]
>>> m3.diff(z)
[0, 0, 1]

Multiplication with other Diofant expressions is applied elementwisely:

>>> (1+x)*m3
[x**3*(x + 1), x*y*(x + 1), z*(x + 1)]

To apply a function to each element of the N-dim array, use ``applyfunc``:

>>> m3.applyfunc(lambda x: x/2)
[x**3/2, x*y/2, z/2]

N-dim arrays can be converted to nested lists by the ``tolist()`` method:

>>> m2.tolist()
[[[0, 1], [2, 3]], [[4, 5], [6, 7]], [[8, 9], [10, 11]]]

If the rank is 2, it is possible to convert them to matrices
with ``tomatrix()``:

>>> m1.tomatrix()
Matrix([
[0, 1,  2,  3],
[4, 5,  6,  7],
[8, 9, 10, 11]])

Products and contractions
-------------------------

Tensor product between arrays `A_{i_1,\ldots,i_n}` and `B_{j_1,\ldots,j_m}`
creates the combined array `P = A \otimes B` defined as

.. math::
    P_{i_1,\ldots,i_n,j_1,\ldots,j_m} :=
    A_{i_1,\ldots,i_n}\cdot B_{j_1,\ldots,j_m}

It is available through ``tensorproduct(...)``:

>>> A = Array([x, y, z, t])
>>> B = Array([1, 2, 3, 4])
>>> tensorproduct(A, B)
[[x, 2*x, 3*x, 4*x], [y, 2*y, 3*y, 4*y], [z, 2*z, 3*z, 4*z],
 [t, 2*t, 3*t, 4*t]]

Tensor product between a rank-1 array and a matrix creates a rank-3 array:

>>> p1 = tensorproduct(A, eye(4))
>>> p1
[[[x, 0, 0, 0], [0, x, 0, 0], [0, 0, x, 0], [0, 0, 0, x]],
 [[y, 0, 0, 0], [0, y, 0, 0], [0, 0, y, 0], [0, 0, 0, y]],
 [[z, 0, 0, 0], [0, z, 0, 0], [0, 0, z, 0], [0, 0, 0, z]],
 [[t, 0, 0, 0], [0, t, 0, 0], [0, 0, t, 0], [0, 0, 0, t]]]

Now, to get back `A_0 \otimes \mathbf{1}` one can access `p_{0,m,n}` by slicing:

>>> p1[0, :, :]
[[x, 0, 0, 0], [0, x, 0, 0], [0, 0, x, 0], [0, 0, 0, x]]

Tensor contraction sums over the specified axes, for example contracting
positions `a` and `b` means

.. math ::
    A_{i_1,\ldots,i_a,\ldots,i_b,\ldots,i_n}
     \implies \sum_k A_{i_1,\ldots,k,\ldots,k,\ldots,i_n}

Remember that Python indexing is zero starting, to contract the a-th and b-th
axes it is therefore necessary to specify `a-1` and `b-1`

>>> C = Array([[x, y], [z, t]])

The matrix trace is equivalent to the contraction of a rank-2 array:

`A_{m,n} \implies \sum_k A_{k,k}`

>>> tensorcontraction(C, (0, 1))
t + x

Matrix product is equivalent to a tensor product of two rank-2 arrays,
followed by a contraction of the 2nd and 3rd axes (in Python indexing
axes number 1, 2).

`A_{m,n}\cdot B_{i,j} \implies \sum_k A_{m, k}\cdot B_{k, j}`

>>> D = Array([[2, 1], [0, -1]])
>>> tensorcontraction(tensorproduct(C, D), (1, 2))
[[2*x, x - y], [2*z, -t + z]]

One may verify that the matrix product is equivalent:

>>> Matrix([[x, y], [z, t]])*Matrix([[2, 1], [0, -1]])
Matrix([
[2*x,  x - y],
[2*z, -t + z]])

or equivalently

>>> C.tomatrix()*D.tomatrix()
Matrix([
[2*x,  x - y],
[2*z, -t + z]])


Derivatives by array
--------------------

The usual derivative operation may be extended to support derivation with
respect to arrays, provided that all elements in the that array are symbols or
expressions suitable for derivations.

The definition of a derivative by an array is as follows: given the array
`A_{i_1, \ldots, i_N}` and the array `X_{j_1, \ldots, j_M}`
the derivative of arrays will return a new array `B` defined by

`B_{j_1,\ldots,j_M,i_1,\ldots,i_N} := \frac{\partial A_{i_1,\ldots,i_N}}{\partial X_{j_1,\ldots,j_M}}`

The function ``derive_by_array`` performs such an operation.  With scalars,
it behaves exactly as the ordinary derivative:

>>> derive_by_array(sin(x*y), x)
y*cos(x*y)

Scalar derived by an array basis:

>>> derive_by_array(sin(x*y), [x, y, z])
[y*cos(x*y), x*cos(x*y), 0]

Deriving array by an array basis: `B^{nm} := \frac{\partial A^m}{\partial x^n}`

>>> basis = [x, y, z]
>>> ax = derive_by_array([exp(x), sin(y*z), t], basis)
>>> ax
[[E**x, 0, 0], [0, z*cos(y*z), 0], [0, y*cos(y*z), 0]]

Contraction of the resulting array: `\sum_m \frac{\partial A^m}{\partial x^m}`

>>> tensorcontraction(ax, (0, 1))
E**x + z*cos(y*z)

"""

from .arrayop import (derive_by_array, permutedims, tensorcontraction,
                      tensorproduct)
from .dense_ndim_array import (DenseNDimArray, ImmutableDenseNDimArray,
                               MutableDenseNDimArray)
from .ndim_array import NDimArray
from .sparse_ndim_array import (ImmutableSparseNDimArray,
                                MutableSparseNDimArray, SparseNDimArray)


Array = ImmutableDenseNDimArray
