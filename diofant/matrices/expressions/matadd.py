from functools import reduce
from operator import add

from strategies import condition, exhaust, do_one

from diofant.core import Add, Basic, sympify
from diofant.functions import adjoint
from diofant.matrices.matrices import MatrixBase
from diofant.matrices.expressions.transpose import transpose
from diofant.core.strategies import rm_id, unpack, flatten, sort, glom
from diofant.matrices.expressions.matexpr import MatrixExpr, ShapeError, ZeroMatrix
from diofant.utilities import default_sort_key, sift


class MatAdd(MatrixExpr):
    """A Sum of Matrix Expressions

    MatAdd inherits from and operates like Diofant Add

    >>> from diofant import MatAdd, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> C = MatrixSymbol('C', 5, 5)
    >>> MatAdd(A, B, C)
    A + B + C
    """
    is_MatAdd = True

    def __new__(cls, *args, **kwargs):
        args = list(map(sympify, args))
        check = kwargs.get('check', True)

        obj = Basic.__new__(cls, *args)
        if check:
            validate(*args)
        return obj

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j):
        return Add(*[arg._entry(i, j) for arg in self.args])

    def _eval_transpose(self):
        return MatAdd(*[transpose(arg) for arg in self.args]).doit()

    def _eval_adjoint(self):
        return MatAdd(*[adjoint(arg) for arg in self.args]).doit()

    def _eval_trace(self):
        from .trace import trace
        return Add(*[trace(arg) for arg in self.args]).doit()

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return canonicalize(MatAdd(*args))


def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")

    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned" % (A, B))


def factor_of(arg):
    return arg.as_coeff_mmul()[0]


def matrix_of(arg):
    return unpack(arg.as_coeff_mmul()[1])


def combine(cnt, mat):
    if cnt == 1:
        return mat
    else:
        return cnt * mat


def merge_explicit(matadd):
    """ Merge explicit MatrixBase arguments

    >>> from diofant import MatrixSymbol, eye, Matrix, MatAdd, pprint
    >>> from diofant.matrices.expressions.matadd import merge_explicit
    >>> A = MatrixSymbol('A', 2, 2)
    >>> B = eye(2)
    >>> C = Matrix([[1, 2], [3, 4]])
    >>> X = MatAdd(A, B, C)
    >>> pprint(X, use_unicode=False)
    A + [1  0] + [1  2]
        [    ]   [    ]
        [0  1]   [3  4]
    >>> pprint(merge_explicit(X), use_unicode=False)
    A + [2  2]
        [    ]
        [3  5]
    """
    groups = sift(matadd.args, lambda arg: isinstance(arg, MatrixBase))
    if len(groups[True]) > 1:
        return MatAdd(*(groups[False] + [reduce(add, groups[True])]))
    else:
        return matadd


rules = (rm_id(lambda x: x == 0 or isinstance(x, ZeroMatrix)),
         unpack,
         flatten,
         glom(matrix_of, factor_of, combine),
         merge_explicit,
         sort(default_sort_key))

canonicalize = exhaust(condition(lambda x: isinstance(x, MatAdd),
                                 do_one(rules)))
