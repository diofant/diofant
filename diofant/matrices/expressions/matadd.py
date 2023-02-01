import functools
import operator

from ...core import Add, Expr
from ...core.logic import _fuzzy_group
from ...core.strategies import (condition, do_one, exhaust, flatten, glom,
                                rm_id, sort, unpack)
from ...core.sympify import sympify
from ...functions import adjoint
from ...utilities import default_sort_key, sift
from ..matrices import MatrixBase, ShapeError
from .matexpr import MatrixExpr, ZeroMatrix
from .transpose import transpose


class MatAdd(MatrixExpr):
    """A Sum of Matrix Expressions

    MatAdd inherits from and operates like Diofant Add

    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> C = MatrixSymbol('C', 5, 5)
    >>> MatAdd(A, B, C)
    A + B + C
    """

    is_MatAdd = True

    def _eval_is_commutative(self):
        return _fuzzy_group((a.is_commutative for a in self.args),
                            quick_exit=True)

    def __new__(cls, *args, **kwargs):
        args = list(map(sympify, args))
        check = kwargs.get('check', True)

        obj = Expr.__new__(cls, *args)
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
        raise TypeError('Mix of Matrix and Scalar symbols')

    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError(f'Matrices {A} and {B} are not aligned')


def factor_of(arg):
    return arg.as_coeff_mmul()[0]


def matrix_of(arg):
    return unpack(arg.as_coeff_mmul()[1])


def combine(cnt, mat):
    if cnt == 1:
        return mat
    return cnt * mat


def merge_explicit(matadd):
    """Merge explicit MatrixBase arguments

    >>> A = MatrixSymbol('A', 2, 2)
    >>> B = eye(2)
    >>> C = Matrix([[1, 2], [3, 4]])
    >>> X = MatAdd(A, B, C)
    >>> pprint(X, use_unicode=False)
        [1  0]   [1  2]
    A + [    ] + [    ]
        [0  1]   [3  4]
    >>> pprint(merge_explicit(X), use_unicode=False)
        [2  2]
    A + [    ]
        [3  5]
    """
    groups = sift(matadd.args, lambda arg: isinstance(arg, MatrixBase))
    if len(groups[True]) > 1:
        return MatAdd(*(groups[False] + [functools.reduce(operator.add, groups[True])]))
    return matadd


rules = (rm_id(lambda x: x == 0 or isinstance(x, ZeroMatrix)),
         unpack,
         flatten,
         glom(matrix_of, factor_of, combine),
         merge_explicit,
         sort(default_sort_key))

canonicalize = exhaust(condition(lambda x: isinstance(x, MatAdd),
                                 do_one(rules)))
