from strategies import do_one, exhaust
from strategies.core import typed

from ...core import Add, Expr, Mul, Number, sympify
from ...core.logic import _fuzzy_group
from ...core.strategies import flatten, rm_id, unpack
from ...functions import adjoint
from ..matrices import MatrixBase, ShapeError
from .matexpr import Identity, MatrixExpr, ZeroMatrix
from .transpose import transpose


class MatMul(MatrixExpr):
    """
    A product of matrix expressions

    Examples
    ========

    >>> A = MatrixSymbol('A', 5, 4)
    >>> B = MatrixSymbol('B', 4, 3)
    >>> C = MatrixSymbol('C', 3, 6)
    >>> MatMul(A, B, C)
    A*B*C

    """

    is_MatMul = True

    def _eval_is_commutative(self):
        return _fuzzy_group(a.is_commutative for a in self.args)

    def __new__(cls, *args, **kwargs):
        check = kwargs.get('check', True)

        args = list(map(sympify, args))
        obj = Expr.__new__(cls, *args)
        factor, matrices = obj.as_coeff_matrices()
        if check:
            validate(*matrices)
        return obj

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return matrices[0].rows, matrices[-1].cols

    def _entry(self, i, j, expand=True):
        coeff, matrices = self.as_coeff_matrices()

        if len(matrices) == 1:  # situation like 2*X, matmul is just X
            return coeff * matrices[0][i, j]

        head, tail = matrices[0], matrices[1:]
        X = head
        Y = MatMul(*tail)

        from ...core import Dummy
        from ...concrete import Sum
        from .. import ImmutableMatrix
        k = Dummy('k', integer=True)
        if X.has(ImmutableMatrix) or Y.has(ImmutableMatrix):
            return coeff*Add(*[X[i, k]*Y[k, j] for k in range(X.cols)])
        result = Sum(coeff*X[i, k]*Y[k, j], (k, 0, X.cols - 1))
        if not X.cols.is_number:
            # Don't waste time in result.doit() if the sum bounds are symbolic
            expand = False
        return result.doit() if expand else result

    def as_coeff_matrices(self):
        scalars = [x for x in self.args if not x.is_Matrix]
        matrices = [x for x in self.args if x.is_Matrix]
        coeff = Mul(*scalars)

        return coeff, matrices

    def as_coeff_mmul(self):
        coeff, matrices = self.as_coeff_matrices()
        return coeff, MatMul(*matrices)

    def _eval_transpose(self):
        return MatMul(*[transpose(arg) for arg in self.args[::-1]]).doit()

    def _eval_adjoint(self):
        return MatMul(*[adjoint(arg) for arg in self.args[::-1]]).doit()

    def _eval_trace(self):
        factor, mmul = self.as_coeff_mmul()
        if factor != 1:
            from .trace import trace
            return factor * trace(mmul.doit())
        else:
            raise NotImplementedError("Can't simplify any further")

    def _eval_determinant(self):
        from .determinant import Determinant
        factor, matrices = self.as_coeff_matrices()
        square_matrices = only_squares(*matrices)
        return factor**self.rows * Mul(*list(map(Determinant, square_matrices)))

    def _eval_inverse(self):
        try:
            return MatMul(*[
                arg.inverse() if isinstance(arg, MatrixExpr) else arg**-1
                for arg in self.args[::-1]]).doit()
        except ShapeError:
            from .inverse import Inverse
            return Inverse(self)

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return canonicalize(MatMul(*args))


def validate(*matrices):
    """Checks for valid shapes for args of MatMul."""
    for i in range(len(matrices)-1):
        A, B = matrices[i:i+2]
        if A.cols != B.rows:
            raise ShapeError("Matrices %s and %s are not aligned" % (A, B))

# Rules


def newmul(*args):
    if args[0] == 1:
        args = args[1:]
    return MatMul(*args)


def any_zeros(mul):
    if any(arg.is_zero or (arg.is_Matrix and arg.is_ZeroMatrix)
           for arg in mul.args):
        matrices = [arg for arg in mul.args if arg.is_Matrix]
        return ZeroMatrix(matrices[0].rows, matrices[-1].cols)
    return mul


def merge_explicit(matmul):
    """Merge explicit MatrixBase arguments

    >>> A = MatrixSymbol('A', 2, 2)
    >>> B = Matrix([[1, 1], [1, 1]])
    >>> C = Matrix([[1, 2], [3, 4]])
    >>> X = MatMul(A, B, C)
    >>> pprint(X, use_unicode=False)
      [1  1] [1  2]
    A*[    ]*[    ]
      [1  1] [3  4]
    >>> pprint(merge_explicit(X), use_unicode=False)
      [4  6]
    A*[    ]
      [4  6]

    >>> X = MatMul(B, A, C)
    >>> pprint(X, use_unicode=False)
    [1  1]   [1  2]
    [    ]*A*[    ]
    [1  1]   [3  4]
    >>> pprint(merge_explicit(X), use_unicode=False)
    [1  1]   [1  2]
    [    ]*A*[    ]
    [1  1]   [3  4]

    """
    if not any(isinstance(arg, MatrixBase) for arg in matmul.args):
        return matmul
    newargs = []
    last = matmul.args[0]
    for arg in matmul.args[1:]:
        if isinstance(arg, (MatrixBase, Number)) and isinstance(last, (MatrixBase, Number)):
            last = last * arg
        else:
            newargs.append(last)
            last = arg
    newargs.append(last)

    return MatMul(*newargs)


def xxinv(mul):
    """X * X.inverse() -> Identity."""
    factor, matrices = mul.as_coeff_matrices()
    for i, (X, Y) in enumerate(zip(matrices[:-1], matrices[1:])):
        try:
            if X.is_square and Y.is_square and X == Y.inverse():
                I = Identity(X.rows)
                return newmul(factor, *(matrices[:i] + [I] + matrices[i+2:]))
        except ValueError:  # Y might not be invertible
            pass

    return mul


def remove_ids(mul):
    """ Remove Identities from a MatMul

    This is a modified version of diofant.core.strategies.rm_id.
    This is necesssary because MatMul may contain both MatrixExprs and Exprs
    as args.

    See Also
    ========

    diofant.core.strategies.rm_id

    """
    # Separate Exprs from MatrixExprs in args
    factor, mmul = mul.as_coeff_mmul()
    # Apply standard rm_id for MatMuls
    result = rm_id(lambda x: x.is_Identity is True)(mmul)
    if result != mmul:
        return newmul(factor, *result.args)  # Recombine and return
    else:
        return mul


def factor_in_front(mul):
    factor, matrices = mul.as_coeff_matrices()
    if factor != 1:
        return newmul(factor, *matrices)
    return mul


rules = (any_zeros, remove_ids, xxinv, unpack, rm_id(lambda x: x == 1),
         merge_explicit, factor_in_front, flatten)

canonicalize = exhaust(typed({MatMul: do_one(rules)}))


def only_squares(*matrices):
    """factor matrices only if they are square."""
    if matrices[0].rows != matrices[-1].cols:
        raise RuntimeError("Invalid matrices being multiplied")
    out = []
    start = 0
    for i, M in enumerate(matrices):
        if M.cols == matrices[start].rows:
            out.append(MatMul(*matrices[start:i+1]).doit())
            start = i+1
    return out
