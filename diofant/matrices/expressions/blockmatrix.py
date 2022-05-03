from ...core import Add, Expr, Integer
from ...core.strategies import (bottom_up, condition, do_one, exhaust, typed,
                                unpack)
from ...core.sympify import sympify
from ...logic import false
from ...utilities import sift
from .determinant import Determinant
from .inverse import Inverse
from .matadd import MatAdd
from .matexpr import Identity, MatrixExpr, ZeroMatrix
from .matmul import MatMul
from .slice import MatrixSlice
from .trace import Trace
from .transpose import Transpose, transpose


class BlockMatrix(MatrixExpr):
    """A BlockMatrix is a Matrix composed of other smaller, submatrices

    The submatrices are stored in a Diofant Matrix object but accessed as part of
    a Matrix Expression

    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m, m)
    >>> Z = MatrixSymbol('Z', n, m)
    >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m, n), Y]])
    >>> B
    Matrix([
    [X, Z],
    [0, Y]])

    >>> C = BlockMatrix([[Identity(n), Z]])
    >>> C
    Matrix([[I, Z]])

    >>> block_collapse(C*B)
    Matrix([[X, Z + Z*Y]])

    """

    def __new__(cls, *args):
        from ..immutable import ImmutableMatrix
        args = map(sympify, args)
        mat = ImmutableMatrix(*args)

        obj = Expr.__new__(cls, mat)
        return obj

    @property
    def shape(self):
        numrows = numcols = Integer(0)
        M = self.blocks
        for i in range(M.shape[0]):
            numrows += M[i, 0].shape[0]
        for i in range(M.shape[1]):
            numcols += M[0, i].shape[1]
        return numrows, numcols

    @property
    def blockshape(self):
        return self.blocks.shape

    @property
    def blocks(self):
        return self.args[0]

    @property
    def rowblocksizes(self):
        return [self.blocks[i, 0].rows for i in range(self.blockshape[0])]

    @property
    def colblocksizes(self):
        return [self.blocks[0, i].cols for i in range(self.blockshape[1])]

    def structurally_equal(self, other):
        return (isinstance(other, BlockMatrix)
                and self.shape == other.shape
                and self.blockshape == other.blockshape
                and self.rowblocksizes == other.rowblocksizes
                and self.colblocksizes == other.colblocksizes)

    def _blockmul(self, other):
        if (isinstance(other, BlockMatrix) and
                self.colblocksizes == other.rowblocksizes):
            return BlockMatrix(self.blocks*other.blocks)

        return self * other

    def _blockadd(self, other):
        if (isinstance(other, BlockMatrix)
                and self.structurally_equal(other)):
            return BlockMatrix(self.blocks + other.blocks)

        return self + other

    def _eval_transpose(self):
        from .. import Matrix

        # Flip all the individual matrices
        matrices = [transpose(matrix) for matrix in self.blocks]
        # Make a copy
        M = Matrix(self.blockshape[0], self.blockshape[1], matrices)
        # Transpose the block structure
        M = M.transpose()
        return BlockMatrix(M)

    def _eval_trace(self):
        if self.rowblocksizes == self.colblocksizes:
            return Add(*[Trace(self.blocks[i, i])
                         for i in range(self.blockshape[0])])
        raise NotImplementedError("Can't perform trace of irregular "
                                  'blockshape')  # pragma: no cover

    def _eval_determinant(self):
        return Determinant(self)

    def transpose(self):
        """Return transpose of matrix.

        Examples
        ========

        >>> X = MatrixSymbol('X', n, n)
        >>> Y = MatrixSymbol('Y', m, m)
        >>> Z = MatrixSymbol('Z', n, m)
        >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m, n), Y]])
        >>> B.transpose()
        Matrix([
        [X.T,  0],
        [Z.T, Y.T]])
        >>> _.transpose()
        Matrix([
        [X, Z],
        [0, Y]])

        """
        return self._eval_transpose()

    def _entry(self, i, j):
        # Find row entry
        for row_block, numrows in enumerate(self.rowblocksizes):  # pragma: no branch
            if (i < numrows) != false:
                break
            i -= numrows
        for col_block, numcols in enumerate(self.colblocksizes):  # pragma: no branch
            if (j < numcols) != false:
                break
            j -= numcols
        return self.blocks[row_block, col_block][i, j]

    @property
    def is_Identity(self):
        if self.blockshape[0] != self.blockshape[1]:
            return False
        for i in range(self.blockshape[0]):
            for j in range(self.blockshape[1]):
                if i == j and not self.blocks[i, j].is_Identity:
                    return False
                if i != j and not self.blocks[i, j].is_ZeroMatrix:
                    return False
        return True

    @property
    def is_structurally_symmetric(self):
        return self.rowblocksizes == self.colblocksizes

    def equals(self, other):
        if self == other:
            return True
        if isinstance(other, BlockMatrix) and self.blocks == other.blocks:
            return True
        return super().equals(other)


class BlockDiagMatrix(BlockMatrix):
    """
    A BlockDiagMatrix is a BlockMatrix with matrices only along the diagonal

    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m, m)
    >>> BlockDiagMatrix(X, Y)
    Matrix([
    [X, 0],
    [0, Y]])

    """

    def __new__(cls, *mats):
        return Expr.__new__(BlockDiagMatrix, *mats)

    @property
    def diag(self):
        return self.args

    @property
    def blocks(self):
        from ..immutable import ImmutableMatrix
        mats = self.args
        data = [[mats[i] if i == j else ZeroMatrix(mats[i].rows, mats[j].cols)
                 for j in range(len(mats))]
                for i in range(len(mats))]
        return ImmutableMatrix(data)

    @property
    def shape(self):
        return (sum(block.rows for block in self.args),
                sum(block.cols for block in self.args))

    @property
    def blockshape(self):
        n = len(self.args)
        return n, n

    @property
    def rowblocksizes(self):
        return [block.rows for block in self.args]

    @property
    def colblocksizes(self):
        return [block.cols for block in self.args]

    def _eval_inverse(self, expand='ignored'):
        return BlockDiagMatrix(*[mat.inverse() for mat in self.args])

    def _blockmul(self, other):
        if (isinstance(other, BlockDiagMatrix) and
                self.colblocksizes == other.rowblocksizes):
            return BlockDiagMatrix(*[a*b for a, b in zip(self.args, other.args)])
        else:
            return BlockMatrix._blockmul(self, other)

    def _blockadd(self, other):
        if (isinstance(other, BlockDiagMatrix) and
                self.blockshape == other.blockshape and
                self.rowblocksizes == other.rowblocksizes and
                self.colblocksizes == other.colblocksizes):
            return BlockDiagMatrix(*[a + b for a, b in zip(self.args, other.args)])
        else:
            return BlockMatrix._blockadd(self, other)


def block_collapse(expr):
    """Evaluates a block matrix expression

    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m, m)
    >>> Z = MatrixSymbol('Z', n, m)
    >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m, n), Y]])
    >>> B
    Matrix([
    [X, Z],
    [0, Y]])

    >>> C = BlockMatrix([[Identity(n), Z]])
    >>> C
    Matrix([[I, Z]])

    >>> block_collapse(C*B)
    Matrix([[X, Z + Z*Y]])

    """
    def hasbm(expr):
        return isinstance(expr, MatrixExpr) and expr.has(BlockMatrix)
    rule = exhaust(
        bottom_up(exhaust(condition(hasbm, typed(
            {MatAdd: do_one([bc_matadd, bc_block_plus_ident]),
             MatMul: do_one([bc_matmul, bc_dist]),
             Transpose: bc_transpose,
             Inverse: bc_inverse,
             BlockMatrix: do_one([bc_unpack, deblock])})))))
    result = rule(expr)
    return result.doit()


def bc_unpack(expr):
    if expr.blockshape == (1, 1):
        return expr.blocks[0, 0]
    return expr


def bc_matadd(expr):
    args = sift(expr.args, lambda M: isinstance(M, BlockMatrix))
    blocks = args[True]
    if not blocks:
        return expr

    nonblocks = args[False]
    block = blocks[0]
    for b in blocks[1:]:
        block = block._blockadd(b)
    if nonblocks:
        return MatAdd(*nonblocks) + block
    else:
        return block


def bc_block_plus_ident(expr):
    idents = [arg for arg in expr.args if arg.is_Identity]
    if not idents:
        return expr

    blocks = [arg for arg in expr.args if isinstance(arg, BlockMatrix)]
    if (blocks and all(b.structurally_equal(blocks[0]) for b in blocks)
            and blocks[0].is_structurally_symmetric):
        block_id = BlockDiagMatrix(*[Identity(k)
                                     for k in blocks[0].rowblocksizes])
        return MatAdd(block_id * len(idents), *blocks).doit()

    return expr


def bc_dist(expr):
    """Turn  a*[X, Y] into [a*X, a*Y]."""
    factor, mat = expr.as_coeff_mmul()
    if factor != 1 and isinstance(unpack(mat), BlockMatrix):
        B = unpack(mat).blocks
        return BlockMatrix([[factor * B[i, j] for j in range(B.cols)]
                            for i in range(B.rows)])
    return expr


def bc_matmul(expr):
    factor, matrices = expr.as_coeff_matrices()

    i = 0
    while i + 1 < len(matrices):
        A, B = matrices[i:i+2]
        if isinstance(A, BlockMatrix) and isinstance(B, BlockMatrix):
            matrices[i] = A._blockmul(B)
            matrices.pop(i+1)
        elif isinstance(A, BlockMatrix):
            matrices[i] = A._blockmul(BlockMatrix([[B]]))
            matrices.pop(i+1)
        elif isinstance(B, BlockMatrix):
            matrices[i] = BlockMatrix([[A]])._blockmul(B)
            matrices.pop(i+1)
        else:
            i += 1
    return MatMul(factor, *matrices).doit()


def bc_transpose(expr):
    return BlockMatrix(block_collapse(expr.arg).blocks.applyfunc(transpose).T)


def bc_inverse(expr):
    return blockinverse_2x2(Inverse(reblock_2x2(expr.arg)))


def blockinverse_2x2(expr):
    # Cite: The Matrix Cookbook Section 9.1.3
    [[A, B],
     [C, D]] = expr.arg.blocks.tolist()

    return BlockMatrix([[+(A - B*D.inverse()*C).inverse(),  (-A).inverse()*B*(D - C*A.inverse()*B).inverse()],
                        [-(D - C*A.inverse()*B).inverse()*C*A.inverse(),     (D - C*A.inverse()*B).inverse()]])


def deblock(B):
    """Flatten a BlockMatrix of BlockMatrices."""
    if not isinstance(B, BlockMatrix) or not B.blocks.has(BlockMatrix):
        return B

    def wrap(x):
        return x if isinstance(x, BlockMatrix) else BlockMatrix([[x]])

    bb = B.blocks.applyfunc(wrap)  # everything is a block

    from .. import Matrix
    MM = Matrix(0, sum(bb[0, i].blocks.shape[1] for i in range(bb.shape[1])), [])
    for row in range(bb.shape[0]):
        M = Matrix(bb[row, 0].blocks)
        for col in range(1, bb.shape[1]):
            M = M.row_join(bb[row, col].blocks)
        MM = MM.col_join(M)

    return BlockMatrix(MM)


def reblock_2x2(B):
    """Reblock a BlockMatrix so that it has 2x2 blocks of block matrices."""
    if not isinstance(B, BlockMatrix) or not all(d > 2 for d in B.blocks.shape):
        return B

    BM = BlockMatrix  # for brevity's sake
    return BM([[    B.blocks[0,  0], BM(B.blocks[0,  1:])],
               [BM(B.blocks[1:, 0]), BM(B.blocks[1:, 1:])]])


def bounds(sizes):
    """Convert sequence of numbers into pairs of low-high pairs

    >>> bounds((1, 10, 50))
    [(0, 1), (1, 11), (11, 61)]

    """
    low = 0
    rv = []
    for size in sizes:
        rv.append((low, low + size))
        low += size
    return rv


def blockcut(expr, rowsizes, colsizes):
    """Cut a matrix expression into Blocks

    >>> M = ImmutableMatrix(4, 4, range(16))
    >>> B = blockcut(M, (1, 3), (1, 3))
    >>> type(B).__name__
    'BlockMatrix'
    >>> ImmutableMatrix(B.blocks[0, 1])
    Matrix([[1, 2, 3]])

    """
    rowbounds = bounds(rowsizes)
    colbounds = bounds(colsizes)
    return BlockMatrix([[MatrixSlice(expr, rowbound, colbound)
                         for colbound in colbounds]
                        for rowbound in rowbounds])
