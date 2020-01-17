"""A module which handles Matrix Expressions."""

from .adjoint import Adjoint
from .blockmatrix import BlockDiagMatrix, BlockMatrix, block_collapse, blockcut
from .determinant import Determinant, det
from .diagonal import DiagonalMatrix, DiagonalOf
from .funcmatrix import FunctionMatrix
from .hadamard import HadamardProduct, hadamard_product
from .inverse import Inverse
from .matadd import MatAdd
from .matexpr import Identity, MatrixExpr, MatrixSymbol, ZeroMatrix
from .matmul import MatMul
from .matpow import MatPow
from .slice import MatrixSlice
from .trace import Trace, trace
from .transpose import Transpose


__all__ = ('Adjoint', 'BlockDiagMatrix', 'BlockMatrix', 'block_collapse',
           'blockcut', 'Determinant', 'det', 'DiagonalMatrix', 'DiagonalOf',
           'FunctionMatrix', 'HadamardProduct', 'hadamard_product', 'Inverse',
           'MatAdd', 'Identity', 'MatrixExpr', 'MatrixSymbol', 'ZeroMatrix',
           'MatMul', 'MatPow', 'MatrixSlice', 'Trace', 'trace', 'Transpose')
