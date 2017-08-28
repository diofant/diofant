""" A module which handles Matrix Expressions """

from .slice import MatrixSlice  # noqa: F401  # noqa: F401
from .blockmatrix import (BlockMatrix, BlockDiagMatrix,  # noqa: F401
                          block_collapse, blockcut)
from .funcmatrix import FunctionMatrix  # noqa: F401
from .inverse import Inverse  # noqa: F401
from .matadd import MatAdd  # noqa: F401
from .matexpr import (Identity, MatrixExpr,  # noqa: F401
                      MatrixSymbol, ZeroMatrix)
from .matmul import MatMul  # noqa: F401
from .matpow import MatPow  # noqa: F401
from .trace import Trace, trace  # noqa: F401
from .determinant import Determinant, det  # noqa: F401
from .transpose import Transpose  # noqa: F401
from .adjoint import Adjoint  # noqa: F401
from .hadamard import hadamard_product, HadamardProduct  # noqa: F401
from .diagonal import DiagonalMatrix, DiagonalOf  # noqa: F401
