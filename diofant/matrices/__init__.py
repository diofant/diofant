"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .matrices import (DeferredVector, ShapeError,  # noqa: F401
                       NonSquareMatrixError, MatrixBase)
from .dense import (GramSchmidt, MutableMatrix,  # noqa: F401
                    MutableDenseMatrix, casoratian, diag, eye, hessian,
                    jordan_cell, list2numpy, matrix2numpy,
                    matrix_multiply_elementwise, ones, randMatrix,
                    rot_axis1, rot_axis2, rot_axis3, symarray, wronskian,
                    zeros, vandermonde)
from .sparse import MutableSparseMatrix, SparseMatrix  # noqa: F401
from .immutable import (ImmutableMatrix, ImmutableSparseMatrix,  # noqa: F401
                        ImmutableDenseMatrix)
from .expressions import (MatrixSlice, BlockDiagMatrix,  # noqa: F401
                          BlockMatrix, FunctionMatrix, Identity, Inverse,
                          MatAdd, MatMul, MatPow, MatrixExpr, MatrixSymbol,
                          Trace, Transpose, ZeroMatrix, blockcut,
                          block_collapse, Adjoint, hadamard_product,
                          HadamardProduct, Determinant, det, DiagonalMatrix,
                          DiagonalOf, trace)

Matrix = MutableMatrix
