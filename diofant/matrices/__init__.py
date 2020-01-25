"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .dense import (GramSchmidt, MutableDenseMatrix, MutableMatrix, casoratian,
                    diag, eye, hessian, jordan_cell, list2numpy, matrix2numpy,
                    matrix_multiply_elementwise, ones, randMatrix, rot_axis1,
                    rot_axis2, rot_axis3, symarray, vandermonde, wronskian,
                    zeros)
from .expressions import (Adjoint, BlockDiagMatrix, BlockMatrix, Determinant,
                          DiagonalMatrix, DiagonalOf, FunctionMatrix,
                          HadamardProduct, Identity, Inverse, MatAdd, MatMul,
                          MatPow, MatrixExpr, MatrixSlice, MatrixSymbol, Trace,
                          Transpose, ZeroMatrix, block_collapse, blockcut, det,
                          hadamard_product, trace)
from .immutable import (ImmutableDenseMatrix, ImmutableMatrix,
                        ImmutableSparseMatrix)
from .matrices import MatrixBase, NonSquareMatrixError, ShapeError
from .sparse import MutableSparseMatrix, SparseMatrix


Matrix = MutableMatrix


__all__ = ('GramSchmidt', 'MutableDenseMatrix', 'MutableMatrix', 'casoratian',
           'diag', 'eye', 'hessian', 'jordan_cell', 'list2numpy',
           'matrix2numpy', 'matrix_multiply_elementwise', 'ones', 'randMatrix',
           'rot_axis1', 'rot_axis2', 'rot_axis3', 'symarray', 'vandermonde',
           'wronskian', 'zeros', 'Adjoint', 'BlockDiagMatrix', 'BlockMatrix',
           'Determinant', 'DiagonalMatrix', 'DiagonalOf', 'FunctionMatrix',
           'HadamardProduct', 'Identity', 'Inverse', 'MatAdd', 'MatMul',
           'MatPow', 'MatrixExpr', 'MatrixSlice', 'MatrixSymbol', 'Trace',
           'Transpose', 'ZeroMatrix', 'block_collapse', 'blockcut', 'det',
           'hadamard_product', 'trace', 'ImmutableDenseMatrix',
           'ImmutableMatrix', 'ImmutableSparseMatrix', 'MatrixBase',
           'NonSquareMatrixError', 'ShapeError', 'MutableSparseMatrix',
           'SparseMatrix', 'Matrix')
