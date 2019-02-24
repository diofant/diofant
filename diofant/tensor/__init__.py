"""A module to manipulate symbolic objects with indices including tensors

"""
from .indexed import IndexedBase, Idx, Indexed  # noqa: F401
from .index_methods import get_contraction_structure, get_indices  # noqa: F401
from .array import (MutableDenseNDimArray, ImmutableDenseNDimArray,  # noqa: F401
                    MutableSparseNDimArray, ImmutableSparseNDimArray, NDimArray,
                    tensorproduct, tensorcontraction, derive_by_array, permutedims,
                    Array, DenseNDimArray, SparseNDimArray)
