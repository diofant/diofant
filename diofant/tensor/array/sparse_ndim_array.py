import functools
import itertools

from ...core import Dict, Expr, Integer, Tuple
from ...core.sympify import sympify
from ...matrices import SparseMatrix
from ...utilities import flatten
from .mutable_ndim_array import MutableNDimArray
from .ndim_array import ImmutableNDimArray, NDimArray


class SparseNDimArray(NDimArray):
    """Sparse N-dim array."""

    def __getitem__(self, index):
        """
        Get an element from a sparse N-dim array.

        Examples
        ========

        >>> a = MutableSparseNDimArray(range(4), (2, 2))
        >>> a
        [[0, 1], [2, 3]]
        >>> a[0, 0]
        0
        >>> a[1, 1]
        3
        >>> a[0]
        0
        >>> a[2]
        2

        Symbolic indexing:
        >>> a[n, m]
        [[0, 1], [2, 3]][n, m]

        Replace `n` and `m` to get element `(0, 0)`:
        >>> a[n, m].subs({n: 0, m: 0})
        0

        """
        syindex = self._check_symbolic_index(index)
        if syindex is not None:
            return syindex

        # `index` is a tuple with one or more slices:
        if isinstance(index, tuple) and any(isinstance(i, slice) for i in index):

            def slice_expand(s, dim):
                if not isinstance(s, slice):
                    return s,
                start, stop, step = s.indices(dim)
                return [start + i*step for i in range((stop-start)//step)]

            sl_factors = [slice_expand(i, dim) for (i, dim) in zip(index, self.shape)]
            eindices = itertools.product(*sl_factors)
            array = [self._sparse_array.get(self._parse_index(i), Integer(0)) for i in eindices]
            nshape = [len(el) for i, el in enumerate(sl_factors) if isinstance(index[i], slice)]
            return type(self)(array, nshape)
        else:
            # `index` is a single slice:
            if isinstance(index, slice):
                start, stop, step = index.indices(self._loop_size)
                retvec = [self._sparse_array.get(ind, Integer(0)) for ind in range(start, stop, step)]
                return retvec
            # `index` is a number or a tuple without any slice:
            else:
                index = self._parse_index(index)
                return self._sparse_array.get(index, Integer(0))

    @classmethod
    def zeros(cls, *shape):
        """Return a sparse N-dim array of zeros."""
        return cls({}, shape)

    def tomatrix(self):
        """
        Converts MutableDenseNDimArray to Matrix. Can convert only 2-dim array, else will raise error.

        Examples
        ========

        >>> a = MutableSparseNDimArray([1 for i in range(9)], (3, 3))
        >>> b = a.tomatrix()
        >>> b
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])

        """
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')

        mat_sparse = {}
        for key, value in self._sparse_array.items():
            mat_sparse[self._get_tuple_index(key)] = value

        return SparseMatrix(self.shape[0], self.shape[1], mat_sparse)

    def __iter__(self):
        def iterator():
            for i in range(self._loop_size):
                yield self[i]
        return iterator()


class ImmutableSparseNDimArray(SparseNDimArray, ImmutableNDimArray):
    """An immutable version of a sparse N-dim array."""

    def __new__(cls, iterable=None, shape=None, **kwargs):

        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple(*(sympify(x, strict=True) for x in shape))
        loop_size = functools.reduce(lambda x, y: x*y, shape) if shape else 0

        # Sparse array:
        if isinstance(flat_list, (dict, Dict)):
            sparse_array = Dict(flat_list)
        else:
            sparse_array = {}
            for i, el in enumerate(flatten(flat_list)):
                if el != 0:
                    sparse_array[i] = sympify(el, strict=True)

        sparse_array = Dict(sparse_array)

        self = Expr.__new__(cls, sparse_array, shape, **kwargs)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size
        self._sparse_array = sparse_array

        return self

    def __setitem__(self, index, value):
        raise TypeError('immutable N-dim array')


class MutableSparseNDimArray(MutableNDimArray, SparseNDimArray):
    """A mutable version of a sparse N-dim array."""

    def __new__(cls, iterable=None, shape=None, **kwargs):

        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        self = object.__new__(cls)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = functools.reduce(lambda x, y: x*y, shape) if shape else 0

        # Sparse array:
        if isinstance(flat_list, (dict, Dict)):
            self._sparse_array = dict(flat_list)
            return self

        self._sparse_array = {}

        for i, el in enumerate(flatten(flat_list)):
            if el != 0:
                self._sparse_array[i] = sympify(el, strict=True)

        return self

    def __setitem__(self, index, value):
        """Allows to set items to MutableDenseNDimArray.

        Examples
        ========

        >>> a = MutableSparseNDimArray.zeros(2, 2)
        >>> a[0, 0] = 1
        >>> a[1, 1] = 1
        >>> a
        [[1, 0], [0, 1]]

        """
        index = self._parse_index(index)
        value = sympify(value, strict=True)

        if value == 0 and index in self._sparse_array:
            self._sparse_array.pop(index)
        else:
            self._sparse_array[index] = value
