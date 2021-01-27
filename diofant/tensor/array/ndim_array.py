import collections

from ...core import Expr, Integer
from ...core.sympify import sympify
from ...logic import true
from ...matrices import MatrixBase
from ...printing.defaults import DefaultPrinting
from ..indexed import Indexed


class NDimArray(DefaultPrinting):
    """N-dim array.

    Examples
    ========

    Create an N-dim array of zeros:

    >>> a = MutableDenseNDimArray.zeros(2, 3, 4)
    >>> a
    [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
     [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]

    Create an N-dim array from a list;

    >>> a = MutableDenseNDimArray([[2, 3], [4, 5]])
    >>> a
    [[2, 3], [4, 5]]

    >>> b = MutableDenseNDimArray([[[1, 2], [3, 4], [5, 6]],
    ...                            [[7, 8], [9, 10], [11, 12]]])
    >>> b
    [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]

    Create an N-dim array from a flat list with dimension shape:

    >>> a = MutableDenseNDimArray([1, 2, 3, 4, 5, 6], (2, 3))
    >>> a
    [[1, 2, 3], [4, 5, 6]]

    Create an N-dim array from a matrix:

    >>> a = Matrix([[1, 2], [3, 4]])
    >>> a
    Matrix([
    [1, 2],
    [3, 4]])
    >>> b = MutableDenseNDimArray(a)
    >>> b
    [[1, 2], [3, 4]]

    Arithmetic operations on N-dim arrays

    >>> a = MutableDenseNDimArray([1, 1, 1, 1], (2, 2))
    >>> b = MutableDenseNDimArray([4, 4, 4, 4], (2, 2))
    >>> c = a + b
    >>> c
    [[5, 5], [5, 5]]
    >>> a - b
    [[-3, -3], [-3, -3]]

    """

    def _parse_index(self, index):

        if isinstance(index, (int, Integer)):
            if index >= self._loop_size:
                raise ValueError('index out of range')
            return index

        if len(index) != self._rank:
            raise ValueError('Wrong number of array axes')

        real_index = 0
        # check if input index can exist in current indexing
        for i in range(self._rank):
            if index[i] >= self.shape[i]:
                raise ValueError('Index ' + str(index) + ' out of border')
            real_index = real_index*self.shape[i] + index[i]

        return real_index

    def _get_tuple_index(self, integer_index):
        index = []
        for i, sh in enumerate(reversed(self.shape)):
            index.append(integer_index % sh)
            integer_index //= sh
        index.reverse()
        return tuple(index)

    def _check_symbolic_index(self, index):
        # Check if any index is symbolic:
        tuple_index = (index if isinstance(index, tuple) else (index,))
        if any((isinstance(i, Expr) and (not i.is_number)) for i in tuple_index):
            for i, nth_dim in zip(tuple_index, self.shape):
                i = sympify(i)
                if ((i < 0) is true) or ((i >= nth_dim) is true):
                    raise ValueError('index out of range')
            return Indexed(self, *tuple_index)

    def _setter_iterable_check(self, value):
        if isinstance(value, (collections.abc.Iterable, MatrixBase, NDimArray)):
            raise NotImplementedError

    @classmethod
    def _scan_iterable_shape(cls, iterable):
        def f(pointer):
            if not isinstance(pointer, collections.abc.Iterable):
                return [pointer], ()

            result = []
            elems, shapes = zip(*[f(i) for i in pointer])
            if len(set(shapes)) != 1:
                raise ValueError('could not determine shape unambiguously')
            for i in elems:
                result.extend(i)
            return result, (len(shapes),)+shapes[0]

        return f(iterable)

    @classmethod
    def _handle_ndarray_creation_inputs(cls, iterable=None, shape=None, **kwargs):

        if shape is None and iterable is None:
            shape = ()
            iterable = ()

        # Construction from another `NDimArray`:
        elif shape is None and isinstance(iterable, NDimArray):
            shape = iterable.shape
            iterable = list(iterable)

        # Construct N-dim array from an iterable (numpy arrays included):
        elif shape is None and isinstance(iterable, collections.abc.Iterable):
            iterable, shape = cls._scan_iterable_shape(iterable)

        # Construct N-dim array from a Matrix:
        elif shape is None and isinstance(iterable, MatrixBase):
            shape = iterable.shape

        # Construct NDimArray(iterable, shape)
        elif shape is not None:
            pass

        else:
            raise TypeError('Data type not understood')

        if isinstance(shape, (int, Integer)):
            shape = shape,

        shape = tuple(shape)

        if any(not isinstance(dim, (int, Integer)) for dim in shape):
            raise TypeError('Shape should contain integers only.')

        if isinstance(iterable, collections.abc.Mapping):
            for k, v in list(iterable.items()):
                if not isinstance(k, collections.abc.Sequence):
                    continue
                new_key = 0
                for i, idx in enumerate(k):
                    new_key = new_key * shape[i] + idx
                iterable[new_key] = iterable[k]
                del iterable[k]

        return shape, iterable

    def __len__(self):
        """Overload common function len(). Returns number of elements in array.

        Examples
        ========

        >>> a = MutableDenseNDimArray.zeros(3, 3)
        >>> a
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        >>> len(a)
        9

        """
        return self._loop_size

    @property
    def shape(self):
        """
        Returns array shape (dimension).

        Examples
        ========

        >>> a = MutableDenseNDimArray.zeros(3, 3)
        >>> a.shape
        (3, 3)

        """
        return self._shape

    def rank(self):
        """
        Returns rank of array.

        Examples
        ========

        >>> a = MutableDenseNDimArray.zeros(3, 4, 5, 6, 3)
        >>> a.rank()
        5

        """
        return self._rank

    def diff(self, *args, **kwargs):
        """
        Calculate the derivative of each element in the array.

        Examples
        ========

        >>> M = ImmutableDenseNDimArray([[x, y], [1, x*y]])
        >>> M.diff(x)
        [[1, 0], [0, y]]

        """
        return type(self)(map(lambda x: x.diff(*args, **kwargs), self), self.shape)

    def applyfunc(self, f):
        """Apply a function to each element of the N-dim array.

        Examples
        ========

        >>> m = ImmutableDenseNDimArray([i*2+j for i in range(2)
        ...                              for j in range(2)], (2, 2))
        >>> m
        [[0, 1], [2, 3]]
        >>> m.applyfunc(lambda i: 2*i)
        [[0, 2], [4, 6]]

        """
        return type(self)(map(f, self), self.shape)

    def __str__(self):
        """Returns string, allows to use standard functions print() and str().

        Examples
        ========

        >>> a = MutableDenseNDimArray.zeros(2, 2)
        >>> a
        [[0, 0], [0, 0]]

        """
        def f(sh, shape_left, i, j):
            if len(shape_left) == 1:
                return '['+', '.join([str(self[e]) for e in range(i, j)])+']'

            sh //= shape_left[0]
            return '[' + ', '.join([f(sh, shape_left[1:], i+e*sh, i+(e+1)*sh) for e in range(shape_left[0])]) + ']'  # + '\n'*len(shape_left)

        return f(self._loop_size, self.shape, 0, self._loop_size)

    def tolist(self):
        """
        Conveting MutableDenseNDimArray to one-dim list

        Examples
        ========

        >>> a = MutableDenseNDimArray([1, 2, 3, 4], (2, 2))
        >>> a
        [[1, 2], [3, 4]]
        >>> b = a.tolist()
        >>> b
        [[1, 2], [3, 4]]

        """

        def f(sh, shape_left, i, j):
            if len(shape_left) == 1:
                return [self[e] for e in range(i, j)]
            result = []
            sh //= shape_left[0]
            for e in range(shape_left[0]):
                result.append(f(sh, shape_left[1:], i+e*sh, i+(e+1)*sh))
            return result

        return f(self._loop_size, self.shape, 0, self._loop_size)

    def __add__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError('array shape mismatch')
        result_list = [i + j for i, j in zip(self, other)]

        return type(self)(result_list, self.shape)

    def __sub__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError('array shape mismatch')
        result_list = [i - j for i, j in zip(self, other)]

        return type(self)(result_list, self.shape)

    def __mul__(self, other):
        if isinstance(other, (collections.abc.Iterable, NDimArray, MatrixBase)):
            raise ValueError('scalar expected, use tensorproduct(...) for tensorial product')
        other = sympify(other)
        result_list = [i*other for i in self]
        return type(self)(result_list, self.shape)

    def __rmul__(self, other):
        if isinstance(other, (collections.abc.Iterable, NDimArray, MatrixBase)):
            raise ValueError('scalar expected, use tensorproduct(...) for tensorial product')
        other = sympify(other)
        result_list = [other*i for i in self]
        return type(self)(result_list, self.shape)

    def __truediv__(self, other):
        if isinstance(other, (collections.abc.Iterable, NDimArray, MatrixBase)):
            raise ValueError('scalar expected')
        other = sympify(other)
        result_list = [i/other for i in self]
        return type(self)(result_list, self.shape)

    def __rtruediv__(self, other):
        return NotImplemented

    def __eq__(self, other):
        """
        Compare NDimArray instances.

        Instances equal if they have same shape and data.

        Examples
        ========

        >>> a = MutableDenseNDimArray.zeros(2, 3)
        >>> b = MutableDenseNDimArray.zeros(2, 3)
        >>> a == b
        True
        >>> c = a.reshape(3, 2)
        >>> c == b
        False
        >>> a[0, 0] = 1
        >>> b[0, 0] = 2
        >>> a == b
        False

        """
        if not isinstance(other, NDimArray):
            return False
        return (self.shape == other.shape) and (list(self) == list(other))

    def _eval_transpose(self):
        from .arrayop import permutedims
        if self.rank() != 2:
            raise ValueError('array rank not 2')
        return permutedims(self, (1, 0))

    def transpose(self):
        return self._eval_transpose()

    def _eval_conjugate(self):
        return self.func([i.conjugate() for i in self], self.shape)

    def conjugate(self):
        return self._eval_conjugate()

    def _eval_adjoint(self):
        return self.transpose().conjugate()

    def adjoint(self):
        return self._eval_adjoint()


class ImmutableNDimArray(NDimArray, Expr):
    """An immutable version of the N-dim array."""

    _op_priority = 11.0

    def _subs(self, old, new, **hints):
        return super()._subs(old, new, **hints)
