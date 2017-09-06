import inspect
from functools import partial

from ..external import import_module
from .printer import Printer


theano = import_module('theano')
if theano:  # pragma: no branch
    ts = theano.scalar
    tt = theano.tensor
    from theano.sandbox import linalg as tlinalg

    from .. import (Add, Mul, Abs, sign, ceiling, floor, log, exp, sqrt,
                    cos, acos, sin, asin, tan, atan, atan2, cosh, acosh,
                    sinh, asinh, tanh, atanh, re, im, arg, erf, gamma,
                    loggamma, Pow, Eq, StrictGreaterThan, StrictLessThan,
                    LessThan, GreaterThan, And, Or, Max, Min, MatAdd,
                    HadamardProduct, Trace, Determinant, Inverse, Transpose,
                    Basic, Piecewise)

    mapping = {
        Add: tt.add,
        Mul: tt.mul,
        Abs: tt.abs_,
        sign: tt.sgn,
        ceiling: tt.ceil,
        floor: tt.floor,
        log: tt.log,
        exp: tt.exp,
        sqrt: tt.sqrt,
        cos: tt.cos,
        acos: tt.arccos,
        sin: tt.sin,
        asin: tt.arcsin,
        tan: tt.tan,
        atan: tt.arctan,
        atan2: tt.arctan2,
        cosh: tt.cosh,
        acosh: tt.arccosh,
        sinh: tt.sinh,
        asinh: tt.arcsinh,
        tanh: tt.tanh,
        atanh: tt.arctanh,
        re: tt.real,
        im: tt.imag,
        arg: tt.angle,
        erf: tt.erf,
        gamma: tt.gamma,
        loggamma: tt.gammaln,
        Pow: tt.pow,
        Eq: tt.eq,
        StrictGreaterThan: tt.gt,
        StrictLessThan: tt.lt,
        LessThan: tt.le,
        GreaterThan: tt.ge,
        And: tt.and_,
        Or: tt.or_,
        Max: tt.maximum,  # Diofant accept >2 inputs, Theano only 2
        Min: tt.minimum,  # Diofant accept >2 inputs, Theano only 2

        # Matrices
        MatAdd: tt.Elemwise(ts.add),
        HadamardProduct: tt.Elemwise(ts.mul),
        Trace: tlinalg.trace,
        Determinant: tlinalg.det,
        Inverse: tlinalg.matrix_inverse,
        Transpose: tt.DimShuffle((False, False), [1, 0]),
    }


class TheanoPrinter(Printer):
    """ Code printer for Theano computations """

    printmethod = "_theano"

    def __init__(self, *args, **kwargs):
        self.cache = kwargs.pop('cache', {})
        super(TheanoPrinter, self).__init__(*args, **kwargs)

    def _print_Symbol(self, s, dtypes={}, broadcastables={}):
        dtype = dtypes.get(s, 'floatX')
        broadcastable = broadcastables.get(s, ())
        key = (s.name, dtype, broadcastable, type(s))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.tensor(name=s.name, dtype=dtype, broadcastable=broadcastable)
            self.cache[key] = value
            return value

    def _print_AppliedUndef(self, s, dtypes={}, broadcastables={}):
        dtype = dtypes.get(s, 'floatX')
        broadcastable = broadcastables.get(s, ())
        name = str(type(s)) + '_' + str(s.args[0])
        key = (name, dtype, broadcastable, type(s), s.args)
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.tensor(name=name, dtype=dtype, broadcastable=broadcastable)
            self.cache[key] = value
            return value

    def _print_Basic(self, expr, **kwargs):
        op = mapping[type(expr)]
        children = [self._print(arg, **kwargs) for arg in expr.args]
        return op(*children)

    def _print_MatrixSymbol(self, X, dtypes={}, **kwargs):
        dtype = dtypes.get(X, 'floatX')
        key = (X.name, dtype, type(X))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.Tensor(dtype, (False, False))(X.name)
            self.cache[key] = value
            return value

    def _print_DenseMatrix(self, X, **kwargs):
        return tt.stacklists([[self._print(arg, **kwargs) for arg in L]
                              for L in X.tolist()])
    _print_ImmutableMatrix = _print_DenseMatrix

    def _print_MatMul(self, expr, **kwargs):
        children = [self._print(arg, **kwargs) for arg in expr.args]
        result = children[0]
        for child in children[1:]:
            result = tt.dot(result, child)
        return result

    def _print_MatrixSlice(self, expr, **kwargs):
        parent = self._print(expr.parent, **kwargs)
        rowslice = self._print(slice(*expr.rowslice), **kwargs)
        colslice = self._print(slice(*expr.colslice), **kwargs)
        return parent[rowslice, colslice]

    def _print_BlockMatrix(self, expr, **kwargs):
        nrows, ncols = expr.blocks.shape
        blocks = [[self._print(expr.blocks[r, c], **kwargs)
                   for c in range(ncols)]
                  for r in range(nrows)]
        return tt.join(0, *[tt.join(1, *row) for row in blocks])

    def _print_slice(self, expr, **kwargs):
        return slice(*[self._print(i, **kwargs)
                       if isinstance(i, Basic) else i
                       for i in (expr.start, expr.stop, expr.step)])

    def _print_Exp1(self, expr, **kwargs):
        return 2.718281828459045

    def _print_Piecewise(self, expr, **kwargs):
        import numpy as np
        e, cond = expr.args[0].args
        if len(expr.args) == 1:
            return tt.switch(self._print(cond, **kwargs),
                             self._print(e, **kwargs),
                             np.nan)
        return tt.switch(self._print(cond, **kwargs),
                         self._print(e, **kwargs),
                         self._print(Piecewise(*expr.args[1:]), **kwargs))

    def _print_Rational(self, expr, **kwargs):
        return tt.true_div(self._print(expr.p, **kwargs),
                           self._print(expr.q, **kwargs))

    def _print_Integer(self, expr, **kwargs):
        return expr.p

    def _print_factorial(self, expr, **kwargs):
        return self._print(gamma(expr.args[0] + 1), **kwargs)

    def _print_Derivative(self, deriv, **kwargs):
        rv = self._print(deriv.expr, **kwargs)
        for var in deriv.variables:
            var = self._print(var, **kwargs)
            rv = tt.Rop(rv, var, tt.ones_like(var))
        return rv

    def emptyPrinter(self, expr):
        return expr

    def doprint(self, expr, **kwargs):
        """Returns printer's representation for expr (as a string)"""
        return self._print(expr, **kwargs)


global_cache = {}


def theano_code(expr, cache=global_cache, **kwargs):
    if not theano:  # pragma: no cover
        raise ImportError("theano is required for theano_code")
    return TheanoPrinter(cache=cache, settings={}).doprint(expr, **kwargs)


def dim_handling(inputs, dim=None, dims={}, broadcastables={}, keys=(),
                 **kwargs):
    """ Handle various input types for dimensions in tensor_wrap

    See Also:
        tensor_wrap
        theano_funciton
    """
    if dim:
        dims = dict(zip(inputs, [dim]*len(inputs)))
    if dims:
        maxdim = max(dims.values())
        broadcastables = {i: (False,)*dims[i] + (True,)*(maxdim-dims[i])
                          for i in inputs}
    return broadcastables


def theano_function(inputs, outputs, dtypes={}, cache=None, **kwargs):
    """ Create Theano function from Diofant expressions """
    if not theano:  # pragma: no cover
        raise ImportError("theano is required for theano_function")
    cache = {} if cache is None else cache
    broadcastables = dim_handling(inputs, **kwargs)

    # Remove keyword arguments corresponding to dim_handling
    dim_names = inspect.getfullargspec(dim_handling)[0]
    theano_kwargs = {k: v for k, v in kwargs.items() if k not in dim_names}

    code = partial(theano_code, cache=cache, dtypes=dtypes,
                   broadcastables=broadcastables)
    tinputs = list(map(code, inputs))
    toutputs = list(map(code, outputs))
    toutputs = toutputs[0] if len(toutputs) == 1 else toutputs
    return theano.function(tinputs, toutputs, **theano_kwargs)
