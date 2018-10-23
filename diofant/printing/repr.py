"""
A Printer for generating executable code.

The most important function here is srepr (that is an exact equivalent of
builtin repr, except for optional arguments) that returns a string so that the
relation eval(srepr(expr))=expr holds in an appropriate environment.
"""

import mpmath.libmp as mlib
from mpmath.libmp import prec_to_dps, repr_dps

from ..core.compatibility import default_sort_key
from ..core.function import AppliedUndef
from .printer import Printer


class ReprPrinter(Printer):
    printmethod = "_diofantrepr"

    _default_settings = {
        "order": None
    }

    def reprify(self, args, sep):
        """
        Prints each item in `args` and joins them with `sep`.
        """
        return sep.join([self.doprint(item) for item in args])

    def emptyPrinter(self, expr):
        """
        The fallback printer.
        """
        if hasattr(expr, "args") and hasattr(expr.args, "__iter__"):
            l = []
            for o in expr.args:
                l.append(self._print(o))
            return expr.__class__.__name__ + '(%s)' % ', '.join(l)
        else:
            return repr(expr)

    def _print_Dict(self, expr):
        l = []
        for o in sorted(expr.args, key=default_sort_key):
            l.append(self._print(o))
        return expr.__class__.__name__ + '(%s)' % ', '.join(l)

    def _print_Add(self, expr, order=None):
        args = expr.as_ordered_terms(order=order or self.order)
        args = map(self._print, args)
        return "Add(%s)" % ", ".join(args)

    def _print_Function(self, expr):
        r = self._print(expr.func)
        r += '(%s)' % ', '.join([self._print(a) for a in expr.args])
        return r

    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef):
            return 'Function(%r)' % (expr.__name__)
        else:
            return expr.__name__

    def _print_RationalConstant(self, expr):
        return 'Rational(%s, %s)' % (expr.numerator, expr.denominator)

    def _print_AtomicExpr(self, expr):
        return str(expr)

    def _print_NumberSymbol(self, expr):
        return str(expr)

    def _print_Integer(self, expr):
        return 'Integer(%i)' % int(expr.numerator)

    def _print_list(self, expr):
        return "[%s]" % self.reprify(expr, ", ")

    def _print_MatrixBase(self, expr):
        # special case for some empty matrices
        if (expr.rows == 0) ^ (expr.cols == 0):
            return '%s(%s, %s, %s)' % (expr.__class__.__name__,
                                       self._print(expr.rows),
                                       self._print(expr.cols),
                                       self._print([]))
        l = []
        for i in range(expr.rows):
            l.append([])
            for j in range(expr.cols):
                l[-1].append(expr[i, j])
        return '%s(%s)' % (expr.__class__.__name__, self._print(l))

    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_BooleanTrue(self, expr):
        return "true"

    def _print_BooleanFalse(self, expr):
        return "false"

    def _print_NaN(self, expr):
        return "nan"

    def _print_Mul(self, expr, order=None):
        terms = expr.args
        args = expr._new_rawargs(*terms).as_ordered_factors()
        args = map(self._print, args)
        return "Mul(%s)" % ", ".join(args)

    def _print_Rational(self, expr):
        return 'Rational(%s, %s)' % (self._print(int(expr.numerator)),
                                     self._print(int(expr.denominator)))

    def _print_Float(self, expr):
        dps = prec_to_dps(expr._prec)
        r = mlib.to_str(expr._mpf_, repr_dps(expr._prec))
        return "%s('%s', dps=%i)" % (expr.__class__.__name__, r, dps)

    def _print_Symbol(self, expr):
        d = expr._assumptions.generator
        if d == {}:
            return "%s(%s)" % (expr.__class__.__name__, self._print(expr.name))
        else:
            attr = ['%s=%s' % (k, v) for k, v in d.items()]
            return "%s(%s, %s)" % (expr.__class__.__name__,
                                   self._print(expr.name), ', '.join(attr))
    _print_Dummy = _print_Symbol
    _print_Wild = _print_Symbol

    def _print_str(self, expr):
        return repr(expr)

    def _print_tuple(self, expr):
        if len(expr) == 1:
            return "(%s,)" % self._print(expr[0])
        else:
            return "(%s)" % self.reprify(expr, ", ")

    def _print_WildFunction(self, expr):
        return "%s('%s')" % (expr.__class__.__name__, expr.name)

    def _print_PolynomialRing(self, ring):
        return "%s(%s, %s, %s)" % (ring.__class__.__name__,
                                   self._print(ring.domain), self._print(ring.symbols), self._print(ring.order))

    def _print_GMPYIntegerRing(self, expr):
        return "%s()" % expr.__class__.__name__
    _print_GMPYRationalField = _print_GMPYIntegerRing
    _print_PythonIntegerRing = _print_GMPYIntegerRing
    _print_PythonRationalField = _print_GMPYIntegerRing

    _print_LexOrder = _print_GMPYIntegerRing
    _print_GradedLexOrder = _print_LexOrder

    def _print_FractionField(self, field):
        return "%s(%s, %s, %s)" % (field.__class__.__name__,
                                   self._print(field.domain), self._print(field.symbols), self._print(field.order))

    def _print_PolyElement(self, poly):
        terms = list(poly.terms())
        terms.sort(key=poly.ring.order, reverse=True)
        return "%s(%s, %s)" % (poly.__class__.__name__, self._print(poly.ring), self._print(terms))

    def _print_FracElement(self, frac):
        numer_terms = list(frac.numer.terms())
        numer_terms.sort(key=frac.field.order, reverse=True)
        denom_terms = list(frac.denom.terms())
        denom_terms.sort(key=frac.field.order, reverse=True)
        numer = self._print(numer_terms)
        denom = self._print(denom_terms)
        return "%s(%s, %s, %s)" % (frac.__class__.__name__, self._print(frac.field), numer, denom)

    def _print_AlgebraicField(self, expr):
        return "AlgebraicField(%s, %s)" % (self._print(expr.domain),
                                           self._print(expr.ext.as_expr()))

    def _print_AlgebraicElement(self, expr):
        return "%s(%s)" % (self._print(expr.parent),
                           self._print(list(map(expr.domain.to_expr, expr.rep))))


def srepr(expr, **settings):
    """return expr in repr form"""
    return ReprPrinter(settings).doprint(expr)
