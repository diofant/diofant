"""
A Printer for generating executable code.

The most important function here is srepr (that is an exact equivalent of
builtin repr, except for optional arguments) that returns a string so that the
relation eval(srepr(expr))=expr holds in an appropriate environment.
"""

from __future__ import annotations

import typing

import mpmath.libmp as mlib
from mpmath.libmp import prec_to_dps, repr_dps

from ..core.function import AppliedUndef
from ..utilities import default_sort_key
from .defaults import DefaultPrinting
from .printer import Printer


class ReprPrinter(Printer):
    """Repr printer."""

    printmethod = '_diofantrepr'

    _default_settings: dict[str, typing.Any] = {
        'order': None
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
        if hasattr(expr, 'args') and hasattr(expr.args, '__iter__'):
            l = []
            for o in expr.args:
                l.append(self._print(o))
            return expr.__class__.__name__ + '(%s)' % ', '.join(l)
        elif hasattr(expr, '__repr__') and not issubclass(expr.__class__,
                                                          DefaultPrinting):
            return repr(expr)
        else:
            return object.__repr__(expr)

    def _print_Dict(self, expr):
        l = []
        for o in sorted(expr.args, key=default_sort_key):
            l.append(self._print(o))
        return expr.__class__.__name__ + '(%s)' % ', '.join(l)

    def _print_Add(self, expr, order=None):
        args = expr.as_ordered_terms(order=order or self.order)
        args = map(self._print, args)
        return 'Add(%s)' % ', '.join(args)

    def _print_Function(self, expr):
        r = self._print(expr.func)
        r += '(%s)' % ', '.join([self._print(a) for a in expr.args])
        return r

    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef):
            return f'Function({expr.__name__!r})'
        else:
            return expr.__name__

    def _print_RationalConstant(self, expr):
        return f'Rational({expr.numerator}, {expr.denominator})'

    def _print_AtomicExpr(self, expr):
        return str(expr)

    def _print_NumberSymbol(self, expr):
        return str(expr)

    def _print_Integer(self, expr):
        return 'Integer(%i)' % int(expr.numerator)

    def _print_list(self, expr):
        return '[%s]' % self.reprify(expr, ', ')

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
        return f'{expr.__class__.__name__}({self._print(l)})'

    def _print_BooleanTrue(self, expr):
        return 'true'

    def _print_BooleanFalse(self, expr):
        return 'false'

    def _print_NaN(self, expr):
        return 'nan'

    def _print_Mul(self, expr, order=None):
        terms = expr.args
        args = expr._new_rawargs(*terms).as_ordered_factors()
        args = map(self._print, args)
        return 'Mul(%s)' % ', '.join(args)

    def _print_Rational(self, expr):
        return 'Rational(%s, %s)' % (self._print(int(expr.numerator)),
                                     self._print(int(expr.denominator)))

    def _print_Float(self, expr):
        dps = prec_to_dps(expr._prec)
        r = mlib.to_str(expr._mpf_, repr_dps(expr._prec))
        return f"{expr.__class__.__name__}('{r}', dps={dps:d})"

    def _print_BaseSymbol(self, expr):
        d = expr._assumptions.generator
        if d == {}:
            return f'{expr.__class__.__name__}({self._print(expr.name)})'
        else:
            attr = [f'{k}={v}' for k, v in d.items()]
            return '%s(%s, %s)' % (expr.__class__.__name__,
                                   self._print(expr.name), ', '.join(attr))

    def _print_str(self, expr):
        return repr(expr)

    def _print_tuple(self, expr):
        if len(expr) == 1:
            return '(%s,)' % self._print(expr[0])
        else:
            return '(%s)' % self.reprify(expr, ', ')

    def _print_WildFunction(self, expr):
        return f"{expr.__class__.__name__}('{expr.name}')"

    def _print_PolynomialRing(self, ring):
        return '%s(%s, %s, %s)' % (ring.__class__.__name__,
                                   self._print(ring.domain),
                                   self._print(ring.symbols),
                                   self._print(ring.order))

    def _print_GMPYIntegerRing(self, expr):
        return f'{expr.__class__.__name__}()'
    _print_GMPYRationalField = _print_GMPYIntegerRing
    _print_PythonIntegerRing = _print_GMPYIntegerRing
    _print_PythonRationalField = _print_GMPYIntegerRing

    _print_LexOrder = _print_GMPYIntegerRing
    _print_GradedLexOrder = _print_LexOrder

    def _print_FractionField(self, field):
        return '%s(%s, %s, %s)' % (field.__class__.__name__,
                                   self._print(field.domain), self._print(field.symbols), self._print(field.order))

    def _print_PolyElement(self, poly):
        terms = list(poly.items())
        terms.sort(key=poly.ring.order, reverse=True)
        return f'{poly.__class__.__name__}({self._print(poly.ring)}, {self._print(terms)})'

    def _print_FracElement(self, frac):
        numer_terms = list(frac.numerator.items())
        numer_terms.sort(key=frac.field.order, reverse=True)
        denom_terms = list(frac.denominator.items())
        denom_terms.sort(key=frac.field.order, reverse=True)
        numer = self._print(numer_terms)
        denom = self._print(denom_terms)
        return f'{frac.__class__.__name__}({self._print(frac.field)}, {numer}, {denom})'

    def _print_AlgebraicField(self, expr):
        return 'AlgebraicField(%s, %s)' % (self._print(expr.domain),
                                           self._print(expr.ext.as_expr()))

    def _print_AlgebraicElement(self, expr):
        return '%s(%s)' % (self._print(expr.parent),
                           self._print(list(map(expr.domain.domain.to_expr, expr.rep.all_coeffs()))))

    def _print_Domain(self, expr):
        return expr.rep


def srepr(expr, **settings):
    """Return expr in repr form."""
    return ReprPrinter(settings).doprint(expr)
