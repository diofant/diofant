import functools

from ..utilities import default_sort_key
from .str import StrPrinter


class LambdaPrinter(StrPrinter):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

    def _print_MatrixBase(self, expr):
        return "%s(%s)" % (expr.__class__.__name__,
                           self._print((expr.tolist())))

    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_Piecewise(self, expr):
        result = []
        i = 0
        for arg in expr.args:
            e = arg.expr
            c = arg.cond
            result.append('((')
            result.append(self._print(e))
            result.append(') if (')
            result.append(self._print(c))
            result.append(') else (')
            i += 1
        result = result[:-1]
        result.append(') else None)')
        result.append(')'*(2*i - 2))
        return ''.join(result)

    def _print_And(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' and ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Or(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' or ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Not(self, expr):
        result = ['(', 'not (', self._print(expr.args[0]), '))']
        return ''.join(result)

    def _print_BooleanTrue(self, expr):
        return "True"

    def _print_BooleanFalse(self, expr):
        return "False"

    def _print_ITE(self, expr):
        result = [
            '((', self._print(expr.args[1]),
            ') if (', self._print(expr.args[0]),
            ') else (', self._print(expr.args[2]), '))'
        ]
        return ''.join(result)

    def _print_Dummy(self, expr):
        return super()._print_Dummy(expr).replace("(", "_lpar_").replace(")", "_rpar_")


class NumPyPrinter(LambdaPrinter):
    """
    Numpy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """

    _default_settings = {
        "order": "none",
        "full_prec": "auto",
    }

    def _print_MatMul(self, expr):
        """Matrix multiplication printer"""
        return '({0})'.format(').dot('.join(self._print(i) for i in expr.args))

    def _print_Piecewise(self, expr):
        """Piecewise function printer"""
        exprs = '[{0}]'.format(','.join(self._print(arg.expr) for arg in expr.args))
        conds = '[{0}]'.format(','.join(self._print(arg.cond) for arg in expr.args))
        # If [default_value, True] is a (expr, cond) sequence in a Piecewise object
        #     it will behave the same as passing the 'default' kwarg to select()
        #     *as long as* it is the last element in expr.args.
        # If this is not the case, it may be triggered prematurely.
        return 'select({0}, {1}, default=nan)'.format(conds, exprs)

    def _print_Relational(self, expr):
        """Relational printer"""
        op = {'==': 'equal',
              '!=': 'not_equal',
              '<': 'less',
              '<=': 'less_equal',
              '>': 'greater',
              '>=': 'greater_equal'}
        return '{op}({lhs}, {rhs})'.format(op=op[expr.rel_op],
                                           lhs=self._print(expr.lhs),
                                           rhs=self._print(expr.rhs))

    def _print_And(self, expr):
        """Logical And printer"""
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to NUMPY_TRANSLATIONS.
        return functools.reduce(lambda x, y: 'logical_and({0}, {1})'.format(self._print(x), self._print(y)), expr.args)

    def _print_Or(self, expr):
        """Logical Or printer"""
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to NUMPY_TRANSLATIONS.
        return functools.reduce(lambda x, y: 'logical_or({0}, {1})'.format(self._print(x), self._print(y)), expr.args)

    def _print_Xor(self, expr):
        """Logical Xor printer"""
        return functools.reduce(lambda x, y: 'logical_xor({0}, {1})'.format(self._print(x), self._print(y)), expr.args)

    def _print_Not(self, expr):
        """Logical Not printer"""
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{0}({1})'.format('logical_not', ','.join(self._print(i) for i in expr.args))

    def _print_Min(self, expr):
        return '{0}(({1}))'.format('amin', ','.join(self._print(i) for i in expr.args))

    def _print_Max(self, expr):
        return '{0}(({1}))'.format('amax', ','.join(self._print(i) for i in expr.args))


class MpmathPrinter(LambdaPrinter):
    """Mpmath printer."""

    def _print_RootOf(self, expr):
        if expr.is_real:
            return ("findroot(lambda %s: %s, %s, "
                    "method='bisection')" % (self._print(expr.poly.gen),
                                             self._print(expr.expr),
                                             self._print(expr.interval.as_tuple())))
        else:
            return ("findroot(lambda %s: %s, mpc%s, "
                    "method='secant')" % (self._print(expr.poly.gen),
                                          self._print(expr.expr),
                                          self._print(expr.interval.center)))

    def _print_Sum(self, expr):
        return "nsum(lambda %s: %s, %s)" % (",".join([self._print(v) for v in expr.variables]),
                                            self._print(expr.function),
                                            ",".join([self._print(v[1:]) for v in expr.limits]))

    def _print_Infinity(self, expr):
        return "inf"

    def _print_Float(self, e):
        # XXX: This does not handle setting mpmath.mp.dps. It is assumed that
        # the caller of the lambdified function will have set it to sufficient
        # precision to match the Floats in the expression.

        # Remove 'mpz' if gmpy is installed.
        args = str(tuple(map(int, e._mpf_)))
        return 'mpf(%s)' % args

    def _print_GoldenRatio(self, expr):
        return "phi"

    def _print_Pow(self, expr):
        if expr.exp.is_Rational:
            n, d = expr.exp.as_numer_denom()
            if d == 1:
                if n >= 0:
                    return "%s**%s" % (self._print(expr.base), n)
                else:
                    return "power(%s, %s)" % (self._print(expr.base), n)
            else:
                if n >= 2:
                    return "root(%s, %s)**%s" % (self._print(expr.base), d, n)
                elif n == 1:
                    return "root(%s, %s)" % (self._print(expr.base), d)
                else:
                    return "power(root(%s, %s), %s)" % (self._print(expr.base),
                                                        d, n)
        else:
            return super()._print_Pow(expr)

    def _print_Rational(self, expr):
        n, d = expr.numerator, expr.denominator
        if d == 1:
            return "%s" % n
        else:
            return "%s*power(%s, -1)" % (n, d)


def lambdarepr(expr, **settings):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter(settings).doprint(expr)
