"""
A Printer for generating readable representation of most diofant classes.
"""

from __future__ import annotations

import typing

import mpmath.libmp as mlib
from mpmath.libmp import prec_to_dps

from ..core import Integer, Mul, Pow, Rational, S, oo
from ..core.mul import _keep_coeff
from ..utilities import default_sort_key
from .defaults import DefaultPrinting
from .precedence import PRECEDENCE, precedence
from .printer import Printer


class StrPrinter(Printer):
    """Str printer."""

    printmethod = '_diofantstr'
    _default_settings: dict[str, typing.Any] = {
        'order': None,
        'full_prec': 'auto',
    }

    _relationals: dict[str, str] = {}

    def parenthesize(self, item, level):
        if precedence(item) <= level:
            return '(%s)' % self._print(item)
        else:
            return self._print(item)

    def stringify(self, args, sep, level=0):
        return sep.join([self.parenthesize(item, level) for item in args])

    def emptyPrinter(self, expr):
        if isinstance(expr, str):
            return expr
        elif hasattr(expr, '__str__') and not issubclass(expr.__class__,
                                                         DefaultPrinting):
            return str(expr)
        else:
            return repr(expr)

    def _print_Add(self, expr, order=None):
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = expr.as_ordered_terms(order=order or self.order)

        PREC = precedence(expr)
        l = []
        for term in terms:
            t = self._print(term)
            if t.startswith('-'):
                sign = '-'
                t = t[1:]
            else:
                sign = '+'
            if precedence(term) < PREC:
                l.extend([sign, f'({t})'])
            else:
                l.extend([sign, t])
        sign = l.pop(0)
        if sign == '+':
            sign = ''
        return sign + ' '.join(l)

    def _print_BooleanTrue(self, expr):
        return 'true'

    def _print_BooleanFalse(self, expr):
        return 'false'

    def _print_Not(self, expr):
        return '~%s' % self.parenthesize(expr.args[0], PRECEDENCE['Not'])

    def _print_And(self, expr):
        return self.stringify(expr.args, ' & ', PRECEDENCE['BitwiseAnd'])

    def _print_Or(self, expr):
        return self.stringify(expr.args, ' | ', PRECEDENCE['BitwiseOr'])

    def _print_Basic(self, expr):
        l = [self._print(o) for o in expr.args]
        return expr.__class__.__name__ + '(%s)' % ', '.join(l)

    def _print_BlockMatrix(self, B):
        return self._print(B.blocks)

    def _print_Catalan(self, expr):
        return 'Catalan'

    def _print_ComplexInfinity(self, expr):
        return 'zoo'

    def _print_Derivative(self, expr):
        return 'Derivative(%s)' % ', '.join(map(self._print, expr.args))

    def _print_dict(self, d):
        keys = sorted(d, key=default_sort_key)
        items = []

        for key in keys:
            item = f'{self._print(key)}: {self._print(d[key])}'
            items.append(item)

        return '{%s}' % ', '.join(items)

    def _print_Dict(self, expr):
        return self._print_dict(expr)

    def _print_RandomDomain(self, d):
        return 'Domain: ' + self._print(d.as_boolean())

    def _print_Dummy(self, expr):
        return '_' + expr.name

    def _print_EulerGamma(self, expr):
        return 'EulerGamma'

    def _print_Exp1(self, expr):
        return 'E'

    def _print_ExprCondPair(self, expr):
        return f'({expr.expr}, {expr.cond})'

    def _print_FiniteSet(self, s):
        s = sorted(s, key=default_sort_key)
        if len(s) > 10:
            printset = s[:3] + ['...'] + s[-3:]
        else:
            printset = s
        return '{' + ', '.join(self._print(el) for el in printset) + '}'

    def _print_Function(self, expr):
        return expr.func.__name__ + '(%s)' % self.stringify(expr.args, ', ')

    def _print_GeometryEntity(self, expr):
        # GeometryEntity is special -- it's base is tuple
        return str(expr)

    def _print_GoldenRatio(self, expr):
        return 'GoldenRatio'

    def _print_ImaginaryUnit(self, expr):
        return 'I'

    def _print_Infinity(self, expr):
        return 'oo'

    def _print_Integral(self, expr):
        def _xab_tostr(xab):
            if len(xab) == 1:
                return self._print(xab[0])
            else:
                return self._print((xab[0],) + tuple(xab[1:]))
        L = ', '.join([_xab_tostr(l) for l in expr.limits])
        return f'Integral({self._print(expr.function)}, {L})'

    def _print_Interval(self, i):
        if i.left_open:
            left = '('
        else:
            left = '['

        if i.right_open:
            right = ')'
        else:
            right = ']'

        return f'{left}{self._print(i.start)}, {self._print(i.end)}{right}'

    def _print_Inverse(self, I):
        return '%s^-1' % self.parenthesize(I.arg, PRECEDENCE['Pow'])

    def _print_Lambda(self, obj):
        args, expr = obj.args
        if len(args) == 1:
            return f'Lambda({args.args[0]}, {expr})'
        else:
            arg_string = ', '.join(self._print(arg) for arg in args)
            return f'Lambda(({arg_string}), {expr})'

    def _print_LatticeOp(self, expr):
        args = sorted(expr.args, key=default_sort_key)
        return expr.func.__name__ + '(%s)' % ', '.join(self._print(arg) for arg in args)

    def _print_Limit(self, expr):
        e, z, z0, dir = expr.args
        if str(dir) == '+':
            return f'Limit({e}, {z}, {z0})'
        else:
            return f"Limit({e}, {z}, {z0}, dir='{dir}')"

    def _print_list(self, expr):
        return '[%s]' % self.stringify(expr, ', ')

    def _print_MatrixBase(self, expr):
        return expr._format_str(self)

    def _print_MatrixElement(self, expr):
        return self._print(expr.parent) + f'[{expr.i}, {expr.j}]'

    def _print_MatrixSlice(self, expr):
        def strslice(x):
            x = list(x)
            if x[2] == 1:
                del x[2]
            if x[1] == x[0] + 1:
                del x[1]
            if x[0] == 0:
                x[0] = ''
            return ':'.join(map(self._print, x))
        return (self._print(expr.parent) + '[' +
                strslice(expr.rowslice) + ', ' +
                strslice(expr.colslice) + ']')

    def _print_Mul(self, expr):

        prec = precedence(expr)

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = '-'
        else:
            sign = ''

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order != 'none':
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        multiple_ones = len([x for x in args if x == 1]) > 1

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            elif item.is_Rational and item is not oo:
                if item.numerator != 1 or multiple_ones:
                    a.append(Rational(item.numerator))
                if item.denominator != 1:
                    b.append(Rational(item.denominator))
            else:
                a.append(item)

        a = a or [Integer(1)]

        a_str = [self.parenthesize(x, prec) for x in a]
        b_str = [self.parenthesize(x, prec) for x in b]

        if len(b) == 0:
            return sign + '*'.join(a_str)
        elif len(b) == 1:
            return sign + '*'.join(a_str) + '/' + b_str[0]
        else:
            return sign + '*'.join(a_str) + '/(%s)' % '*'.join(b_str)

    def _print_MatMul(self, expr):
        return '*'.join([self.parenthesize(arg, precedence(expr))
                         for arg in expr.args])

    def _print_HadamardProduct(self, expr):
        return '.*'.join([self.parenthesize(arg, precedence(expr))
                          for arg in expr.args])

    def _print_MatAdd(self, expr):
        return ' + '.join([self.parenthesize(arg, precedence(expr))
                           for arg in expr.args])

    def _print_NaN(self, expr):
        return 'nan'

    def _print_NegativeInfinity(self, expr):
        return '-oo'

    def _print_Order(self, expr):
        if all(p == 0 for p in expr.point) or not len(expr.variables):
            if len(expr.variables) <= 1:
                return 'O(%s)' % self._print(expr.expr)
            else:
                return 'O(%s)' % self.stringify((expr.expr,) + expr.variables, ', ', 0)
        else:
            return 'O(%s)' % self.stringify(expr.args, ', ', 0)

    def _print_Cycle(self, expr):
        """We want it to print as Cycle in doctests for which a repr is required.

        With __repr__ defined in Cycle, interactive output gives Cycle form but
        during doctests, the dict's __repr__ form is used. Defining this _print
        function solves that problem.

        >>> Cycle(1, 2)  # will print as a dict without this method
        Cycle(1, 2)

        """
        return expr.__repr__()

    def _print_Permutation(self, expr):
        from ..combinatorics import Cycle, Permutation
        if Permutation.print_cyclic:
            if not expr.size:
                return 'Permutation()'
            # before taking Cycle notation, see if the last element is
            # a singleton and move it to the head of the string
            s = Cycle(expr)(expr.size - 1).__repr__()[len('Cycle'):]
            last = s.rfind('(')
            if not last == 0 and ',' not in s[last:]:
                s = s[last:] + s[:last]
            return f'Permutation{s}'
        else:
            s = expr.support()
            if not s:
                if expr.size < 5:
                    return 'Permutation(%s)' % str(expr.array_form)
                return f'Permutation([], size={expr.size})'
            trim = str(expr.array_form[:s[-1] + 1]) + f', size={expr.size}'
            use = full = str(expr.array_form)
            if len(trim) < len(full):
                use = trim
            return f'Permutation({use})'

    def _print_TensorIndex(self, expr):
        return expr._print()

    def _print_TensorHead(self, expr):
        return expr._print()

    def _print_Tensor(self, expr):
        return expr._print()

    def _print_TensMul(self, expr):
        return expr._print()

    def _print_TensAdd(self, expr):
        return expr._print()

    def _print_PermutationGroup(self, expr):
        p = ['    %s' % str(a) for a in expr.args]
        return 'PermutationGroup([\n%s])' % ',\n'.join(p)

    def _print_Pi(self, expr):
        return 'pi'

    def _print_PolyElement(self, poly):
        return poly._str(self, PRECEDENCE, '%s**%d', '*')

    def _print_FracElement(self, frac):
        if frac.denominator == 1:
            return self._print(frac.numerator)
        else:
            numer = self.parenthesize(frac.numerator, PRECEDENCE['Add'])
            denom = self.parenthesize(frac.denominator, PRECEDENCE['Atom']-1)
            return numer + '/' + denom

    def _print_Poly(self, expr):
        terms, gens = [], expr.gens

        for monom, coeff in expr.terms():
            s_monom = []

            for i, exp in enumerate(monom):
                if exp > 0:
                    if exp == 1:
                        s_monom.append(self._print(gens[i]))
                    else:
                        s_monom.append(self.parenthesize(gens[i],
                                                         PRECEDENCE['Atom'] - 1) + f'**{exp:d}')

            s_monom = '*'.join(s_monom)

            if coeff.is_Add:
                if s_monom:
                    s_coeff = '(' + self._print(coeff) + ')'
                else:
                    s_coeff = self._print(coeff)
            else:
                if s_monom:
                    if coeff == 1:
                        terms.extend(['+', s_monom])
                        continue

                    if coeff == -1:
                        terms.extend(['-', s_monom])
                        continue

                s_coeff = self._print(coeff)

            if not s_monom:
                s_term = s_coeff
            else:
                s_term = s_coeff + '*' + s_monom

            if s_term.startswith('-'):
                terms.extend(['-', s_term[1:]])
            else:
                terms.extend(['+', s_term])

        if not terms:
            terms.extend(['+', '0'])

        modifier = terms.pop(0)

        if modifier == '-':
            terms[0] = '-' + terms[0]

        format = expr.__class__.__name__ + '(%s, %s'

        from ..polys.polyerrors import PolynomialError

        try:
            format += ', modulus=%s' % expr.get_modulus()
        except PolynomialError:
            format += f", domain='{expr.domain}'"

        format += ')'

        return format % (' '.join(terms), ', '.join(self._print(s) for s in expr.gens))

    def _print_ProductSet(self, p):
        return ' x '.join(self._print(set) for set in p.sets)

    def _print_AlgebraicElement(self, expr):
        return self._print(expr.parent.to_expr(expr))

    def _print_ModularInteger(self, expr):
        return f'{expr.rep} mod {expr.parent.characteristic}'

    def _print_GaloisFieldElement(self, expr):
        from ..domains import ZZ_python
        return 'GF(%s, %s)(%s)' % (expr.parent.characteristic,
                                   expr.mod.set_domain(ZZ_python).all_coeffs(),
                                   int(expr))

    def _print_Pow(self, expr, rational=False):
        PREC = precedence(expr)

        if not expr.exp.is_Float and expr.exp == Rational(1, 2) and not rational:
            return 'sqrt(%s)' % self._print(expr.base)

        if expr.is_commutative and not expr.exp.is_Float:
            if -expr.exp == Rational(1, 2) and not rational:
                return '1/sqrt(%s)' % self._print(expr.base)
            if expr.exp == -1:
                return '1/%s' % self.parenthesize(expr.base, PREC)

        e = self.parenthesize(expr.exp, PREC)
        if (self.printmethod == '_diofantrepr' and
                expr.exp.is_Rational and not expr.exp.is_Integer):
            # The parenthesized exp should be '(Rational(a, b))' so strip
            # parens, but just check to be sure.
            return f'{self.parenthesize(expr.base, PREC)}**{e[1:-1]}'
        return f'{self.parenthesize(expr.base, PREC)}**{e}'

    def _print_Mod(self, expr):
        PREC = precedence(expr)

        a, b = expr.args
        return f'{self.parenthesize(a, PREC)}%{self.parenthesize(b, PREC)}'

    def _print_MatPow(self, expr):
        PREC = precedence(expr)
        return '%s**%s' % (self.parenthesize(expr.base, PREC),
                           self.parenthesize(expr.exp, PREC))

    def _print_ImmutableDenseNDimArray(self, expr):
        return str(expr)

    _print_ImmutableSparseNDimArray = _print_ImmutableDenseNDimArray

    _print_MutableDenseNDimArray = _print_ImmutableDenseNDimArray
    _print_MutableSparseNDimArray = _print_ImmutableDenseNDimArray

    def _print_Integer(self, expr):
        return str(expr.numerator)

    def _print_Rational(self, expr):
        if expr.denominator == 1:
            return str(expr.numerator)
        else:
            return f'{expr.numerator}/{expr.denominator}'

    def _print_Float(self, expr):
        prec = expr._prec
        if prec < 5:
            dps = 0
        else:
            dps = prec_to_dps(expr._prec)
        if self._settings['full_prec'] is True:
            strip = False
        elif self._settings['full_prec'] is False:
            strip = True
        elif self._settings['full_prec'] == 'auto':
            strip = self._print_level > 1
        else:
            raise NotImplementedError
        rv = mlib.to_str(expr._mpf_, dps, strip_zeros=strip)
        if rv.startswith('-.0'):
            rv = '-0.' + rv[3:]
        elif rv.startswith('.0'):
            rv = '0.' + rv[2:]
        elif rv.startswith('+'):
            rv = rv[1:]
        return rv

    def _print_Relational(self, expr):

        charmap = {
            '==': 'Eq',
            '!=': 'Ne',
        }

        if expr.rel_op in charmap:
            return f'{charmap[expr.rel_op]}({expr.lhs}, {expr.rhs})'

        return '%s %s %s' % (self.parenthesize(expr.lhs, precedence(expr)),
                             self._relationals.get(expr.rel_op) or expr.rel_op,
                             self.parenthesize(expr.rhs, precedence(expr)))

    def _print_RootOf(self, expr):
        p = self._print_Add(expr.expr, order='lex')
        if expr.free_symbols:
            return f'RootOf({p}, {expr.poly.gen}, {expr.index:d})'
        else:
            return f'RootOf({p}, {expr.index:d})'

    def _print_RootSum(self, expr):
        args = [self._print_Add(expr.expr, order='lex')]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun))

        return 'RootSum(%s)' % ', '.join(args)

    def _print_GroebnerBasis(self, basis):
        cls = basis.__class__.__name__

        exprs = [ self._print_Add(arg, order=basis.order)
                  for arg in basis.exprs ]
        exprs = '[%s]' % ', '.join(exprs)

        gens = [ self._print(gen) for gen in basis.gens ]
        domain = "domain='%s'" % self._print(basis.domain)
        order = "order='%s'" % self._print(basis.order)

        args = [exprs] + gens + [domain, order]

        return '%s(%s)' % (cls, ', '.join(args))

    def _print_frozenset(self, s):
        items = sorted(s, key=default_sort_key)

        args = ', '.join(self._print(item) for item in items)
        if args:
            args = '{%s}' % args
        return f'{type(s).__name__}({args})'

    def _print_set(self, expr):
        if expr == set():
            return 'set()'
        return '{%s}' % self.stringify(sorted(expr, key=default_sort_key), ', ')

    def _print_Sum(self, expr):
        def _xab_tostr(xab):
            return self._print((xab[0],) + tuple(xab[1:]))
        L = ', '.join([_xab_tostr(l) for l in expr.limits])
        return f'Sum({self._print(expr.function)}, {L})'

    def _print_Symbol(self, expr):
        return expr.name
    _print_BaseSymbol = _print_Symbol
    _print_MatrixSymbol = _print_Symbol
    _print_RandomSymbol = _print_Symbol

    def _print_Identity(self, expr):
        return 'I'

    def _print_ZeroMatrix(self, expr):
        return '0'

    def _print_str(self, expr):
        return expr

    def _print_tuple(self, expr):
        if len(expr) == 1:
            return '(%s,)' % self._print(expr[0])
        else:
            return '(%s)' % self.stringify(expr, ', ')

    def _print_Tuple(self, expr):
        return self._print_tuple(expr)

    def _print_Monomial(self, expr):
        if expr.gens:
            return '*'.join(['%s**%s' % (gen, exp)
                             for gen, exp in zip(expr.gens, expr)])
        else:
            return self._print_tuple(expr)

    def _print_Transpose(self, T):
        return '%s.T' % self.parenthesize(T.arg, PRECEDENCE['Pow'])

    def _print_Union(self, expr):
        return ' U '.join(self._print(set) for set in expr.args)

    def _print_Complement(self, expr):
        return r' \ '.join(self._print(set) for set in expr.args)

    def _print_Wild(self, expr):
        return expr.name + '_'

    def _print_WildFunction(self, expr):
        return expr.name + '_'

    def _print_Zero(self, expr):
        return '0'

    def _print_BaseScalarField(self, field):
        return field._coord_sys._names[field._index]

    def _print_BaseVectorField(self, field):
        return f'e_{field._coord_sys._names[field._index]}'

    def _print_Differential(self, diff):
        field = diff._form_field
        if hasattr(field, '_coord_sys'):
            return f'd{field._coord_sys._names[field._index]}'
        else:
            return 'd(%s)' % self._print(field)

    def _print_Tr(self, expr):
        # TODO : Handle indices
        return '%s(%s)' % ('Tr', self._print(expr.args[0]))

    def _print_Domain(self, expr):
        return expr.rep


def sstr(expr, **settings):
    """Returns the expression as a string.

    For large expressions where speed is a concern, use the setting
    order='none'.

    Examples
    ========

    >>> sstr(Eq(a + b, 0))
    'Eq(a + b, 0)'

    """
    p = StrPrinter(settings)
    s = p.doprint(expr)

    return s


class StrReprPrinter(StrPrinter):
    """(internal) -- see sstrrepr"""

    def _print_str(self, s):
        return repr(s)


def sstrrepr(expr, **settings):
    """Return expr in mixed str/repr form.

    i.e. strings are returned in repr form with quotes, and everything else
    is returned in str form.

    This function could be useful for hooking into sys.displayhook
    """
    p = StrReprPrinter(settings)
    s = p.doprint(expr)

    return s
