import itertools

from ...core import Add, Equality, Integer, Mul, Pow, Rational, S, Symbol, oo
from ...core.function import _coeff_isneg
from ...logic import true
from ...utilities import default_sort_key, group
from ..conventions import requires_partial
from ..printer import Printer
from ..str import sstr
from .pretty_symbology import (annotated, greek_unicode, hobj, pretty_atom,
                               pretty_symbol, pretty_use_unicode, vobj, xobj,
                               xsym)
from .stringpict import prettyForm, stringPict


# rename for usage from outside
pprint_use_unicode = pretty_use_unicode


class PrettyPrinter(Printer):
    """Printer, which converts an expression into 2D ASCII-art figure."""

    printmethod = '_pretty'

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'use_unicode': None,
        'wrap_line': True,
        'num_columns': None,
    }

    def __init__(self, settings=None):
        Printer.__init__(self, settings)
        self.emptyPrinter = lambda x: prettyForm(str(x))

    @property
    def _use_unicode(self):
        if self._settings['use_unicode']:
            return True
        else:
            return pretty_use_unicode()

    def doprint(self, expr):
        return self._print(expr).render(**self._settings)

    # empty op so _print(stringPict) returns the same
    def _print_stringPict(self, e):
        return e

    def _print_atan2(self, e):
        pform = prettyForm(*self._print_seq(e.args).parens())
        pform = prettyForm(*pform.left('atan2'))
        return pform

    def _print_Symbol(self, e):
        symb = pretty_symbol(e.name)
        return prettyForm(symb)
    _print_Dummy = _print_Symbol
    _print_RandomSymbol = _print_Symbol

    def _print_Float(self, e):
        # we will use StrPrinter's Float printer, but we need to handle the
        # full_prec ourselves, according to the self._print_level
        full_prec = self._settings['full_prec']
        if full_prec == 'auto':
            full_prec = self._print_level == 1
        return prettyForm(sstr(e, full_prec=full_prec))

    def _print_Atom(self, e):
        try:
            # print atoms like Exp1 or Pi
            return prettyForm(pretty_atom(e.__class__.__name__))
        except KeyError:
            return self.emptyPrinter(e)

    # Infinity inherits from Number, so we have to override _print_XXX order
    _print_Infinity = _print_Atom
    _print_NegativeInfinity = _print_Atom
    _print_EmptySet = _print_Atom
    _print_Naturals = _print_Atom
    _print_Naturals0 = _print_Atom
    _print_Integers = _print_Atom
    _print_Rationals = _print_Atom
    _print_Reals = _print_Atom

    def _print_subfactorial(self, e):
        x = e.args[0]
        pform = self._print(x)
        # Add parentheses if needed
        if not ((x.is_Integer and x.is_nonnegative) or x.is_Symbol):
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('!'))
        return pform

    def _print_factorial(self, e):
        x = e.args[0]
        pform = self._print(x)
        # Add parentheses if needed
        if not ((x.is_Integer and x.is_nonnegative) or x.is_Symbol):
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.right('!'))
        return pform

    def _print_factorial2(self, e):
        x = e.args[0]
        pform = self._print(x)
        # Add parentheses if needed
        if not ((x.is_Integer and x.is_nonnegative) or x.is_Symbol):
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.right('!!'))
        return pform

    def _print_binomial(self, e):
        n, k = e.args

        n_pform = self._print(n)
        k_pform = self._print(k)

        bar = ' '*max(n_pform.width(), k_pform.width())

        pform = prettyForm(*k_pform.above(bar))
        pform = prettyForm(*pform.above(n_pform))
        pform = prettyForm(*pform.parens('(', ')'))

        pform.baseline = (pform.baseline + 1)//2

        return pform

    def _print_Relational(self, e):
        op = prettyForm(' ' + xsym(e.rel_op) + ' ')

        l = self._print(e.lhs)
        r = self._print(e.rhs)
        pform = prettyForm(*stringPict.next(l, op, r))
        return pform

    def _print_Not(self, e):
        from ...logic import Equivalent, Implies
        if self._use_unicode:
            arg = e.args[0]
            pform = self._print(arg)
            if isinstance(arg, Equivalent):
                return self._print_Equivalent(arg, altchar='\N{NOT IDENTICAL TO}')
            if isinstance(arg, Implies):
                return self._print_Implies(arg, altchar='\N{RIGHTWARDS ARROW WITH STROKE}')

            if arg.is_Boolean and not arg.is_Not:
                pform = prettyForm(*pform.parens())

            return prettyForm(*pform.left('\N{NOT SIGN}'))
        else:
            return self._print_Function(e)

    def __print_Boolean(self, e, char, sort=True):
        args = e.args
        if sort:
            args = sorted(e.args, key=default_sort_key)
        arg = args[0]
        pform = self._print(arg)

        if arg.is_Boolean and (not arg.is_Not and not arg.is_Atom):
            pform = prettyForm(*pform.parens())

        for arg in args[1:]:
            pform_arg = self._print(arg)

            if arg.is_Boolean and (not arg.is_Not and not arg.is_Atom):
                pform_arg = prettyForm(*pform_arg.parens())

            pform = prettyForm(*pform.right(f' {char} '))
            pform = prettyForm(*pform.right(pform_arg))

        return pform

    def _print_And(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, '\N{LOGICAL AND}')
        else:
            return self._print_Function(e, sort=True)

    def _print_Or(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, '\N{LOGICAL OR}')
        else:
            return self._print_Function(e, sort=True)

    def _print_Xor(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, '\N{XOR}')
        else:
            return self._print_Function(e, sort=True)

    def _print_Nand(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, '\N{NAND}')
        else:
            return self._print_Function(e, sort=True)

    def _print_Nor(self, e):
        if self._use_unicode:
            return self.__print_Boolean(e, '\N{NOR}')
        else:
            return self._print_Function(e, sort=True)

    def _print_Implies(self, e, altchar=None):
        if self._use_unicode:
            return self.__print_Boolean(e, altchar or '\N{RIGHTWARDS ARROW}', sort=False)
        else:
            return self._print_Function(e)

    def _print_Equivalent(self, e, altchar=None):
        if self._use_unicode:
            return self.__print_Boolean(e, altchar or '\N{IDENTICAL TO}')
        else:
            return self._print_Function(e, sort=True)

    def _print_conjugate(self, e):
        pform = self._print(e.args[0])
        return prettyForm( *pform.above( hobj('_', pform.width())) )

    def _print_Abs(self, e):
        pform = self._print(e.args[0])
        pform = prettyForm(*pform.parens('|', '|'))
        return pform
    _print_Determinant = _print_Abs

    def _print_floor(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lfloor', 'rfloor'))
            return pform
        else:
            return self._print_Function(e)

    def _print_ceiling(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens('lceil', 'rceil'))
            return pform
        else:
            return self._print_Function(e)

    def _print_Derivative(self, deriv):
        if requires_partial(deriv) and self._use_unicode:
            deriv_symbol = '\N{PARTIAL DIFFERENTIAL}'
        else:
            deriv_symbol = r'd'
        syms = list(reversed(deriv.variables))
        x = None

        for sym, num in group(syms, multiple=False):
            s = self._print(sym)
            ds = prettyForm(*s.left(deriv_symbol))

            if num > 1:
                ds = ds**prettyForm(str(num))

            if x is None:
                x = ds
            else:
                x = prettyForm(*x.right(' '))
                x = prettyForm(*x.right(ds))

        f = prettyForm(
            binding=prettyForm.FUNC, *self._print(deriv.expr).parens())

        pform = prettyForm(deriv_symbol)

        if len(syms) > 1:
            pform = pform**prettyForm(str(len(syms)))

        pform = prettyForm(*pform.below(stringPict.LINE, x))
        pform.baseline = pform.baseline + 1
        pform = prettyForm(*stringPict.next(pform, f))
        pform.binding = prettyForm.MUL

        return pform

    def _print_Integral(self, integral):
        f = integral.function

        # Add parentheses if arg involves addition of terms and
        # create a pretty form for the argument
        prettyF = self._print(f)
        # XXX generalize parens
        if f.is_Add:
            prettyF = prettyForm(*prettyF.parens())

        # dx dy dz ...
        arg = prettyF
        for x in integral.limits:
            prettyArg = self._print(x[0])
            # XXX qparens (parens if needs-parens)
            if prettyArg.width() > 1:
                prettyArg = prettyForm(*prettyArg.parens())

            arg = prettyForm(*arg.right(' d', prettyArg))

        # \int \int \int ...
        firstterm = True
        s = None
        for lim in integral.limits:
            x = lim[0]
            # Create bar based on the height of the argument
            h = arg.height()
            H = h + 2

            # XXX hack!
            ascii_mode = not self._use_unicode
            if ascii_mode:
                H += 2

            vint = vobj('int', H)

            # Construct the pretty form with the integral sign and the argument
            pform = prettyForm(vint)
            pform.baseline = arg.baseline + (
                H - h)//2    # covering the whole argument

            if len(lim) > 1:
                # Create pretty forms for endpoints, if definite integral.
                # Do not print empty endpoints.
                if len(lim) == 2:
                    prettyA = prettyForm('')
                    prettyB = self._print(lim[1])
                if len(lim) == 3:
                    prettyA = self._print(lim[1])
                    prettyB = self._print(lim[2])

                if ascii_mode:  # XXX hack
                    # Add spacing so that endpoint can more easily be
                    # identified with the correct integral sign
                    spc = max(1, 3 - prettyB.width())
                    prettyB = prettyForm(*prettyB.left(' ' * spc))

                    spc = max(1, 4 - prettyA.width())
                    prettyA = prettyForm(*prettyA.right(' ' * spc))

                pform = prettyForm(*pform.above(prettyB))
                pform = prettyForm(*pform.below(prettyA))

            if not ascii_mode:  # XXX hack
                pform = prettyForm(*pform.right(' '))

            if firstterm:
                s = pform   # first term
                firstterm = False
            else:
                s = prettyForm(*s.left(pform))

        pform = prettyForm(*arg.left(s))
        pform.binding = prettyForm.MUL
        return pform

    def _print_Product(self, expr):
        func = expr.term
        pretty_func = self._print(func)

        horizontal_chr = xobj('_', 1)
        corner_chr = xobj('_', 1)
        vertical_chr = xobj('|', 1)

        if self._use_unicode:
            # use unicode corners
            horizontal_chr = xobj('-', 1)
            corner_chr = '\N{BOX DRAWINGS LIGHT DOWN AND HORIZONTAL}'

        func_height = pretty_func.height()

        first = True
        max_upper = 0
        sign_height = 0

        for lim in expr.limits:
            width = (func_height + 2) * 5 // 3 - 2
            sign_lines = []
            sign_lines.append(corner_chr + (horizontal_chr*width) + corner_chr)
            for i in range(func_height + 1):
                sign_lines.append(vertical_chr + (' '*width) + vertical_chr)

            pretty_sign = stringPict('')
            pretty_sign = prettyForm(*pretty_sign.stack(*sign_lines))

            pretty_upper = self._print(lim[2])
            pretty_lower = self._print(Equality(lim[0], lim[1]))

            max_upper = max(max_upper, pretty_upper.height())

            if first:
                sign_height = pretty_sign.height()

            pretty_sign = prettyForm(*pretty_sign.above(pretty_upper))
            pretty_sign = prettyForm(*pretty_sign.below(pretty_lower))

            if first:
                pretty_func.baseline = 0
                first = False

            height = pretty_sign.height()
            padding = stringPict('')
            padding = prettyForm(*padding.stack(*[' ']*(height - 1)))
            pretty_sign = prettyForm(*pretty_sign.right(padding))

            pretty_func = prettyForm(*pretty_sign.right(pretty_func))

        pretty_func.baseline = max_upper + sign_height//2
        pretty_func.binding = prettyForm.MUL
        return pretty_func

    def _print_Sum(self, expr):
        ascii_mode = not self._use_unicode

        def asum(hrequired, lower, upper, use_ascii):
            h = max(hrequired, 2)
            d = h//2
            w = d + 1
            more = hrequired % 2

            lines = []
            if use_ascii:
                lines.append('_'*w + ' ')
                lines.append(r'\%s`' % (' '*(w - 1)))
                for i in range(1, d):
                    lines.append('%s\\%s' % (' '*i, ' '*(w - i)))
                if more:
                    lines.append('%s)%s' % (' '*d, ' '*(w - d)))
                for i in reversed(range(1, d)):
                    lines.append('%s/%s' % (' '*i, ' '*(w - i)))
                lines.append('/' + '_'*(w - 1) + ',')
                return d, h + more, lines, 0
            else:
                w = w + more
                d = d + more
                vsum = vobj('sum', 4)
                lines.append('_'*w)
                for i in range(d):
                    lines.append('%s%s%s' % (' '*i, vsum[2], ' '*(w - i - 1)))
                for i in reversed(range(d)):
                    lines.append('%s%s%s' % (' '*i, vsum[4], ' '*(w - i - 1)))
                lines.append(vsum[8]*w)
                return d, h + 2*more, lines, more

        f = expr.function

        prettyF = self._print(f)

        if f.is_Add:  # add parens
            prettyF = prettyForm(*prettyF.parens())

        H = prettyF.height() + 2

        # \sum \sum \sum ...
        first = True
        max_upper = 0
        sign_height = 0

        for lim in expr.limits:
            assert len(lim) == 3  # no support for indefinite sums yet
            prettyUpper = self._print(lim[2])
            prettyLower = self._print(Equality(lim[0], lim[1],
                                               evaluate=False))

            max_upper = max(max_upper, prettyUpper.height())

            # Create sum sign based on the height of the argument
            d, h, slines, adjustment = asum(
                H, prettyLower.width(), prettyUpper.width(), ascii_mode)
            prettySign = stringPict('')
            prettySign = prettyForm(*prettySign.stack(*slines))

            if first:
                sign_height = prettySign.height()

            prettySign = prettyForm(*prettySign.above(prettyUpper))
            prettySign = prettyForm(*prettySign.below(prettyLower))

            if first:
                # change F baseline so it centers on the sign
                prettyF.baseline -= d - (prettyF.height()//2 -
                                         prettyF.baseline) - adjustment
                first = False

            # put padding to the right
            pad = stringPict('')
            pad = prettyForm(*pad.stack(*[' ']*h))
            prettySign = prettyForm(*prettySign.right(pad))
            # put the present prettyF to the right
            prettyF = prettyForm(*prettySign.right(prettyF))

        prettyF.baseline = max_upper + sign_height//2
        prettyF.binding = prettyForm.MUL
        return prettyF

    def _print_Limit(self, l):
        e, z, z0, dir = l.args

        E = self._print(e)
        if e.is_Add or e.is_Relational:
            E = prettyForm(*E.parens())
        Lim = prettyForm('lim')

        LimArg = self._print(z)
        if self._use_unicode:
            LimArg = prettyForm(*LimArg.right('\N{BOX DRAWINGS LIGHT HORIZONTAL}\N{RIGHTWARDS ARROW}'))
        else:
            LimArg = prettyForm(*LimArg.right('->'))
        LimArg = prettyForm(*LimArg.right(self._print(z0)))

        if str(dir) == 'real' or z0 in (oo, -oo):
            dir = ''
        else:
            if self._use_unicode:
                dir = '\N{SUPERSCRIPT PLUS SIGN}' if str(dir) == '+' else '\N{SUPERSCRIPT MINUS}'

        LimArg = prettyForm(*LimArg.right(self._print(dir)))

        Lim = prettyForm(*Lim.below(LimArg))
        Lim = prettyForm(*Lim.right(E))

        return Lim

    def _print_matrix_contents(self, e):
        """
        This method factors out what is essentially grid printing.
        """
        M = e   # matrix
        Ms = {}  # i,j -> pretty(M[i,j])
        for i in range(M.rows):
            for j in range(M.cols):
                Ms[i, j] = self._print(M[i, j])

        # h- and v- spacers
        hsep = 2
        vsep = 1

        # max width for columns
        maxw = [-1] * M.cols

        for j in range(M.cols):
            maxw[j] = max([Ms[i, j].width() for i in range(M.rows)] or [0])

        # drawing result
        D = None

        for i in range(M.rows):

            D_row = None
            for j in range(M.cols):
                s = Ms[i, j]

                # reshape s to maxw
                # XXX this should be generalized, and go to stringPict.reshape ?
                assert s.width() <= maxw[j]

                # hcenter it, +0.5 to the right                        2
                # ( it's better to align formula starts for say 0 and r )
                # XXX this is not good in all cases -- maybe introduce vbaseline?
                wdelta = maxw[j] - s.width()
                wleft = wdelta // 2
                wright = wdelta - wleft

                s = prettyForm(*s.right(' '*wright))
                s = prettyForm(*s.left(' '*wleft))

                # we don't need vcenter cells -- this is automatically done in
                # a pretty way because when their baselines are taking into
                # account in .right()

                if D_row is None:
                    D_row = s   # first box in a row
                    continue

                D_row = prettyForm(*D_row.right(' '*hsep))  # h-spacer
                D_row = prettyForm(*D_row.right(s))

            if D is None:
                D = D_row       # first row in a picture
                continue

            # v-spacer
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))

            D = prettyForm(*D.below(D_row))

        if D is None:
            D = prettyForm('')  # Empty Matrix

        return D

    def _print_MatrixBase(self, e):
        D = self._print_matrix_contents(e)
        D.baseline = D.height()//2
        D = prettyForm(*D.parens('[', ']'))
        return D

    def _print_Trace(self, e):
        D = self._print(e.arg)
        D = prettyForm(*D.parens('(', ')'))
        D.baseline = D.height()//2
        D = prettyForm(*D.left('\n'*0 + 'tr'))
        return D

    def _print_MatrixElement(self, expr):
        from ...matrices import MatrixSymbol
        if (isinstance(expr.parent, MatrixSymbol)
                and expr.i.is_number and expr.j.is_number):
            return self._print(
                Symbol(expr.parent.name + f'_{expr.i:d}{expr.j:d}'))
        else:
            prettyFunc = self._print(expr.parent)
            prettyIndices = self._print_seq((expr.i, expr.j),
                                            delimiter=', ').parens(left='[',
                                                                   right=']')[0]
            pform = prettyForm(binding=prettyForm.FUNC,
                               *stringPict.next(prettyFunc, prettyIndices))

            # store pform parts so it can be reassembled e.g. when powered
            pform.prettyFunc = prettyFunc
            pform.prettyArgs = prettyIndices

            return pform

    def _print_Transpose(self, expr):
        pform = self._print(expr.arg)
        from ...matrices import MatrixSymbol
        if not isinstance(expr.arg, MatrixSymbol):
            pform = prettyForm(*pform.parens())
        pform = pform**(prettyForm('T'))
        return pform

    def _print_Adjoint(self, expr):
        pform = self._print(expr.arg)
        if self._use_unicode:
            dag = prettyForm('\N{DAGGER}')
        else:
            dag = prettyForm('+')
        from ...matrices import MatrixSymbol
        if not isinstance(expr.arg, MatrixSymbol):
            pform = prettyForm(*pform.parens())
        pform = pform**dag
        return pform

    def _print_MatAdd(self, expr):
        return self._print_seq(expr.args, None, None, ' + ')

    def _print_MatMul(self, expr):
        args = list(expr.args)
        from ...matrices import HadamardProduct, MatAdd
        for i, a in enumerate(args):
            if (isinstance(a, (Add, MatAdd, HadamardProduct))
                    and len(expr.args) > 1):
                args[i] = prettyForm(*self._print(a).parens())
            else:
                args[i] = self._print(a)

        return prettyForm.__mul__(*args)

    def _print_MatPow(self, expr):
        pform = self._print(expr.base)
        from ...matrices import MatrixSymbol
        if not isinstance(expr.base, MatrixSymbol):
            pform = prettyForm(*pform.parens())
        pform = pform**(self._print(expr.exp))
        return pform

    _print_MatrixSymbol = _print_Symbol

    def _print_BasisDependent(self, expr):
        from ...vector import Vector

        if not self._use_unicode:
            raise NotImplementedError('ASCII pretty printing of BasisDependent is not implemented')

        if expr == expr.zero:
            return prettyForm(expr.zero._pretty_form)
        o1 = []
        vectstrs = []
        if isinstance(expr, Vector):
            items = expr.separate().items()
        else:
            items = [(0, expr)]
        for system, vect in items:
            inneritems = list(vect.components.items())
            inneritems.sort(key=lambda x: x[0].__str__())
            for k, v in inneritems:
                # If the coef of the basis vector is 1 - we skip the 1
                if v == 1:
                    o1.append('' + k._pretty_form)

                # Same for -1
                elif v == -1:
                    o1.append('(-1) ' + k._pretty_form)

                else:
                    # We always wrap the measure numbers in #parentheses
                    arg_str = self._print(v).parens()[0]

                    o1.append(arg_str + ' ' + k._pretty_form)
                vectstrs.append(k._pretty_form)

        # Fixing the newlines
        lengths = []
        strs = ['']
        for i, partstr in enumerate(o1):
            # XXX: What is this hack?
            if '\n' in partstr:
                tempstr = partstr
                tempstr = tempstr.replace(vectstrs[i], '')
                tempstr = tempstr.replace('\N{RIGHT PARENTHESIS UPPER HOOK}',
                                          '\N{RIGHT PARENTHESIS UPPER HOOK}'
                                          + ' ' + vectstrs[i])
                o1[i] = tempstr
        o1 = [x.split('\n') for x in o1]
        n_newlines = max(len(x) for x in o1)
        for parts in o1:
            lengths.append(len(parts[0]))
            for j in range(n_newlines):
                if j+1 <= len(parts):
                    if j >= len(strs):
                        strs.append(' ' * (sum(lengths[:-1]) +
                                           3*(len(lengths) - 1)))
                    if j == 0:
                        strs[0] += parts[0] + ' + '
                    else:
                        strs[j] += parts[j] + ' ' * (lengths[-1] -
                                                     len(parts[j]) + 3)
                else:
                    if j >= len(strs):
                        strs.append(' ' * (sum(lengths[:-1]) +
                                           3*(len(lengths) - 1)))
                    strs[j] += ' '*(lengths[-1]+3)

        return prettyForm('\n'.join([s[:-3] for s in strs]))

    def _print_NDimArray(self, expr):
        from ...matrices import ImmutableMatrix

        if expr.rank() == 0:
            return self._print(expr[()])

        level_str = [[]] + [[] for i in range(expr.rank())]
        shape_ranges = [list(range(i)) for i in expr.shape]
        for outer_i in itertools.product(*shape_ranges):
            level_str[-1].append(expr[outer_i])
            even = True
            for back_outer_i in range(expr.rank()-1, -1, -1):
                if len(level_str[back_outer_i+1]) < expr.shape[back_outer_i]:
                    break
                if even:
                    level_str[back_outer_i].append(level_str[back_outer_i+1])
                else:
                    level_str[back_outer_i].append(ImmutableMatrix(level_str[back_outer_i+1]))
                    if len(level_str[back_outer_i + 1]) == 1:
                        level_str[back_outer_i][-1] = ImmutableMatrix([[level_str[back_outer_i][-1]]])
                even = not even
                level_str[back_outer_i+1] = []

        out_expr = level_str[0][0]
        if expr.rank() % 2 == 1:
            out_expr = ImmutableMatrix([out_expr])

        return self._print(out_expr)

    _print_ImmutableDenseNDimArray = _print_NDimArray
    _print_ImmutableSparseNDimArray = _print_NDimArray
    _print_MutableDenseNDimArray = _print_NDimArray
    _print_MutableSparseNDimArray = _print_NDimArray

    def _print_Piecewise(self, pexpr):

        P = {}
        for n, ec in enumerate(pexpr.args):
            P[n, 0] = self._print(ec.expr)
            if ec.cond == true:
                P[n, 1] = prettyForm('otherwise')
            else:
                P[n, 1] = prettyForm(
                    *prettyForm('for ').right(self._print(ec.cond)))
        hsep = 2
        vsep = 1
        len_args = len(pexpr.args)

        # max widths
        maxw = [max(P[i, j].width() for i in range(len_args))
                for j in range(2)]

        # FIXME: Refactor this code and matrix into some tabular environment.
        # drawing result
        D = None

        for i in range(len_args):
            D_row = None
            for j in range(2):
                p = P[i, j]
                assert p.width() <= maxw[j]

                wdelta = maxw[j] - p.width()
                wleft = wdelta // 2
                wright = wdelta - wleft

                p = prettyForm(*p.right(' '*wright))
                p = prettyForm(*p.left(' '*wleft))

                if D_row is None:
                    D_row = p
                    continue

                D_row = prettyForm(*D_row.right(' '*hsep))  # h-spacer
                D_row = prettyForm(*D_row.right(p))
            if D is None:
                D = D_row       # first row in a picture
                continue

            # v-spacer
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))

            D = prettyForm(*D.below(D_row))

        D = prettyForm(*D.parens('{', ''))
        D.baseline = D.height()//2
        D.binding = prettyForm.OPEN
        return D

    def _hprint_vec(self, v):
        D = None

        for a in v:
            p = a
            if D is None:
                D = p
            else:
                D = prettyForm(*D.right(', '))
                D = prettyForm(*D.right(p))
        if D is None:
            D = stringPict(' ')

        return D

    def _hprint_vseparator(self, p1, p2):
        tmp = prettyForm(*p1.right(p2))
        sep = stringPict(vobj('|', tmp.height()), baseline=tmp.baseline)
        return prettyForm(*p1.right(sep, p2))

    def _print_hyper(self, e):
        # FIXME refactor Matrix, Piecewise, and this into a tabular environment
        ap = [self._print(a) for a in e.ap]
        bq = [self._print(b) for b in e.bq]

        P = self._print(e.argument)
        P.baseline = P.height()//2

        # Drawing result - first create the ap, bq vectors
        D = None
        for v in [ap, bq]:
            D_row = self._hprint_vec(v)
            if D is None:
                D = D_row       # first row in a picture
            else:
                D = prettyForm(*D.below(' '))
                D = prettyForm(*D.below(D_row))

        # make sure that the argument `z' is centred vertically
        D.baseline = D.height()//2

        # insert horizontal separator
        P = prettyForm(*P.left(' '))
        D = prettyForm(*D.right(' '))

        # insert separating `|`
        D = self._hprint_vseparator(D, P)

        # add parens
        D = prettyForm(*D.parens('(', ')'))

        # create the F symbol
        above = D.height()//2 - 1
        below = D.height() - above - 1

        sz, t, b, add, img = annotated('F')
        F = prettyForm('\n' * (above - t) + img + '\n' * (below - b),
                       baseline=above + sz)
        add = (sz + 1)//2

        F = prettyForm(*F.left(self._print(len(e.ap))))
        F = prettyForm(*F.right(self._print(len(e.bq))))
        F.baseline = above + add

        D = prettyForm(*F.right(' ', D))

        return D

    def _print_meijerg(self, e):
        # FIXME refactor Matrix, Piecewise, and this into a tabular environment

        v = {}
        v[(0, 0)] = [self._print(a) for a in e.an]
        v[(0, 1)] = [self._print(a) for a in e.aother]
        v[(1, 0)] = [self._print(b) for b in e.bm]
        v[(1, 1)] = [self._print(b) for b in e.bother]

        P = self._print(e.argument)
        P.baseline = P.height()//2

        vp = {}
        for idx in v:
            vp[idx] = self._hprint_vec(v[idx])

        for i in range(2):
            maxw = max(vp[(0, i)].width(), vp[(1, i)].width())
            for j in range(2):
                s = vp[(j, i)]
                left = (maxw - s.width()) // 2
                right = maxw - left - s.width()
                s = prettyForm(*s.left(' ' * left))
                s = prettyForm(*s.right(' ' * right))
                vp[(j, i)] = s

        D1 = prettyForm(*vp[(0, 0)].right('  ', vp[(0, 1)]))
        D1 = prettyForm(*D1.below(' '))
        D2 = prettyForm(*vp[(1, 0)].right('  ', vp[(1, 1)]))
        D = prettyForm(*D1.below(D2))

        # make sure that the argument `z' is centred vertically
        D.baseline = D.height()//2

        # insert horizontal separator
        P = prettyForm(*P.left(' '))
        D = prettyForm(*D.right(' '))

        # insert separating `|`
        D = self._hprint_vseparator(D, P)

        # add parens
        D = prettyForm(*D.parens('(', ')'))

        # create the G symbol
        above = D.height()//2 - 1
        below = D.height() - above - 1

        sz, t, b, add, img = annotated('G')
        F = prettyForm('\n' * (above - t) + img + '\n' * (below - b),
                       baseline=above + sz)

        pp = self._print(len(e.ap))
        pq = self._print(len(e.bq))
        pm = self._print(len(e.bm))
        pn = self._print(len(e.an))

        def adjust(p1, p2):
            diff = p1.width() - p2.width()
            if diff == 0:
                return p1, p2
            elif diff > 0:
                return p1, prettyForm(*p2.left(' '*diff))
            else:
                return prettyForm(*p1.left(' '*-diff)), p2
        pp, pm = adjust(pp, pm)
        pq, pn = adjust(pq, pn)
        pu = prettyForm(*pm.right(', ', pn))
        pl = prettyForm(*pp.right(', ', pq))

        ht = F.baseline - above - 2
        pu = prettyForm(*pu.below('\n'*ht))
        p = prettyForm(*pu.below(pl))

        F.baseline = above
        F = prettyForm(*F.right(p))

        F.baseline = above + add

        D = prettyForm(*F.right(' ', D))

        return D

    def _print_Function(self, e, sort=False):
        # XXX works only for applied functions
        func = e.func
        args = e.args
        if sort:
            args = sorted(args, key=default_sort_key)

        func_name = func.__name__

        prettyFunc = self._print(Symbol(func_name))
        prettyArgs = prettyForm(*self._print_seq(args).parens())

        pform = prettyForm(
            binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

    def _print_GeometryEntity(self, expr):
        # GeometryEntity is based on Tuple but should not print like a Tuple
        return self.emptyPrinter(expr)

    def _print_Lambda(self, e):
        vars, expr = e.args
        if self._use_unicode:
            arrow = ' \N{RIGHTWARDS ARROW FROM BAR} '
        else:
            arrow = ' -> '
        if len(vars) == 1:
            var_form = self._print(vars[0])
        else:
            var_form = self._print(tuple(vars))

        return prettyForm(*stringPict.next(var_form, arrow, self._print(expr)), binding=8)

    def _print_Order(self, expr):
        pform = self._print(expr.expr)
        if ((expr.point and any(p != 0 for p in expr.point)) or
                len(expr.variables) > 1):
            pform = prettyForm(*pform.right('; '))
            if len(expr.variables) > 1:
                pform = prettyForm(*pform.right(self._print(expr.variables)))
            else:
                pform = prettyForm(*pform.right(self._print(expr.variables[0])))
            if self._use_unicode:
                pform = prettyForm(*pform.right(' \N{RIGHTWARDS ARROW} '))
            else:
                pform = prettyForm(*pform.right(' -> '))
            if len(expr.point) > 1:
                pform = prettyForm(*pform.right(self._print(expr.point)))
            else:
                pform = prettyForm(*pform.right(self._print(expr.point[0])))
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('O'))
        return pform

    def _print_gamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek_unicode['Gamma']))
            return pform
        else:
            return self._print_Function(e)

    def _print_uppergamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.right(', ', self._print(e.args[1])))
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek_unicode['Gamma']))
            return pform
        else:
            return self._print_Function(e)

    def _print_lowergamma(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.right(', ', self._print(e.args[1])))
            pform = prettyForm(*pform.parens())
            pform = prettyForm(*pform.left(greek_unicode['gamma']))
            return pform
        else:
            return self._print_Function(e)

    def _print_DiracDelta(self, e):
        if self._use_unicode:
            pform = self._print(e.args[0])
            pform = prettyForm(*pform.parens())
            delta = self._print(greek_unicode['delta'])
            if len(e.args) > 1 and e.args[1] != 0:
                k = prettyForm(*self._print(e.args[1]).parens())
                delta = delta**k
            pform = prettyForm(*pform.left(delta))
            return pform
        else:
            return self._print_Function(e)

    def _print_expint(self, e):
        from ...core import Function
        if e.args[0].is_Integer and self._use_unicode:
            return self._print_Function(Function(f'E_{e.args[0]}')(e.args[1]))
        return self._print_Function(e)

    def _print_Chi(self, e):
        # This needs a special case since otherwise it comes out as greek
        # letter chi...
        prettyFunc = prettyForm('Chi')
        prettyArgs = prettyForm(*self._print_seq(e.args).parens())

        pform = prettyForm(
            binding=prettyForm.FUNC, *stringPict.next(prettyFunc, prettyArgs))

        # store pform parts so it can be reassembled e.g. when powered
        pform.prettyFunc = prettyFunc
        pform.prettyArgs = prettyArgs

        return pform

    def _print_elliptic_e(self, e):
        pforma0 = self._print(e.args[0])
        if len(e.args) == 1:
            pform = pforma0
        else:
            pforma1 = self._print(e.args[1])
            pform = self._hprint_vseparator(pforma0, pforma1)
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('E'))
        return pform

    def _print_elliptic_k(self, e):
        pform = self._print(e.args[0])
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('K'))
        return pform

    def _print_elliptic_f(self, e):
        pforma0 = self._print(e.args[0])
        pforma1 = self._print(e.args[1])
        pform = self._hprint_vseparator(pforma0, pforma1)
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left('F'))
        return pform

    def _print_elliptic_pi(self, e):
        name = greek_unicode['Pi'] if self._use_unicode else 'Pi'
        pforma0 = self._print(e.args[0])
        pforma1 = self._print(e.args[1])
        if len(e.args) == 2:
            pform = self._hprint_vseparator(pforma0, pforma1)
        else:
            pforma2 = self._print(e.args[2])
            pforma = self._hprint_vseparator(pforma1, pforma2)
            pforma = prettyForm(*pforma.left('; '))
            pform = prettyForm(*pforma.left(pforma0))
        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left(name))
        return pform

    def _print_Mod(self, expr):
        pform = self._print(expr.args[0])
        if pform.binding > prettyForm.MUL:
            pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.right(' mod '))
        pform = prettyForm(*pform.right(self._print(expr.args[1])))
        pform.binding = prettyForm.OPEN
        return pform

    def _print_GoldenRatio(self, expr):
        if self._use_unicode:
            return prettyForm(pretty_symbol('phi'))
        return self._print(Symbol('GoldenRatio'))

    def _print_EulerGamma(self, expr):
        if self._use_unicode:
            return prettyForm(pretty_symbol('gamma'))
        return self._print(Symbol('EulerGamma'))

    def _print_Add(self, expr, order=None):
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = expr.as_ordered_terms(order=order or self.order)
        pforms, indices = [], []

        def pretty_negative(pform, index):
            """Prepend a minus sign to a pretty form."""
            # TODO: Move this code to prettyForm
            if index == 0:
                if pform.height() > 1:
                    pform_neg = '- '
                else:
                    pform_neg = '-'
            else:
                pform_neg = ' - '

            if pform.binding > prettyForm.NEG:
                p = stringPict(*pform.parens())
            else:
                p = pform
            p = stringPict.next(pform_neg, p)
            # Lower the binding to NEG, even if it was higher. Otherwise, it
            # will print as a + ( - (b)), instead of a - (b).
            return prettyForm(binding=prettyForm.NEG, *p)

        for i, term in enumerate(terms):
            if term.is_Mul and _coeff_isneg(term):
                coeff, other = term.as_coeff_mul(rational=False)
                pform = self._print(Mul(-coeff, *other, evaluate=False))
                pforms.append(pretty_negative(pform, i))
            elif term.is_Rational and term.denominator > 1:
                pforms.append(None)
                indices.append(i)
            elif term.is_Number and term < 0:
                pform = self._print(-term)
                pforms.append(pretty_negative(pform, i))
            elif term.is_Relational:
                pforms.append(prettyForm(*self._print(term).parens()))
            else:
                pforms.append(self._print(term))

        if indices:
            large = True

            for pform in pforms:
                if pform is not None and pform.height() > 1:
                    break
            else:
                large = False

            for i in indices:
                term, negative = terms[i], False

                if term < 0:
                    term, negative = -term, True

                if large:
                    pform = prettyForm(str(term.numerator))/prettyForm(str(term.denominator))
                else:
                    pform = self._print(term)

                if negative:
                    pform = pretty_negative(pform, i)

                pforms[i] = pform

        return prettyForm.__add__(*pforms)

    def _print_Mul(self, product):
        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order != 'none':
            args = product.as_ordered_factors()
        else:
            args = product.args

        multiple_ones = len([x for x in args if x == 1]) > 1

        # Gather terms for numerator/denominator
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

        from ...concrete import Product, Sum
        from ...functions import Piecewise
        from ...integrals import Integral

        # Convert to pretty forms. Add parens to Add instances if there
        # is more than one term in the numer/denom
        for i in range(len(a)):
            if (a[i].is_Add and len(a) > 1) or (i != len(a) - 1 and
                                                isinstance(a[i], (Integral, Piecewise, Product, Sum))):
                a[i] = prettyForm(*self._print(a[i]).parens())
            elif a[i].is_Relational:
                a[i] = prettyForm(*self._print(a[i]).parens())
            else:
                a[i] = self._print(a[i])

        for i in range(len(b)):
            if (b[i].is_Add and len(b) > 1) or (i != len(b) - 1 and
                                                isinstance(b[i], (Integral, Piecewise, Product, Sum))):
                b[i] = prettyForm(*self._print(b[i]).parens())
            else:
                b[i] = self._print(b[i])

        # Construct a pretty form
        if len(b) == 0:
            return prettyForm.__mul__(*a)
        else:
            if len(a) == 0:
                a.append( self._print(Integer(1)) )
            return prettyForm.__mul__(*a)/prettyForm.__mul__(*b)

    # A helper function for _print_Pow to print x**(1/n)
    def _print_nth_root(self, base, expt):
        bpretty = self._print(base)

        # Construct root sign, start with the \/ shape
        _zZ = xobj('/', 1)
        rootsign = xobj('\\', 1) + _zZ
        # Make exponent number to put above it
        if expt.is_Rational:
            exp = str(expt.denominator)
            if exp == '2':
                exp = ''
        else:
            exp = str(self._print(expt.args[0]))
        exp = exp.ljust(2)
        if len(exp) > 2:
            rootsign = ' '*(len(exp) - 2) + rootsign
        # Stack the exponent
        rootsign = stringPict(exp + '\n' + rootsign)
        rootsign.baseline = 0
        # Diagonal: length is one less than height of base
        linelength = bpretty.height() - 1
        diagonal = stringPict('\n'.join(
            ' '*(linelength - i - 1) + _zZ + ' '*i
            for i in range(linelength)
        ))
        # Put baseline just below lowest line: next to exp
        diagonal.baseline = linelength - 1
        # Make the root symbol
        rootsign = prettyForm(*rootsign.right(diagonal))
        # Det the baseline to match contents to fix the height
        # but if the height of bpretty is one, the rootsign must be one higher
        rootsign.baseline = max(1, bpretty.baseline)
        # build result
        s = prettyForm(hobj('_', 2 + bpretty.width()))
        s = prettyForm(*bpretty.above(s))
        s = prettyForm(*s.left(rootsign))
        return s

    def _print_Pow(self, power):
        from ...series import Limit
        from ...simplify import fraction
        b, e = power.as_base_exp()
        if power.is_commutative and not e.is_Float:
            if e == -1:
                return prettyForm('1')/self._print(b)
            n, d = fraction(e)
            if n == 1 and d.is_Atom and not e.is_Integer:
                return self._print_nth_root(b, e)
            if e.is_Rational and e < 0:
                return prettyForm('1')/self._print(Pow(b, -e, evaluate=False))

        if b.is_Relational or isinstance(b, Limit):
            return prettyForm(*self._print(b).parens()).__pow__(self._print(e))

        return self._print(b)**self._print(e)

    def __print_numer_denom(self, p, q):
        if q == 1:
            if p < 0:
                return prettyForm(str(p), binding=prettyForm.NEG)
            else:
                return prettyForm(str(p))
        elif abs(p) >= 10 and abs(q) >= 10:
            # If more than one digit in numer and denom, print larger fraction
            if p < 0:
                return prettyForm(str(p), binding=prettyForm.NEG)/prettyForm(str(q))
                # Old printing method:
                # pform = prettyForm(str(-p))/prettyForm(str(q))
                # return prettyForm(binding=prettyForm.NEG, *pform.left('- '))
            else:
                return prettyForm(str(p))/prettyForm(str(q))

    def _print_Rational(self, expr):
        result = self.__print_numer_denom(expr.numerator, expr.denominator)

        if result is not None:
            return result
        else:
            return self.emptyPrinter(expr)

    def _print_ProductSet(self, p):
        prod_char = '\N{MULTIPLICATION SIGN}' if self._use_unicode else 'x'
        return self._print_seq(p.sets, None, None, f' {prod_char} ',
                               parenthesize=lambda set: set.is_Union or
                               set.is_Intersection or set.is_ProductSet)

    def _print_FiniteSet(self, s):
        items = sorted(s.args, key=default_sort_key)
        return self._print_seq(items, '{', '}', ', ' )

    def _print_set(self, l):
        return self._print_seq(sorted(l, key=default_sort_key), '{', '}')

    def _print_Range(self, s):

        if self._use_unicode:
            dots = '\N{HORIZONTAL ELLIPSIS}'
        else:
            dots = '...'

        if s.start == -oo:
            it = iter(s)
            printset = s.start, dots, s._last_element - s.step, s._last_element
        elif s.stop is oo or len(s) > 4:
            it = iter(s)
            printset = next(it), next(it), dots, s._last_element
        else:
            printset = tuple(s)

        return self._print_seq(printset, '{', '}', ', ' )

    def _print_Interval(self, i):
        if i.left_open:
            left = '('
        else:
            left = '['

        if i.right_open:
            right = ')'
        else:
            right = ']'

        return self._print_seq(i.args[:2], left, right)

    def _print_Intersection(self, u):

        delimiter = ' %s ' % pretty_atom('Intersection', 'n')

        return self._print_seq(u.args, None, None, delimiter,
                               parenthesize=lambda set: set.is_ProductSet or
                               set.is_Union or set.is_Complement)

    def _print_Union(self, u):

        union_delimiter = ' %s ' % pretty_atom('Union', 'U')

        return self._print_seq(u.args, None, None, union_delimiter,
                               parenthesize=lambda set: set.is_ProductSet or
                               set.is_Intersection or set.is_Complement)

    def _print_SymmetricDifference(self, u):
        if not self._use_unicode:
            raise NotImplementedError('ASCII pretty printing of SymmetricDifference is not implemented')

        sym_delimeter = ' %s ' % pretty_atom('SymmetricDifference')

        return self._print_seq(u.args, None, None, sym_delimeter)

    def _print_Complement(self, u):

        delimiter = r' \ '

        return self._print_seq(u.args, None, None, delimiter,
                               parenthesize=lambda set: set.is_ProductSet or set.is_Intersection
                               or set.is_Union)

    def _print_Contains(self, e):
        var, set = e.args
        if self._use_unicode:
            el = ' \N{ELEMENT OF} '
            return prettyForm(*stringPict.next(self._print(var),
                                               el, self._print(set)), binding=8)
        else:
            return prettyForm(sstr(e))

    def _print_seq(self, seq, left=None, right=None, delimiter=', ',
                   parenthesize=lambda x: False):
        s = None

        for item in seq:
            pform = self._print(item)

            if parenthesize(item):
                pform = prettyForm(*pform.parens())
            if s is None:
                # first element
                s = pform
            else:
                s = prettyForm(*stringPict.next(s, delimiter))
                s = prettyForm(*stringPict.next(s, pform))

        if s is None:
            s = stringPict('')

        s = prettyForm(*s.parens(left, right, ifascii_nougly=True))
        return s

    def join(self, delimiter, args):
        pform = None

        for arg in args:
            if pform is None:
                pform = arg
            else:
                pform = prettyForm(*pform.right(delimiter))
                pform = prettyForm(*pform.right(arg))

        if pform is None:
            return prettyForm('')
        else:
            return pform

    def _print_list(self, l):
        return self._print_seq(l, '[', ']')

    def _print_tuple(self, t):
        if len(t) == 1:
            ptuple = prettyForm(*stringPict.next(self._print(t[0]), ','))
            return prettyForm(*ptuple.parens('(', ')', ifascii_nougly=True))
        else:
            return self._print_seq(t, '(', ')')

    def _print_Tuple(self, expr):
        return self._print_tuple(expr)

    def _print_dict(self, d):
        keys = sorted(d, key=default_sort_key)
        items = []

        for k in keys:
            K = self._print(k)
            V = self._print(d[k])
            s = prettyForm(*stringPict.next(K, ': ', V))

            items.append(s)

        return self._print_seq(items, '{', '}')

    def _print_Dict(self, d):
        return self._print_dict(d)

    def _print_frozenset(self, s):
        items = sorted(s, key=default_sort_key)
        pretty = self._print_seq(items, '{', '}')
        pretty = prettyForm(*pretty.parens('(', ')', ifascii_nougly=True))
        pretty = prettyForm(*stringPict.next(type(s).__name__, pretty))
        return pretty

    def _print_PolyElement(self, poly):
        return prettyForm(sstr(poly))

    def _print_FracElement(self, frac):
        return prettyForm(sstr(frac))

    def _print_AlgebraicElement(self, expr):
        return self._print(expr.parent.to_expr(expr))

    def _print_RootOf(self, expr):
        args = [self._print_Add(expr.expr, order='lex')]
        if expr.free_symbols:
            args += expr.args[1:]
        else:
            args += [expr.index]
        pform = prettyForm(*self._print_seq(args).parens())
        pform = prettyForm(*pform.left('RootOf'))
        return pform

    def _print_RootSum(self, expr):
        args = [self._print_Add(expr.expr, order='lex')]

        if expr.fun is not S.IdentityFunction:
            args.append(self._print(expr.fun))

        pform = prettyForm(*self._print_seq(args).parens())
        pform = prettyForm(*pform.left('RootSum'))

        return pform

    def _print_FiniteField(self, expr):
        if self._use_unicode:
            form = '\N{MATHEMATICAL DOUBLE-STRUCK CAPITAL F}_%d'
        else:
            form = 'GF(%d)'

        return prettyForm(pretty_symbol(form % expr.order))

    def _print_IntegerRing(self, expr):
        if self._use_unicode:
            return prettyForm('\N{DOUBLE-STRUCK CAPITAL Z}')
        else:
            return prettyForm('ZZ')

    def _print_RationalField(self, expr):
        if self._use_unicode:
            return prettyForm('\N{DOUBLE-STRUCK CAPITAL Q}')
        else:
            return prettyForm('QQ')

    def _print_RealField(self, domain):
        if self._use_unicode:
            prefix = '\N{DOUBLE-STRUCK CAPITAL R}'
        else:
            prefix = 'RR'

        if domain.has_default_precision:
            return prettyForm(prefix)
        else:
            return self._print(pretty_symbol(prefix + '_' + str(domain.precision)))

    def _print_PolynomialRing(self, expr):
        args = list(expr.symbols)

        if not expr.order.is_default:
            order = prettyForm(*prettyForm('order=').right(self._print(expr.order)))
            args.append(order)

        pform = self._print_seq(args, '[', ']')
        pform = prettyForm(*pform.left(self._print(expr.domain)))

        return pform

    def _print_FractionField(self, expr):
        args = list(expr.symbols)

        if not expr.order.is_default:
            order = prettyForm(*prettyForm('order=').right(self._print(expr.order)))
            args.append(order)

        pform = self._print_seq(args, '(', ')')
        pform = prettyForm(*pform.left(self._print(expr.domain)))

        return pform

    def _print_GroebnerBasis(self, basis):
        exprs = [ self._print_Add(arg, order=basis.order)
                  for arg in basis.exprs ]
        exprs = prettyForm(*self.join(', ', exprs).parens(left='[', right=']'))

        gens = [ self._print(gen) for gen in basis.gens ]

        domain = prettyForm(
            *prettyForm('domain=').right(self._print(basis.domain)))
        order = prettyForm(
            *prettyForm('order=').right(self._print(basis.order)))

        pform = self.join(', ', [exprs] + gens + [domain, order])

        pform = prettyForm(*pform.parens())
        pform = prettyForm(*pform.left(basis.__class__.__name__))

        return pform

    def _print_Subs(self, e):
        pform = self._print(e.expr)
        pform = prettyForm(*pform.parens())

        h = pform.height() if pform.height() > 1 else 2
        rvert = stringPict(vobj('|', h), baseline=pform.baseline)
        pform = prettyForm(*pform.right(rvert))

        b = pform.baseline
        pform.baseline = pform.height() - 1
        pform = prettyForm(*pform.right(self._print_seq([
            self._print_seq((self._print(v[0]), xsym('=='), self._print(v[1])),
                            delimiter='') for v in zip(e.variables, e.point) ])))

        pform.baseline = b
        return pform

    def _print_euler(self, e):
        pform = prettyForm('E')
        arg = self._print(e.args[0])
        pform_arg = prettyForm(' '*arg.width())
        pform_arg = prettyForm(*pform_arg.below(arg))
        pform = prettyForm(*pform.right(pform_arg))
        return pform

    def _print_catalan(self, e):
        pform = prettyForm('C')
        arg = self._print(e.args[0])
        pform_arg = prettyForm(' '*arg.width())
        pform_arg = prettyForm(*pform_arg.below(arg))
        pform = prettyForm(*pform.right(pform_arg))
        return pform

    def _print_KroneckerDelta(self, e):
        pform = self._print(e.args[0])
        pform = prettyForm(*pform.right((prettyForm(','))))
        pform = prettyForm(*pform.right((self._print(e.args[1]))))
        if self._use_unicode:
            a = stringPict(pretty_symbol('delta'))
        else:
            a = stringPict('d')
        b = pform
        top = stringPict(*b.left(' '*a.width()))
        bot = stringPict(*a.right(' '*b.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))

    def _print_RandomDomain(self, d):
        pform = self._print('Domain: ')
        pform = prettyForm(*pform.right(self._print(d.as_boolean())))
        return pform

    def _print_BaseScalarField(self, field):
        string = field._coord_sys._names[field._index]
        return self._print(pretty_symbol(string))

    def _print_BaseVectorField(self, field):
        s = '\N{PARTIAL DIFFERENTIAL}' + '_' + field._coord_sys._names[field._index]
        return self._print(pretty_symbol(s))

    def _print_Tr(self, p):
        # TODO: Handle indices
        pform = self._print(p.args[0])
        pform = prettyForm(*pform.left(f'{p.__class__.__name__}('))
        pform = prettyForm(*pform.right(')'))
        return pform


def pretty(expr, **settings):
    """Returns a string containing the prettified form of expr.

    For information on keyword arguments see pretty_print function.

    """
    pp = PrettyPrinter(settings)

    # XXX: this is an ugly hack, but at least it works
    use_unicode = pp._settings['use_unicode']
    uflag = pretty_use_unicode(use_unicode)

    try:
        return pp.doprint(expr)
    finally:
        pretty_use_unicode(uflag)


def pretty_print(expr, **settings):
    """Prints expr in pretty form.

    pprint is just a shortcut for this function.


    Parameters
    ==========

    expr : expression
        the expression to print
    wrap_line : bool, optional
        line wrapping enabled/disabled, defaults to True
    num_columns : int or None, optional
        number of columns before line breaking (default to None which reads
        the terminal width), useful when using Diofant without terminal.
    use_unicode : bool or None, optional
        use unicode characters, such as the Greek letter pi instead of
        the string pi.
    full_prec : bool or string, optional
        use full precision. Default to "auto"
    order : bool or string, optional
        set to 'none' for long expressions if slow; default is None

    """
    print(pretty(expr, **settings))


#:
pprint = pretty_print
