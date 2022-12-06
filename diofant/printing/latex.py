"""
A Printer which converts an expression into its LaTeX equivalent.
"""

from __future__ import annotations

import itertools
import re
import typing

import mpmath.libmp as mlib
from mpmath.libmp import prec_to_dps

from ..core import Add, Integer, Mod, Symbol, oo
from ..core.alphabets import greeks
from ..core.function import _coeff_isneg
from ..core.operations import AssocOp
from ..core.relational import Relational
from ..sets import Reals
from ..utilities import default_sort_key, has_variety
from .conventions import requires_partial, split_super_sub
from .precedence import PRECEDENCE, precedence
from .printer import Printer


# Hand-picked functions which can be used directly in both LaTeX and MathJax
# Complete list at http://www.mathjax.org/docs/1.1/tex.html#supported-latex-commands
# This variable only contains those functions which diofant uses.
accepted_latex_functions = ['arcsin', 'arccos', 'arctan', 'sin', 'cos', 'tan',
                            'sinh', 'cosh', 'tanh', 'sqrt', 'ln', 'log', 'sec',
                            'csc',  'cot', 'coth', 're', 'im', 'frac', 'root',
                            'arg']

tex_greek_dictionary = {
    'Alpha': 'A',
    'Beta': 'B',
    'Gamma': r'\Gamma',
    'Delta': r'\Delta',
    'Epsilon': 'E',
    'Zeta': 'Z',
    'Eta': 'H',
    'Theta': r'\Theta',
    'Iota': 'I',
    'Kappa': 'K',
    'Lambda': r'\Lambda',
    'Mu': 'M',
    'Nu': 'N',
    'Xi': r'\Xi',
    'omicron': 'o',
    'Omicron': 'O',
    'Pi': r'\Pi',
    'Rho': 'P',
    'Sigma': r'\Sigma',
    'Tau': 'T',
    'Upsilon': r'\Upsilon',
    'Phi': r'\Phi',
    'Chi': 'X',
    'Psi': r'\Psi',
    'Omega': r'\Omega',
    'lamda': r'\lambda',
    'Lamda': r'\Lambda',
    'khi': r'\chi',
    'Khi': r'X',
    'varepsilon': r'\varepsilon',
    'varkappa': r'\varkappa',
    'varphi': r'\varphi',
    'varpi': r'\varpi',
    'varrho': r'\varrho',
    'varsigma': r'\varsigma',
    'vartheta': r'\vartheta',
}

other_symbols = {'aleph', 'beth', 'daleth', 'gimel', 'ell', 'eth', 'hbar',
                 'hslash', 'mho', 'wp', }

# Variable name modifiers
modifier_dict = {
    # Accents
    'mathring': lambda s: r'\mathring{'+s+r'}',
    'ddddot': lambda s: r'\ddddot{'+s+r'}',
    'dddot': lambda s: r'\dddot{'+s+r'}',
    'ddot': lambda s: r'\ddot{'+s+r'}',
    'dot': lambda s: r'\dot{'+s+r'}',
    'check': lambda s: r'\check{'+s+r'}',
    'breve': lambda s: r'\breve{'+s+r'}',
    'acute': lambda s: r'\acute{'+s+r'}',
    'grave': lambda s: r'\grave{'+s+r'}',
    'tilde': lambda s: r'\tilde{'+s+r'}',
    'hat': lambda s: r'\hat{'+s+r'}',
    'bar': lambda s: r'\bar{'+s+r'}',
    'vec': lambda s: r'\vec{'+s+r'}',
    'prime': lambda s: '{'+s+"}'",
    'prm': lambda s: '{'+s+"}'",
    # Faces
    'bold': lambda s: r'\boldsymbol{'+s+r'}',
    'bm': lambda s: r'\boldsymbol{'+s+r'}',
    # Brackets
    'norm': lambda s: r'\left\|{'+s+r'}\right\|',
    'avg': lambda s: r'\left\langle{'+s+r'}\right\rangle',
    'abs': lambda s: r'\left|{'+s+r'}\right|',
    'mag': lambda s: r'\left|{'+s+r'}\right|',
}

greek_letters_set = frozenset(greeks)


_between_two_numbers_p = (
    re.compile(r'[0-9][} ]*$'),  # search
    re.compile(r'[{ ]*[-+0-9]'),  # match
)


class LatexPrinter(Printer):
    """LaTex printer."""

    printmethod = '_latex'

    _default_settings: dict[str, typing.Any] = {
        'order': None,
        'mode': 'plain',
        'itex': False,
        'fold_frac_powers': False,
        'fold_func_brackets': False,
        'fold_short_frac': None,
        'long_frac_ratio': 2,
        'mul_symbol': None,
        'inv_trig_style': 'abbreviated',
        'mat_str': None,
        'mat_delim': '[',
        'symbol_names': {},
    }

    def __init__(self, settings=None):
        Printer.__init__(self, settings)

        valid_modes = ['inline', 'plain', 'equation', 'equation*']
        if self._settings['mode'] not in valid_modes:
            raise ValueError("'mode' must be one of 'inline', 'plain', "
                             "'equation' or 'equation*'")
        if self._settings['inv_trig_style'] not in ['power', 'full',
                                                    'abbreviated']:
            raise ValueError("'inv_trig_style' must be one of 'power', 'full'"
                             "or 'abbreviated'")

        if self._settings['fold_short_frac'] is None and \
                self._settings['mode'] == 'inline':
            self._settings['fold_short_frac'] = True

        mul_symbol_table = {
            None: r' ',
            'ldot': r' \,.\, ',
            'dot': r' \cdot ',
            'times': r' \times '
        }

        self._settings['mul_symbol_latex'] = \
            mul_symbol_table[self._settings['mul_symbol']]

        self._settings['mul_symbol_latex_numbers'] = \
            mul_symbol_table[self._settings['mul_symbol'] or 'dot']

        self._delim_dict = {'(': ')', '[': ']'}

    def parenthesize(self, item, level):
        if precedence(item) <= level:
            return r'\left(%s\right)' % self._print(item)
        else:
            return self._print(item)

    def doprint(self, expr):
        tex = Printer.doprint(self, expr)

        if self._settings['mode'] == 'plain':
            return tex
        elif self._settings['mode'] == 'inline':
            return r'$%s$' % tex
        elif self._settings['itex']:
            return r'$$%s$$' % tex
        else:
            env_str = self._settings['mode']
            return r'\begin{%s}%s\end{%s}' % (env_str, tex, env_str)

    def _needs_mul_brackets(self, expr, first=False, last=False):
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of a Mul, False otherwise. This is True for Add,
        but also for some container objects that would not need brackets
        when appearing last in a Mul, e.g. an Integral. ``last=True``
        specifies that this expr is the last to appear in a Mul.
        ``first=True`` specifies that this expr is the first to appear in a Mul.
        """
        from ..concrete import Product, Sum
        from ..integrals import Integral

        if expr.is_Add:
            return True
        elif expr.is_Relational:
            return True
        elif expr.is_Mul:
            if not first and _coeff_isneg(expr):
                return True

        if expr.is_Piecewise:
            return True

        if expr.has(Mod):
            return True

        if not last and any(expr.has(x) for x in (Integral, Product, Sum)):
            return True

        return False

    def _needs_add_brackets(self, expr):
        """
        Returns True if the expression needs to be wrapped in brackets when
        printed as part of an Add, False otherwise.  This is False for most
        things.
        """
        if expr.is_Relational:
            return True
        if expr.has(Mod):
            return True
        return False

    def _do_exponent(self, expr, exp):
        if exp is not None:
            return r'\left(%s\right)^{%s}' % (expr, exp)
        else:
            return expr

    def _print_bool(self, e):
        return r'\mathrm{%s}' % e

    _print_BooleanTrue = _print_bool
    _print_BooleanFalse = _print_bool

    def _print_NoneType(self, e):
        return r'\mathrm{%s}' % e

    def _print_Add(self, expr, order=None):
        if self.order == 'none':
            terms = list(expr.args)
        else:
            terms = expr.as_ordered_terms(order=order or self.order)

        tex = ''
        for i, term in enumerate(terms):
            if i == 0:
                pass
            elif _coeff_isneg(term):
                tex += ' - '
                term = -term
            else:
                tex += ' + '
            term_tex = self._print(term)
            if self._needs_add_brackets(term):
                term_tex = r'\left(%s\right)' % term_tex
            tex += term_tex

        return tex

    def _print_Cycle(self, expr):
        from ..combinatorics import Permutation
        if not expr or (isinstance(expr, Permutation) and not expr.list()):
            return r'\left( \right)'
        expr = Permutation(expr)
        expr_perm = expr.cyclic_form
        siz = expr.size
        if expr.array_form[-1] == siz - 1:
            expr_perm = expr_perm + [[siz - 1]]
        term_tex = ''
        for i in expr_perm:
            term_tex += str(i).replace(',', r'\;')
        term_tex = term_tex.replace('[', r'\left( ')
        term_tex = term_tex.replace(']', r'\right)')
        return term_tex

    _print_Permutation = _print_Cycle

    def _print_Float(self, expr):
        # Based off of that in StrPrinter
        dps = prec_to_dps(expr._prec)
        str_real = mlib.to_str(expr._mpf_, dps, strip_zeros=True)

        # Must always have a mul symbol (as 2.5 10^{20} just looks odd)
        # thus we use the number separator
        separator = self._settings['mul_symbol_latex_numbers']

        if 'e' in str_real:
            (mant, exp) = str_real.split('e')

            if exp[0] == '+':
                exp = exp[1:]

            return r'%s%s10^{%s}' % (mant, separator, exp)
        elif str_real == '+inf':
            return r'\infty'
        elif str_real == '-inf':
            return r'- \infty'
        else:
            return str_real

    def _print_Mul(self, expr):
        include_parens = False
        if _coeff_isneg(expr):
            expr = -expr
            tex = '- '
            if expr.is_Add:
                tex += '('
                include_parens = True
        else:
            tex = ''

        from ..simplify import fraction
        numer, denom = fraction(expr, exact=True)
        separator = self._settings['mul_symbol_latex']
        numbersep = self._settings['mul_symbol_latex_numbers']

        def convert(expr):
            if not expr.is_Mul:
                return str(self._print(expr))
            else:
                _tex = last_term_tex = ''

                if self.order != 'none':
                    args = expr.as_ordered_factors()
                else:
                    args = expr.args

                for i, term in enumerate(args):
                    term_tex = self._print(term)

                    if self._needs_mul_brackets(term, first=(i == 0),
                                                last=(i == len(args) - 1)):
                        term_tex = r'\left(%s\right)' % term_tex

                    if _between_two_numbers_p[0].search(last_term_tex) and \
                            _between_two_numbers_p[1].match(term_tex):
                        # between two numbers
                        _tex += numbersep
                    elif _tex:
                        _tex += separator

                    _tex += term_tex
                    last_term_tex = term_tex
                return _tex

        if denom == 1:
            # use the original expression here, since fraction() may have
            # altered it when producing numer and denom
            tex += convert(expr)
        else:
            snumer = convert(numer)
            sdenom = convert(denom)
            ldenom = len(sdenom.split())
            ratio = self._settings['long_frac_ratio']
            if self._settings['fold_short_frac'] \
                    and ldenom <= 2 and '^' not in sdenom:
                # handle short fractions
                if self._needs_mul_brackets(numer, last=False):
                    tex += r'\left(%s\right) / %s' % (snumer, sdenom)
                else:
                    tex += r'%s / %s' % (snumer, sdenom)
            elif len(snumer.split()) > ratio*ldenom:
                # handle long fractions
                if self._needs_mul_brackets(numer, last=True):
                    tex += r'\frac{1}{%s}%s\left(%s\right)' \
                        % (sdenom, separator, snumer)
                elif numer.is_Mul:
                    # split a long numerator
                    a = Integer(1)
                    b = Integer(1)
                    for x in numer.args:
                        if self._needs_mul_brackets(x, last=False) or \
                                len(convert(a*x).split()) > ratio*ldenom or \
                                (b.is_commutative is x.is_commutative is False):
                            b *= x
                        else:
                            a *= x
                    if self._needs_mul_brackets(b, last=True):
                        tex += r'\frac{%s}{%s}%s\left(%s\right)' \
                            % (convert(a), sdenom, separator, convert(b))
                    else:
                        tex += r'\frac{%s}{%s}%s%s' \
                            % (convert(a), sdenom, separator, convert(b))
                else:
                    tex += r'\frac{1}{%s}%s%s' % (sdenom, separator, snumer)
            else:
                tex += r'\frac{%s}{%s}' % (snumer, sdenom)

        if include_parens:
            tex += ')'

        return tex

    def _print_Pow(self, expr):
        from ..integrals import Integral

        # Treat root(x, n) as special case
        if expr.exp.is_Rational and abs(expr.exp.numerator) == 1 and expr.exp.denominator != 1:
            base = self._print(expr.base)
            expq = expr.exp.denominator

            if expq == 2:
                tex = r'\sqrt{%s}' % base
            elif self._settings['itex']:
                tex = r'\root{%d}{%s}' % (expq, base)
            else:
                tex = r'\sqrt[%d]{%s}' % (expq, base)

            if expr.exp.is_negative:
                return r'\frac{1}{%s}' % tex
            else:
                return tex
        elif self._settings['fold_frac_powers'] \
            and expr.exp.is_Rational \
                and expr.exp.denominator != 1:
            base, p, q = self.parenthesize(expr.base, PRECEDENCE['Pow']), expr.exp.numerator, expr.exp.denominator
            if expr.base.is_Function:
                return self._print(expr.base, f'{p}/{q}')
            return r'%s^{%s/%s}' % (base, p, q)
        elif expr.exp.is_Rational and expr.exp.is_negative and expr.base.is_commutative:
            # Things like 1/x
            return self._print_Mul(expr)
        else:
            if expr.base.is_Function:
                return self._print(expr.base, self._print(expr.exp))
            else:
                tex = r'%s^{%s}'
                if expr.base.is_Float or isinstance(expr.base, Integral):
                    tex = r'\left(%s\right)^{%s}'

                return tex % (self.parenthesize(expr.base, PRECEDENCE['Pow']),
                              self._print(expr.exp))

    def _print_Sum(self, expr):
        if len(expr.limits) == 1:
            tex = r'\sum_{%s=%s}^{%s} ' % \
                tuple(self._print(i) for i in expr.limits[0])
        else:
            def _format_ineq(l):
                return r'%s \leq %s \leq %s' % \
                    tuple(self._print(s) for s in (l[1], l[0], l[2]))

            tex = r'\sum_{\substack{%s}} ' % \
                str.join('\\\\', [_format_ineq(l) for l in expr.limits])

        if isinstance(expr.function, Add):
            tex += r'\left(%s\right)' % self._print(expr.function)
        else:
            tex += self._print(expr.function)

        return tex

    def _print_Product(self, expr):
        if len(expr.limits) == 1:
            tex = r'\prod_{%s=%s}^{%s} ' % \
                tuple(self._print(i) for i in expr.limits[0])
        else:
            def _format_ineq(l):
                return r'%s \leq %s \leq %s' % \
                    tuple(self._print(s) for s in (l[1], l[0], l[2]))

            tex = r'\prod_{\substack{%s}} ' % \
                str.join('\\\\', [_format_ineq(l) for l in expr.limits])

        if isinstance(expr.function, Add):
            tex += r'\left(%s\right)' % self._print(expr.function)
        else:
            tex += self._print(expr.function)

        return tex

    def _print_BasisDependent(self, expr):
        from ..vector import Vector

        o1 = []
        if expr == expr.zero:
            return expr.zero._latex_form
        if isinstance(expr, Vector):
            items = expr.separate().items()
        else:
            items = [(0, expr)]

        for _, vect in items:
            inneritems = list(vect.components.items())
            inneritems.sort(key=lambda x: str(x[0]))
            for k, v in inneritems:
                if v == 1:
                    o1.append(' + ' + k._latex_form)
                elif v == -1:
                    o1.append(' - ' + k._latex_form)
                else:
                    arg_str = '(' + LatexPrinter().doprint(v) + ')'
                    o1.append(' + ' + arg_str + k._latex_form)

        outstr = (''.join(o1))
        if outstr[1] != '-':
            outstr = outstr[3:]
        else:
            outstr = outstr[1:]
        return outstr

    def _print_Indexed(self, expr):
        tex = self._print(expr.base)+'_{%s}' % ','.join(
            map(self._print, expr.indices))
        return tex

    def _print_IndexedBase(self, expr):
        return self._print(expr.label)

    def _print_Derivative(self, expr):
        dim = len(expr.variables)
        if requires_partial(expr):
            diff_symbol = r'\partial'
        else:
            diff_symbol = r'd'

        if dim == 1:
            tex = r'\frac{%s}{%s %s}' % (diff_symbol, diff_symbol,
                                         self._print(expr.variables[0]))
        else:
            multiplicity, i, tex = [], 1, ''
            current = expr.variables[0]

            for symbol in expr.variables[1:]:
                if symbol == current:
                    i = i + 1
                else:
                    multiplicity.append((current, i))
                    current, i = symbol, 1

            multiplicity.append((current, i))

            for x, i in multiplicity:
                if i == 1:
                    tex += r'%s %s' % (diff_symbol, self._print(x))
                else:
                    tex += r'%s %s^{%s}' % (diff_symbol, self._print(x), i)

            tex = r'\frac{%s^{%s}}{%s} ' % (diff_symbol, dim, tex)

        if isinstance(expr.expr, AssocOp):
            return r'%s\left(%s\right)' % (tex, self._print(expr.expr))
        else:
            return r'%s %s' % (tex, self._print(expr.expr))

    def _print_Subs(self, subs):
        expr, old, new = subs.expr, subs.variables, subs.point
        latex_expr = self._print(expr)
        latex_old = (self._print(e) for e in old)
        latex_new = (self._print(e) for e in new)  # pragma: no branch
        latex_subs = r'\\ '.join(
            e[0] + '=' + e[1] for e in zip(latex_old, latex_new))
        return r'\left. %s \right|_{\substack{ %s }}' % (latex_expr, latex_subs)

    def _print_Integral(self, expr):
        tex, symbols = '', []

        # Only up to \iiiint exists
        if len(expr.limits) <= 4 and all(len(lim) == 1 for lim in expr.limits):
            # Use len(expr.limits)-1 so that syntax highlighters don't think
            # \" is an escaped quote
            tex = r'\i' + 'i'*(len(expr.limits) - 1) + 'nt'
            symbols = [r'\, d%s' % self._print(symbol[0])
                       for symbol in expr.limits]

        else:
            for lim in reversed(expr.limits):
                symbol = lim[0]
                tex += r'\int'

                if len(lim) > 1:
                    if self._settings['mode'] in ['equation', 'equation*'] \
                            and not self._settings['itex']:
                        tex += r'\limits'

                    if len(lim) == 3:
                        tex += '_{%s}^{%s}' % (self._print(lim[1]),
                                               self._print(lim[2]))
                    if len(lim) == 2:
                        tex += '^{%s}' % (self._print(lim[1]))

                symbols.insert(0, r'\, d%s' % self._print(symbol))

        tmpl = r'%s \left(%s\right)%s' if expr.function.is_Add else r'%s %s%s'
        return tmpl % (tex, str(self._print(expr.function)), ''.join(symbols))

    def _print_Limit(self, expr):
        e, z, z0, dir = expr.args

        tex = r'\lim_{%s \to ' % self._print(z)
        if dir == Reals or z0 in (oo, -oo):
            tex += r'%s}' % self._print(z0)
        elif dir in [1, -1]:
            tex += r'%s^%s}' % (self._print(z0), self._print('+' if dir == -1 else '-'))
        else:
            raise NotImplementedError

        if isinstance(e, (AssocOp, Relational)):
            return r'%s\left(%s\right)' % (tex, self._print(e))
        else:
            return r'%s %s' % (tex, self._print(e))

    def _hprint_Function(self, func):
        r"""
        Logic to decide how to render a function to latex
          - if it is a recognized latex name, use the appropriate latex command
          - if it is a single letter, just use that letter
          - if it is a longer name, then put \operatorname{} around it and be
            mindful of undercores in the name
        """
        func = self._deal_with_super_sub(func)

        if func in accepted_latex_functions:
            name = r'\%s' % func
        elif len(func) == 1 or func.startswith('\\'):
            name = func
        else:
            name = r'\operatorname{%s}' % func
        return name

    def _print_Function(self, expr, exp=None):
        r"""
        Render functions to LaTeX, handling functions that LaTeX knows about
        e.g., sin, cos, ... by using the proper LaTeX command (\sin, \cos, ...).
        For single-letter function names, render them as regular LaTeX math
        symbols. For multi-letter function names that LaTeX does not know
        about, (e.g., Li, sech) use \operatorname{} so that the function name
        is rendered in Roman font and LaTeX handles spacing properly.

        expr is the expression involving the function
        exp is an exponent
        """
        func = expr.func.__name__

        args = [str(self._print(arg)) for arg in expr.args]
        # How inverse trig functions should be displayed, formats are:
        # abbreviated: asin, full: arcsin, power: sin^-1
        inv_trig_style = self._settings['inv_trig_style']
        # If we are dealing with a power-style inverse trig function
        inv_trig_power_case = False
        # If it is applicable to fold the argument brackets
        can_fold_brackets = (self._settings['fold_func_brackets'] and
                             len(args) == 1)

        inv_trig_table = ['asin', 'acos', 'atan', 'acot']

        # If the function is an inverse trig function, handle the style
        if func in inv_trig_table:
            if inv_trig_style == 'abbreviated':
                pass
            elif inv_trig_style == 'full':
                func = 'arc' + func[1:]
            else:  # inv_trig_style == 'power':
                func = func[1:]
                inv_trig_power_case = True

                # Can never fold brackets if we're raised to a power
                if exp is not None:
                    can_fold_brackets = False

        if inv_trig_power_case:
            name = r'\%s^{-1}' % func
        elif exp is not None:
            name = r'%s^{%s}' % (self._hprint_Function(func), exp)
        else:
            name = self._hprint_Function(func)

        if can_fold_brackets:
            if func in accepted_latex_functions:
                # Wrap argument safely to avoid parse-time conflicts
                # with the function name itself
                name += r' {%s}'
            else:
                name += r'%s'
        else:
            name += r'{\left (%s \right )}'

        if inv_trig_power_case and exp is not None:
            name += r'^{%s}' % exp

        return name % ','.join(args)

    def _print_UndefinedFunction(self, expr):
        return self._hprint_Function(str(expr))

    def _print_FunctionClass(self, expr):
        if hasattr(expr, '_latex_no_arg'):
            return expr._latex_no_arg(self)

        return self._hprint_Function(str(expr))

    def _print_Lambda(self, expr):
        symbols, expr = expr.args

        if len(symbols) == 1:
            symbols = self._print(symbols[0])
        else:
            symbols = self._print(tuple(symbols))

        tex = r'\left( %s \mapsto %s \right)' % (symbols, self._print(expr))

        return tex

    def _print_Min(self, expr, exp=None):
        args = sorted(expr.args, key=default_sort_key)
        texargs = [r'%s' % self._print(symbol) for symbol in args]
        tex = r'\min\left(%s\right)' % ', '.join(texargs)

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_Max(self, expr, exp=None):
        args = sorted(expr.args, key=default_sort_key)
        texargs = [r'%s' % self._print(symbol) for symbol in args]
        tex = r'\max\left(%s\right)' % ', '.join(texargs)

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_floor(self, expr, exp=None):
        tex = r'\lfloor{%s}\rfloor' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_ceiling(self, expr, exp=None):
        tex = r'\lceil{%s}\rceil' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_Abs(self, expr, exp=None):
        tex = r'\left|{%s}\right|' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex
    _print_Determinant = _print_Abs

    def _print_re(self, expr, exp=None):
        tex = r'\Re{%s}' % self.parenthesize(expr.args[0], PRECEDENCE['Func'])
        return self._do_exponent(tex, exp)

    def _print_im(self, expr, exp=None):
        tex = r'\Im{%s}' % self.parenthesize(expr.args[0], PRECEDENCE['Func'])
        return self._do_exponent(tex, exp)

    def _print_Not(self, e):
        from ..logic import Equivalent, Implies
        if isinstance(e.args[0], Equivalent):
            return self._print_Equivalent(e.args[0], r'\not\equiv')
        if isinstance(e.args[0], Implies):
            return self._print_Implies(e.args[0], r'\not\Rightarrow')
        if e.args[0].is_Boolean:
            return r'\neg (%s)' % self._print(e.args[0])
        else:
            return r'\neg %s' % self._print(e.args[0])

    def _print_LogOp(self, args, char):
        arg = args[0]
        if arg.is_Boolean and (not arg.is_Not and not arg.is_Atom):
            tex = r'\left(%s\right)' % self._print(arg)
        else:
            tex = r'%s' % self._print(arg)

        for arg in args[1:]:
            if arg.is_Boolean and (not arg.is_Not and not arg.is_Atom):
                tex += r' %s \left(%s\right)' % (char, self._print(arg))
            else:
                tex += r' %s %s' % (char, self._print(arg))

        return tex

    def _print_And(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r'\wedge')

    def _print_Or(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r'\vee')

    def _print_Xor(self, e):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, r'\veebar')

    def _print_Implies(self, e, altchar=None):
        return self._print_LogOp(e.args, altchar or r'\Rightarrow')

    def _print_Equivalent(self, e, altchar=None):
        args = sorted(e.args, key=default_sort_key)
        return self._print_LogOp(args, altchar or r'\equiv')

    def _print_conjugate(self, expr, exp=None):
        tex = r'\overline{%s}' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_polar_lift(self, expr, exp=None):
        func = r'\operatorname{polar\_lift}'
        arg = r'{\left (%s \right )}' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}%s' % (func, exp, arg)
        else:
            return r'%s%s' % (func, arg)

    def _print_elliptic_k(self, expr, exp=None):
        tex = r'\left(%s\right)' % self._print(expr.args[0])
        if exp is not None:
            return r'K^{%s}%s' % (exp, tex)
        else:
            return r'K%s' % tex

    def _print_elliptic_f(self, expr, exp=None):
        tex = r'\left(%s\middle| %s\right)' % \
            (self._print(expr.args[0]), self._print(expr.args[1]))
        if exp is not None:
            return r'F^{%s}%s' % (exp, tex)
        else:
            return r'F%s' % tex

    def _print_elliptic_e(self, expr, exp=None):
        if len(expr.args) == 2:
            tex = r'\left(%s\middle| %s\right)' % \
                (self._print(expr.args[0]), self._print(expr.args[1]))
        else:
            tex = r'\left(%s\right)' % self._print(expr.args[0])
        if exp is not None:
            return r'E^{%s}%s' % (exp, tex)
        else:
            return r'E%s' % tex

    def _print_elliptic_pi(self, expr, exp=None):
        if len(expr.args) == 3:
            tex = r'\left(%s; %s\middle| %s\right)' % \
                (self._print(expr.args[0]), self._print(expr.args[1]),
                 self._print(expr.args[2]))
        else:
            tex = r'\left(%s\middle| %s\right)' % \
                (self._print(expr.args[0]), self._print(expr.args[1]))
        if exp is not None:
            return r'\Pi^{%s}%s' % (exp, tex)
        else:
            return r'\Pi%s' % tex

    def _print_uppergamma(self, expr, exp=None):
        tex = r'\left(%s, %s\right)' % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r'\Gamma^{%s}%s' % (exp, tex)
        else:
            return r'\Gamma%s' % tex

    def _print_lowergamma(self, expr, exp=None):
        tex = r'\left(%s, %s\right)' % (self._print(expr.args[0]),
                                        self._print(expr.args[1]))

        if exp is not None:
            return r'\gamma^{%s}%s' % (exp, tex)
        else:
            return r'\gamma%s' % tex

    def _print_expint(self, expr, exp=None):
        tex = r'\left(%s\right)' % self._print(expr.args[1])
        nu = self._print(expr.args[0])

        if exp is not None:
            return r'\operatorname{E}_{%s}^{%s}%s' % (nu, exp, tex)
        else:
            return r'\operatorname{E}_{%s}%s' % (nu, tex)

    def _print_fresnels(self, expr, exp=None):
        tex = r'\left(%s\right)' % self._print(expr.args[0])

        if exp is not None:
            return r'S^{%s}%s' % (exp, tex)
        else:
            return r'S%s' % tex

    def _print_fresnelc(self, expr, exp=None):
        tex = r'\left(%s\right)' % self._print(expr.args[0])

        if exp is not None:
            return r'C^{%s}%s' % (exp, tex)
        else:
            return r'C%s' % tex

    def _print_subfactorial(self, expr, exp=None):
        tex = r'!%s' % self.parenthesize(expr.args[0], PRECEDENCE['Func'])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_factorial(self, expr, exp=None):
        tex = r'%s!' % self.parenthesize(expr.args[0], PRECEDENCE['Func'])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_factorial2(self, expr, exp=None):
        tex = r'%s!!' % self.parenthesize(expr.args[0], PRECEDENCE['Func'])

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_binomial(self, expr, exp=None):
        tex = r'{\binom{%s}{%s}}' % (self._print(expr.args[0]),
                                     self._print(expr.args[1]))

        if exp is not None:
            return r'%s^{%s}' % (tex, exp)
        else:
            return tex

    def _print_RisingFactorial(self, expr, exp=None):
        n, k = expr.args
        base = r'%s' % self.parenthesize(n, PRECEDENCE['Func'])

        tex = r'{%s}^{\left(%s\right)}' % (base, self._print(k))

        return self._do_exponent(tex, exp)

    def _print_FallingFactorial(self, expr, exp=None):
        n, k = expr.args
        sub = r'%s' % self.parenthesize(k, PRECEDENCE['Func'])

        tex = r'{\left(%s\right)}_{%s}' % (self._print(n), sub)

        return self._do_exponent(tex, exp)

    def _hprint_BesselBase(self, expr, exp, sym):
        tex = r'%s' % sym

        need_exp = False
        if exp is not None:
            if tex.find('^') == -1:
                tex = r'%s^{%s}' % (tex, self._print(exp))
            else:
                need_exp = True

        tex = r'%s_{%s}\left(%s\right)' % (tex, self._print(expr.order),
                                           self._print(expr.argument))

        if need_exp:
            tex = self._do_exponent(tex, exp)
        return tex

    def _hprint_vec(self, vec):
        if len(vec) == 0:
            return ''
        s = ''
        for i in vec[:-1]:
            s += '%s, ' % self._print(i)
        s += self._print(vec[-1])
        return s

    def _print_besselj(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'J')

    def _print_besseli(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'I')

    def _print_besselk(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'K')

    def _print_bessely(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'Y')

    def _print_yn(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'y')

    def _print_jn(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'j')

    def _print_hankel1(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'H^{(1)}')

    def _print_hankel2(self, expr, exp=None):
        return self._hprint_BesselBase(expr, exp, 'H^{(2)}')

    def _hprint_airy(self, expr, exp=None, notation=''):
        tex = r'\left(%s\right)' % self._print(expr.args[0])

        if exp is not None:
            return r'%s^{%s}%s' % (notation, exp, tex)
        else:
            return r'%s%s' % (notation, tex)

    def _hprint_airy_prime(self, expr, exp=None, notation=''):
        tex = r'\left(%s\right)' % self._print(expr.args[0])

        if exp is not None:
            return r'{%s^\prime}^{%s}%s' % (notation, exp, tex)
        else:
            return r'%s^\prime%s' % (notation, tex)

    def _print_airyai(self, expr, exp=None):
        return self._hprint_airy(expr, exp, 'Ai')

    def _print_airybi(self, expr, exp=None):
        return self._hprint_airy(expr, exp, 'Bi')

    def _print_airyaiprime(self, expr, exp=None):
        return self._hprint_airy_prime(expr, exp, 'Ai')

    def _print_airybiprime(self, expr, exp=None):
        return self._hprint_airy_prime(expr, exp, 'Bi')

    def _print_hyper(self, expr, exp=None):
        tex = r'{{}_{%s}F_{%s}\left(\begin{matrix} %s \\ %s \end{matrix}' \
              r'\middle| {%s} \right)}' % \
            (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._hprint_vec(expr.ap), self._hprint_vec(expr.bq),
              self._print(expr.argument))

        if exp is not None:
            tex = r'{%s}^{%s}' % (tex, self._print(exp))
        return tex

    def _print_meijerg(self, expr, exp=None):
        tex = r'{G_{%s, %s}^{%s, %s}\left(\begin{matrix} %s & %s \\' \
              r'%s & %s \end{matrix} \middle| {%s} \right)}' % \
            (self._print(len(expr.ap)), self._print(len(expr.bq)),
              self._print(len(expr.bm)), self._print(len(expr.an)),
              self._hprint_vec(expr.an), self._hprint_vec(expr.aother),
              self._hprint_vec(expr.bm), self._hprint_vec(expr.bother),
              self._print(expr.argument))

        if exp is not None:
            tex = r'{%s}^{%s}' % (tex, self._print(exp))
        return tex

    def _print_dirichlet_eta(self, expr, exp=None):
        tex = r'\left(%s\right)' % self._print(expr.args[0])
        if exp is not None:
            return r'\eta^{%s}%s' % (self._print(exp), tex)
        return r'\eta%s' % tex

    def _print_zeta(self, expr, exp=None):
        if len(expr.args) == 2:
            tex = r'\left(%s, %s\right)' % tuple(map(self._print, expr.args))
        else:
            tex = r'\left(%s\right)' % self._print(expr.args[0])
        if exp is not None:
            return r'\zeta^{%s}%s' % (self._print(exp), tex)
        return r'\zeta%s' % tex

    def _print_lerchphi(self, expr, exp=None):
        tex = r'\left(%s, %s, %s\right)' % tuple(map(self._print, expr.args))
        if exp is None:
            return r'\Phi%s' % tex
        return r'\Phi^{%s}%s' % (self._print(exp), tex)

    def _print_polylog(self, expr, exp=None):
        s, z = map(self._print, expr.args)
        tex = r'\left(%s\right)' % z
        if exp is None:
            return r'\operatorname{Li}_{%s}%s' % (s, tex)
        return r'\operatorname{Li}_{%s}^{%s}%s' % (s, self._print(exp), tex)

    def _print_jacobi(self, expr, exp=None):
        n, a, b, x = map(self._print, expr.args)
        tex = r'P_{%s}^{\left(%s,%s\right)}\left(%s\right)' % (n, a, b, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_gegenbauer(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r'C_{%s}^{\left(%s\right)}\left(%s\right)' % (n, a, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_chebyshevt(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r'T_{%s}\left(%s\right)' % (n, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_chebyshevu(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r'U_{%s}\left(%s\right)' % (n, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_legendre(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r'P_{%s}\left(%s\right)' % (n, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_assoc_legendre(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r'P_{%s}^{\left(%s\right)}\left(%s\right)' % (n, a, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_hermite(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r'H_{%s}\left(%s\right)' % (n, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_laguerre(self, expr, exp=None):
        n, x = map(self._print, expr.args)
        tex = r'L_{%s}\left(%s\right)' % (n, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_assoc_laguerre(self, expr, exp=None):
        n, a, x = map(self._print, expr.args)
        tex = r'L_{%s}^{\left(%s\right)}\left(%s\right)' % (n, a, x)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_Ynm(self, expr, exp=None):
        n, m, theta, phi = map(self._print, expr.args)
        tex = r'Y_{%s}^{%s}\left(%s,%s\right)' % (n, m, theta, phi)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_Znm(self, expr, exp=None):
        n, m, theta, phi = map(self._print, expr.args)
        tex = r'Z_{%s}^{%s}\left(%s,%s\right)' % (n, m, theta, phi)
        if exp is not None:
            tex = r'\left(' + tex + r'\right)^{%s}' % (self._print(exp))
        return tex

    def _print_Rational(self, expr):
        if expr.denominator != 1:
            sign = ''
            p = expr.numerator
            if expr.numerator < 0:
                sign = '- '
                p = -p
            if self._settings['fold_short_frac']:
                return r'%s%d / %d' % (sign, p, expr.denominator)
            else:
                return r'%s\frac{%d}{%d}' % (sign, p, expr.denominator)
        else:
            return self._print(expr.numerator)

    def _print_Order(self, expr):
        s = self._print(expr.expr)
        if expr.point and any(p != 0 for p in expr.point) or \
           len(expr.variables) > 1:
            s += '; '
            if len(expr.variables) > 1:
                s += self._print(expr.variables)
            else:
                s += self._print(expr.variables[0])
            s += r'\rightarrow{}'
            if len(expr.point) > 1:
                s += self._print(expr.point)
            else:
                s += self._print(expr.point[0])
        return r'\mathcal{O}\left(%s\right)' % s

    def _print_Symbol(self, expr):
        if expr in self._settings['symbol_names']:
            return self._settings['symbol_names'][expr]

        return self._deal_with_super_sub(expr.name) if \
            '\\' not in expr.name else expr.name
    _print_Wild = _print_Symbol
    _print_Dummy = _print_Symbol
    _print_RandomSymbol = _print_Symbol
    _print_MatrixSymbol = _print_Symbol

    def _deal_with_super_sub(self, string):
        if '{' in string:
            return string

        name, supers, subs = split_super_sub(string)

        name = translate(name)
        supers = [translate(sup) for sup in supers]
        subs = [translate(sub) for sub in subs]

        # glue all items together:
        if len(supers) > 0:
            name += '^{%s}' % ' '.join(supers)
        if len(subs) > 0:
            name += '_{%s}' % ' '.join(subs)

        return name

    def _print_Relational(self, expr):
        if self._settings['itex']:
            gt = r'\gt'
            lt = r'\lt'
        else:
            gt = '>'
            lt = '<'

        charmap = {
            '==': '=',
            '>': gt,
            '<': lt,
            '>=': r'\geq',
            '<=': r'\leq',
            '!=': r'\neq',
        }

        return f'{self._print(expr.lhs)} {charmap[expr.rel_op]} {self._print(expr.rhs)}'

    def _print_Piecewise(self, expr):
        ecpairs = [r'%s & \text{for}\: %s' % (self._print(e), self._print(c))
                   for e, c in expr.args[:-1]]
        ecpairs.append(r'%s & \text{otherwise}' % self._print(expr.args[-1].expr))
        tex = r'\begin{cases} %s \end{cases}'
        return tex % r' \\'.join(ecpairs)

    def _print_MatrixBase(self, expr):
        lines = []

        for line in range(expr.rows):  # horrible, should be 'rows'
            lines.append(' & '.join([self._print(i) for i in expr[line, :]]))

        mat_str = self._settings['mat_str']
        if mat_str is None:
            if self._settings['mode'] == 'inline':
                mat_str = 'smallmatrix'
            else:
                if (expr.cols <= 10) is True:
                    mat_str = 'matrix'
                else:
                    mat_str = 'array'

        out_str = r'\begin{%MATSTR%}%s\end{%MATSTR%}'
        out_str = out_str.replace('%MATSTR%', mat_str)
        if mat_str == 'array':
            out_str = out_str.replace('%s', '{' + 'c'*expr.cols + '}%s')
        if self._settings['mat_delim']:
            left_delim = self._settings['mat_delim']
            right_delim = self._delim_dict[left_delim]
            out_str = r'\left' + left_delim + out_str + \
                      r'\right' + right_delim
        return out_str % r'\\'.join(lines)

    def _print_MatrixSlice(self, expr):
        def latexslice(x):
            x = list(x)
            if x[2] == 1:
                del x[2]
            if x[1] == x[0] + 1:
                del x[1]
            if x[0] == 0:
                x[0] = ''
            return ':'.join(map(self._print, x))
        return (self._print(expr.parent) + r'\left[' +
                latexslice(expr.rowslice) + ', ' +
                latexslice(expr.colslice) + r'\right]')

    def _print_BlockMatrix(self, expr):
        return self._print(expr.blocks)

    def _print_Transpose(self, expr):
        mat = expr.arg
        from ..matrices import MatrixSymbol
        if not isinstance(mat, MatrixSymbol):
            return r'\left(%s\right)^T' % self._print(mat)
        else:
            return '%s^T' % self._print(mat)

    def _print_Adjoint(self, expr):
        mat = expr.arg
        from ..matrices import MatrixSymbol
        if not isinstance(mat, MatrixSymbol):
            return r'\left(%s\right)^\dag' % self._print(mat)
        else:
            return r'%s^\dag' % self._print(mat)

    def _print_MatAdd(self, expr):
        terms = list(expr.args)
        tex = ''
        for i, term in enumerate(terms):
            if i == 0:
                pass
            elif _coeff_isneg(term):
                tex += ' - '
                term = -term
            else:
                tex += ' + '
            tex += self._print(term)
        return tex

    def _print_MatMul(self, expr):
        from ..core import Add
        from ..matrices import HadamardProduct, MatAdd

        def parens(x):
            if isinstance(x, (Add, MatAdd, HadamardProduct)):
                return r'\left(%s\right)' % self._print(x)
            return self._print(x)
        return ' '.join(map(parens, expr.args))

    def _print_Mod(self, expr, exp=None):
        if exp is not None:
            return r'\left(%s\bmod{%s}\right)^{%s}' % (self.parenthesize(expr.args[0],
                                                                         PRECEDENCE['Mul']),
                                                       self._print(expr.args[1]), self._print(exp))
        return r'%s\bmod{%s}' % (self.parenthesize(expr.args[0],
                                                   PRECEDENCE['Mul']),
                                 self._print(expr.args[1]))

    def _print_HadamardProduct(self, expr):
        from ..core import Add
        from ..matrices import MatAdd, MatMul

        def parens(x):
            if isinstance(x, (Add, MatAdd, MatMul)):
                return r'\left(%s\right)' % self._print(x)
            return self._print(x)
        return r' \circ '.join(map(parens, expr.args))

    def _print_MatPow(self, expr):
        base, exp = expr.base, expr.exp
        from ..matrices import MatrixSymbol
        if not isinstance(base, MatrixSymbol):
            return r'\left(%s\right)^{%s}' % (self._print(base), self._print(exp))
        else:
            return '%s^{%s}' % (self._print(base), self._print(exp))

    def _print_ZeroMatrix(self, Z):
        return r'\mathbb{0}'

    def _print_Identity(self, I):
        return r'\mathbb{I}'

    def _print_NDimArray(self, expr):
        mat_str = self._settings['mat_str']
        if mat_str is None:
            if self._settings['mode'] == 'inline':
                mat_str = 'smallmatrix'
            else:
                if (expr.rank() == 0) or (expr.shape[-1] <= 10):
                    mat_str = 'matrix'
                else:
                    mat_str = 'array'
        block_str = r'\begin{%MATSTR%}%s\end{%MATSTR%}'
        block_str = block_str.replace('%MATSTR%', mat_str)
        if self._settings['mat_delim']:
            left_delim = self._settings['mat_delim']
            right_delim = self._delim_dict[left_delim]
            block_str = r'\left' + left_delim + block_str + r'\right' + right_delim

        if expr.rank() == 0:
            return block_str % ''

        level_str = [[]] + [[] for i in range(expr.rank())]
        shape_ranges = [list(range(i)) for i in expr.shape]
        for outer_i in itertools.product(*shape_ranges):
            level_str[-1].append(self._print(expr[outer_i]))
            even = True
            for back_outer_i in range(expr.rank()-1, -1, -1):
                if len(level_str[back_outer_i+1]) < expr.shape[back_outer_i]:
                    break
                if even:
                    level_str[back_outer_i].append(r' & '.join(level_str[back_outer_i+1]))
                else:
                    level_str[back_outer_i].append(block_str % (r'\\'.join(level_str[back_outer_i+1])))
                    if len(level_str[back_outer_i+1]) == 1:
                        level_str[back_outer_i][-1] = r'\left[' + level_str[back_outer_i][-1] + r'\right]'
                even = not even
                level_str[back_outer_i+1] = []

        out_str = level_str[0][0]

        if expr.rank() % 2 == 1:
            out_str = block_str % out_str

        return out_str

    _print_ImmutableDenseNDimArray = _print_NDimArray
    _print_ImmutableSparseNDimArray = _print_NDimArray
    _print_MutableDenseNDimArray = _print_NDimArray
    _print_MutableSparseNDimArray = _print_NDimArray

    def _print_tuple(self, expr):
        return r'\left ( %s\right )' % \
            r', \quad '.join([self._print(i) for i in expr])

    def _print_Tuple(self, expr):
        return self._print_tuple(expr)

    def _print_list(self, expr):
        return r'\left [ %s\right ]' % \
            r', \quad '.join([self._print(i) for i in expr])

    def _print_dict(self, d):
        keys = sorted(d, key=default_sort_key)
        items = []

        for key in keys:
            val = d[key]
            items.append(f'{self._print(key)} : {self._print(val)}')

        return r'\left \{ %s\right \}' % r', \quad '.join(items)

    def _print_Dict(self, expr):
        return self._print_dict(expr)

    def _print_Infinity(self, expr):
        return r'\infty'

    def _print_NegativeInfinity(self, expr):
        return r'-\infty'

    def _print_NaN(self, expr):
        return r'\mathrm{NaN}'

    def _print_ComplexInfinity(self, expr):
        return r'\tilde{\infty}'

    def _print_Exp1(self, expr):
        return r'e'

    def _print_Pi(self, expr):
        return r'\pi'

    def _print_GoldenRatio(self, expr):
        return r'\phi'

    def _print_EulerGamma(self, expr):
        return r'\gamma'

    def _print_ImaginaryUnit(self, expr):
        return r'i'

    def _print_DiracDelta(self, expr, exp=None):
        if len(expr.args) == 1 or expr.args[1] == 0:
            tex = r'\delta\left(%s\right)' % self._print(expr.args[0])
        else:
            tex = r'\delta^{\left( %s \right)}\left( %s \right)' % (
                self._print(expr.args[1]), self._print(expr.args[0]))
        if exp:
            tex = r'\left(%s\right)^{%s}' % (tex, exp)
        return tex

    def _print_Heaviside(self, expr, exp=None):
        tex = r'\theta\left(%s\right)' % self._print(expr.args[0])
        if exp:
            tex = r'\left(%s\right)^{%s}' % (tex, exp)
        return tex

    def _print_KroneckerDelta(self, expr, exp=None):
        i = self._print(expr.args[0])
        j = self._print(expr.args[1])
        if expr.args[0].is_Atom and expr.args[1].is_Atom:
            tex = r'\delta_{%s %s}' % (i, j)
        else:
            tex = r'\delta_{%s, %s}' % (i, j)
        if exp:
            tex = r'\left(%s\right)^{%s}' % (tex, exp)
        return tex

    def _print_LeviCivita(self, expr, exp=None):
        indices = map(self._print, expr.args)
        if all(x.is_Atom for x in expr.args):
            tex = r'\varepsilon_{%s}' % ' '.join(indices)
        else:
            tex = r'\varepsilon_{%s}' % ', '.join(indices)
        if exp:
            tex = r'\left(%s\right)^{%s}' % (tex, exp)
        return tex

    def _print_ProductSet(self, p):
        if len(p.sets) > 1 and not has_variety(p.sets):
            return self._print(p.sets[0]) + '^%d' % len(p.sets)
        else:
            return r' \times '.join(self._print(set) for set in p.sets)

    def _print_RandomDomain(self, d):
        return 'Domain: ' + self._print(d.as_boolean())

    def _print_FiniteSet(self, s):
        items = sorted(s.args, key=default_sort_key)
        return self._print_set(items)

    def _print_set(self, s):
        items = sorted(s, key=default_sort_key)
        items = ', '.join(map(self._print, items))
        return r'\left\{%s\right\}' % items

    _print_frozenset = _print_set

    def _print_Range(self, s):
        if len(s) > 4:
            it = iter(s)
            printset = next(it), next(it), r'\ldots', s._last_element
        else:
            printset = tuple(s)

        return (r'\left\{'
                + r', '.join(self._print(el) for el in printset)
                + r'\right\}')

    def _print_Interval(self, i):
        if i.left_open:
            left = '('
        else:
            left = '['

        if i.right_open:
            right = ')'
        else:
            right = ']'

        return r'\left%s%s, %s\right%s' % (left, self._print(i.start),
                                           self._print(i.end), right)

    def _print_Union(self, u):
        return r' \cup '.join([self._print(i) for i in u.args])

    def _print_Complement(self, u):
        return r' \setminus '.join([self._print(i) for i in u.args])

    def _print_Intersection(self, u):
        return r' \cap '.join([self._print(i) for i in u.args])

    def _print_SymmetricDifference(self, u):
        return r' \triangle '.join([self._print(i) for i in u.args])

    def _print_EmptySet(self, e):
        return r'\emptyset'

    def _print_Naturals(self, n):
        return r'\mathbb{N}'

    def _print_Naturals0(self, n):
        return r'\mathbb{N}_0'

    def _print_Integers(self, i):
        return r'\mathbb{Z}'

    def _print_Reals(self, i):
        return r'\mathbb{R}'

    def _print_ExtendedReals(self, i):
        return r'\overline{\mathbb{R}}'

    def _print_Rationals(self, i):
        return r'\mathbb{Q}'

    def _print_ImageSet(self, s):
        return r'\left\{%s\; |\; %s \in %s\right\}' % (
            self._print(s.lamda.expr),
            ', '.join([self._print(var) for var in s.lamda.variables]),
            self._print(s.base_set))

    def _print_Contains(self, e):
        return r'%s \in %s' % tuple(self._print(a) for a in e.args)

    def _print_FiniteField(self, expr):
        return r'\mathbb{F}_{%s}' % expr.order

    def _print_IntegerRing(self, expr):
        return r'\mathbb{Z}'

    def _print_RationalField(self, expr):
        return r'\mathbb{Q}'

    def _print_RealField(self, expr):
        return r'\mathbb{R}'

    def _print_ComplexField(self, expr):
        return r'\mathbb{C}'

    def _print_PolynomialRing(self, expr):
        domain = self._print(expr.domain)
        symbols = ', '.join(map(self._print, expr.symbols))
        return r'%s\left[%s\right]' % (domain, symbols)

    def _print_FractionField(self, expr):
        domain = self._print(expr.domain)
        symbols = ', '.join(map(self._print, expr.symbols))
        return r'%s\left(%s\right)' % (domain, symbols)

    def _print_Poly(self, poly):
        cls = poly.__class__.__name__
        expr = self._print(poly.as_expr())
        gens = list(map(self._print, poly.gens))
        domain = 'domain=%s' % self._print(poly.domain)
        args = ', '.join([expr] + gens + [domain])
        return r'\operatorname{%s}{\left( %s \right)}' % (cls, args)

    def _print_RootOf(self, root):
        cls = root.__class__.__name__
        expr = self._print(root.expr)
        if root.free_symbols:
            return r'\operatorname{%s} {\left(%s, %s, %d\right)}' % (cls, expr, *root.args[1:])
        else:
            return r'\operatorname{%s} {\left(%s, %d\right)}' % (cls, expr, root.index)

    def _print_RootSum(self, expr):
        cls = expr.__class__.__name__
        args = [self._print(expr.expr)]
        args.append(self._print(expr.fun))
        return r'\operatorname{%s} {\left(%s\right)}' % (cls, ', '.join(args))

    def _print_PolyElement(self, poly):
        mul_symbol = self._settings['mul_symbol_latex']
        return poly._str(self, PRECEDENCE, '{%s}^{%d}', mul_symbol)

    def _print_FracElement(self, frac):
        if frac.denominator == 1:
            return self._print(frac.numerator)
        else:
            numer = self._print(frac.numerator)
            denom = self._print(frac.denominator)
            return r'\frac{%s}{%s}' % (numer, denom)

    def _print_euler(self, expr):
        return r'E_{%s}' % self._print(expr.args[0])

    def _print_catalan(self, expr):
        return r'C_{%s}' % self._print(expr.args[0])

    def _print_MellinTransform(self, expr):
        return r'\mathcal{M}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseMellinTransform(self, expr):
        return r'\mathcal{M}^{-1}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_LaplaceTransform(self, expr):
        return r'\mathcal{L}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseLaplaceTransform(self, expr):
        return r'\mathcal{L}^{-1}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_FourierTransform(self, expr):
        return r'\mathcal{F}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseFourierTransform(self, expr):
        return r'\mathcal{F}^{-1}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_SineTransform(self, expr):
        return r'\mathcal{SIN}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseSineTransform(self, expr):
        return r'\mathcal{SIN}^{-1}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_CosineTransform(self, expr):
        return r'\mathcal{COS}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_InverseCosineTransform(self, expr):
        return r'\mathcal{COS}^{-1}_{%s}\left[%s\right]\left(%s\right)' % (self._print(expr.args[1]), self._print(expr.args[0]), self._print(expr.args[2]))

    def _print_BaseScalarField(self, field):
        string = field._coord_sys._names[field._index]
        return r'\boldsymbol{\mathrm{%s}}' % self._print(Symbol(string))

    def _print_BaseVectorField(self, field):
        string = field._coord_sys._names[field._index]
        return r'\partial_{%s}' % self._print(Symbol(string))

    def _print_Differential(self, diff):
        field = diff._form_field
        if hasattr(field, '_coord_sys'):
            string = field._coord_sys._names[field._index]
            return r'\mathrm{d}%s' % self._print(Symbol(string))
        else:
            return 'd(%s)' % self._print(field)

    def _print_Tr(self, p):
        # Todo: Handle indices
        contents = self._print(p.args[0])
        return r'\mbox{Tr}\left(%s\right)' % contents

    def _print_totient(self, expr, exp=None):
        if exp is not None:
            return r'\left(\phi\left(%s\right)\right)^{%s}' % (self._print(expr.args[0]),
                                                               self._print(exp))
        return r'\phi\left(%s\right)' % self._print(expr.args[0])

    def _print_divisor_sigma(self, expr, exp=None):
        if len(expr.args) == 2:
            tex = r'_%s\left(%s\right)' % tuple(map(self._print,
                                                    (expr.args[1], expr.args[0])))
        else:
            tex = r'\left(%s\right)' % self._print(expr.args[0])
        if exp is not None:
            return r'\sigma^{%s}%s' % (self._print(exp), tex)
        return r'\sigma%s' % tex


def translate(s):
    r"""
    Check for a modifier ending the string.  If present, convert the
    modifier to latex and translate the rest recursively.

    Given a description of a Greek letter or other special character,
    return the appropriate latex.

    Let everything else pass as given.

    >>> translate('alphahatdotprime')
    "{\\dot{\\hat{\\alpha}}}'"
    """
    # Process the rest
    tex = tex_greek_dictionary.get(s)
    if tex:
        return tex
    elif s.lower() in greek_letters_set:
        return '\\' + s.lower()
    elif s in other_symbols:
        return '\\' + s
    else:
        # Process modifiers, if any, and recurse
        for key in sorted(modifier_dict, key=len, reverse=True):
            if s.lower().endswith(key) and len(s) > len(key):
                return modifier_dict[key](translate(s[:-len(key)]))
        return s


def latex(expr, **settings):
    r"""
    Convert the given expression to LaTeX representation.

    >>> from diofant.abc import mu, r, tau

    >>> print(latex((2*tau)**Rational(7, 2)))
    8 \sqrt{2} \tau^{\frac{7}{2}}

    Not using a print statement for printing, results in double backslashes for
    latex commands since that's the way Python escapes backslashes in strings.

    >>> latex((2*tau)**Rational(7, 2))
    '8 \\sqrt{2} \\tau^{\\frac{7}{2}}'

    order: Any of the supported monomial orderings (currently "lex", "grlex", or
    "grevlex") and "none". This parameter does nothing for Mul objects.
    For very large expressions, set the 'order' keyword to 'none' if
    speed is a concern.

    mode: Specifies how the generated code will be delimited. 'mode' can be one
    of 'plain', 'inline', 'equation' or 'equation*'.  If 'mode' is set to
    'plain', then the resulting code will not be delimited at all (this is the
    default). If 'mode' is set to 'inline' then inline LaTeX $ $ will be used.
    If 'mode' is set to 'equation' or 'equation*', the resulting code will be
    enclosed in the 'equation' or 'equation*' environment (remember to import
    'amsmath' for 'equation*'), unless the 'itex' option is set. In the latter
    case, the ``$$ $$`` syntax is used.

    >>> print(latex((2*mu)**Rational(7, 2), mode='plain'))
    8 \sqrt{2} \mu^{\frac{7}{2}}

    >>> print(latex((2*tau)**Rational(7, 2), mode='inline'))
    $8 \sqrt{2} \tau^{7 / 2}$

    >>> print(latex((2*mu)**Rational(7, 2), mode='equation*'))
    \begin{equation*}8 \sqrt{2} \mu^{\frac{7}{2}}\end{equation*}

    >>> print(latex((2*mu)**Rational(7, 2), mode='equation'))
    \begin{equation}8 \sqrt{2} \mu^{\frac{7}{2}}\end{equation}

    itex: Specifies if itex-specific syntax is used, including emitting ``$$ $$``.

    >>> print(latex((2*mu)**Rational(7, 2), mode='equation', itex=True))
    $$8 \sqrt{2} \mu^{\frac{7}{2}}$$

    fold_frac_powers: Emit "^{p/q}" instead of "^{\frac{p}{q}}" for fractional
    powers.

    >>> print(latex((2*tau)**Rational(7, 2), fold_frac_powers=True))
    8 \sqrt{2} \tau^{7/2}

    fold_func_brackets: Fold function brackets where applicable.

    >>> print(latex((2*tau)**sin(Rational(7, 2))))
    \left(2 \tau\right)^{\sin{\left (\frac{7}{2} \right )}}
    >>> print(latex((2*tau)**sin(Rational(7, 2)), fold_func_brackets=True))
    \left(2 \tau\right)^{\sin {\frac{7}{2}}}

    fold_short_frac: Emit "p / q" instead of "\frac{p}{q}" when the
    denominator is simple enough (at most two terms and no powers).
    The default value is `True` for inline mode, False otherwise.

    >>> print(latex(3*x**2/y))
    \frac{3 x^{2}}{y}
    >>> print(latex(3*x**2/y, fold_short_frac=True))
    3 x^{2} / y

    long_frac_ratio: The allowed ratio of the width of the numerator to the
    width of the denominator before we start breaking off long fractions.
    The default value is 2.

    >>> print(latex(Integral(r, r)/2/pi, long_frac_ratio=2))
    \frac{\int r\, dr}{2 \pi}
    >>> print(latex(Integral(r, r)/2/pi, long_frac_ratio=0))
    \frac{1}{2 \pi} \int r\, dr

    mul_symbol: The symbol to use for multiplication. Can be one of None,
    "ldot", "dot", or "times".

    >>> print(latex((2*tau)**sin(Rational(7, 2)), mul_symbol='times'))
    \left(2 \times \tau\right)^{\sin{\left (\frac{7}{2} \right )}}

    inv_trig_style: How inverse trig functions should be displayed. Can be one
    of "abbreviated", "full", or "power". Defaults to "abbreviated".

    >>> print(latex(asin(Rational(7, 2))))
    \operatorname{asin}{\left (\frac{7}{2} \right )}
    >>> print(latex(asin(Rational(7, 2)), inv_trig_style='full'))
    \arcsin{\left (\frac{7}{2} \right )}
    >>> print(latex(asin(Rational(7, 2)), inv_trig_style='power'))
    \sin^{-1}{\left (\frac{7}{2} \right )}

    mat_str: Which matrix environment string to emit. "smallmatrix", "matrix",
    "array", etc. Defaults to "smallmatrix" for inline mode, "matrix" for
    matrices of no more than 10 columns, and "array" otherwise.

    >>> print(latex(Matrix(2, 1, [x, y])))
    \left[\begin{matrix}x\\y\end{matrix}\right]

    >>> print(latex(Matrix(2, 1, [x, y]), mat_str='array'))
    \left[\begin{array}{c}x\\y\end{array}\right]

    mat_delim: The delimiter to wrap around matrices. Can be one of "[", "(",
    or the empty string. Defaults to "[".

    >>> print(latex(Matrix(2, 1, [x, y]), mat_delim='('))
    \left(\begin{matrix}x\\y\end{matrix}\right)

    symbol_names: Dictionary of symbols and the custom strings they should be
    emitted as.

    >>> print(latex(x**2, symbol_names={x: 'x_i'}))
    x_i^{2}

    ``latex`` also supports the builtin container types list, tuple, and
    dictionary.

    >>> print(latex([2/x, y], mode='inline'))
    $\left [ 2 / x, \quad y\right ]$

    """
    return LatexPrinter(settings).doprint(expr)
