"""
C code printer

The CCodePrinter converts single diofant expressions into single C expressions,
using the functions defined in math.h where possible.

A complete code generator, which uses ccode extensively, can be found in
diofant.utilities.codegen. The codegen module can be used to generate complete
source code files that are compilable without further modifications.
"""

from __future__ import annotations

import typing

from ..core import Integer
from ..logic import true
from .codeprinter import Assignment, CodePrinter
from .precedence import precedence


# dictionary mapping diofant function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
    'Abs': [(lambda x: not x.is_integer, 'fabs')],
    'gamma': 'tgamma',
    'sin': 'sin',
    'cos': 'cos',
    'tan': 'tan',
    'asin': 'asin',
    'acos': 'acos',
    'atan': 'atan',
    'atan2': 'atan2',
    'exp': 'exp',
    'log': 'log',
    'erf': 'erf',
    'sinh': 'sinh',
    'cosh': 'cosh',
    'tanh': 'tanh',
    'asinh': 'asinh',
    'acosh': 'acosh',
    'atanh': 'atanh',
    'floor': 'floor',
    'ceiling': 'ceil',
}

# These are the core reserved words in the C language. Taken from:
# http://crasseux.com/books/ctutorial/Reserved-words-in-C.html

reserved_words = ['auto',
                  'if',
                  'break',
                  'int',
                  'case',
                  'long',
                  'char',
                  'register',
                  'continue',
                  'return',
                  'default',
                  'short',
                  'do',
                  'sizeof',
                  'double',
                  'static',
                  'else',
                  'struct',
                  'entry',
                  'switch',
                  'extern',
                  'typedef',
                  'float',
                  'union',
                  'for',
                  'unsigned',
                  'goto',
                  'while',
                  'enum',
                  'void',
                  'const',
                  'signed',
                  'volatile']


class CCodePrinter(CodePrinter):
    """A printer to convert python expressions to strings of c code."""

    printmethod = '_ccode'
    language = 'C'

    _default_settings: dict[str, typing.Any] = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'contract': True,
        'dereference': set(),
        'error_on_reserved': False,
        'reserved_word_suffix': '_',
    }

    def __init__(self, settings={}):
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)
        self._dereference = set(settings.get('dereference', []))
        self.reserved_words = set(reserved_words)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return f'{codestring};'

    def _get_comment(self, text):
        return f'// {text}'

    def _declare_number_const(self, name, value):
        return f'double const {name} = {value};'

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for i in range(rows) for j in range(cols))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        loopstart = 'for (int %(var)s=%(start)s; %(var)s<%(end)s; %(var)s++){'
        for i in indices:
            # C arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'var': self._print(i.label),
                'start': self._print(i.lower),
                'end': self._print(i.upper + 1)})
            close_lines.append('}')
        return open_lines, close_lines

    def _print_Pow(self, expr):
        if 'Pow' in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)
        if expr.exp == -1:
            return f'1.0/{self.parenthesize(expr.base, PREC)}'
        elif expr.exp == 0.5:
            return f'sqrt({self._print(expr.base)})'
        else:
            return f'pow({self._print(expr.base)}, {self._print(expr.exp)})'

    def _print_Rational(self, expr):
        p, q = int(expr.numerator), int(expr.denominator)
        return f'{p:d}.0L/{q:d}.0L'

    def _print_Indexed(self, expr):
        # calculate index for 1d array
        dims = expr.shape
        elem = Integer(0)
        offset = Integer(1)
        for i in reversed(range(expr.rank)):
            elem += expr.indices[i]*offset
            offset *= dims[i]
        return f'{self._print(expr.base.label)}[{self._print(elem)}]'

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _print_Exp1(self, expr):
        return 'M_E'

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

    def _print_Piecewise(self, expr):
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append(f'if ({self._print(c)}) {{')
                elif i == len(expr.args) - 1 and c == true:
                    lines.append('else {')
                else:
                    lines.append(f'else if ({self._print(c)}) {{')
                code0 = self._print(e)
                lines.append(code0)
                lines.append('}')
            return '\n'.join(lines)
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = [f'(({self._print(c)}) ? (\n{self._print(e)}\n)\n'
                       for e, c in expr.args[:-1]]
            last_line = f': (\n{self._print(expr.args[-1].expr)}\n)'
            return ': '.join(ecpairs) + last_line + ' '.join([')'*len(ecpairs)])

    def _print_ITE(self, expr):
        from ..functions import Piecewise
        _piecewise = Piecewise((expr.args[1], expr.args[0]), (expr.args[2], True))
        return self._print(_piecewise)

    def _print_MatrixElement(self, expr):
        return f'{expr.parent}[{expr.j + expr.i * expr.parent.shape[1]}]'

    def _print_Symbol(self, expr):

        name = super()._print_Symbol(expr)

        if expr in self._dereference:
            return f'(*{name})'
        else:
            return name

    def _print_sign(self, func):
        e = self._print(func.args[0])
        return f'((({e}) > 0) - (({e}) < 0))'

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines."""
        if isinstance(code, str):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = '   '
        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [line.lstrip(' \t') for line in code]

        increase = [int(any(map(line.endswith, inc_token))) for line in code]
        decrease = [int(any(map(line.startswith, dec_token)))
                    for line in code]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if line in ('', '\n'):
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append(f'{tab * level}{line}')
            level += increase[n]
        return pretty


def ccode(expr, assign_to=None, **settings):
    """Converts an expr to a string of c code

    Parameters
    ==========

    expr : Expr
        A diofant expression to be converted.
    assign_to : optional
        When given, the argument is used as the name of the variable to which
        the expression is assigned. Can be a string, ``Symbol``,
        ``MatrixSymbol``, or ``Indexed`` type. This is helpful in case of
        line-wrapping, or for expressions that generate multi-line statements.
    precision : integer, optional
        The precision for numbers such as pi [default=15].
    user_functions : dict, optional
        A dictionary where the keys are string representations of either
        ``FunctionClass`` or ``UndefinedFunction`` instances and the values
        are their desired C string representations. Alternatively, the
        dictionary value can be a list of tuples i.e. [(argument_test,
        cfunction_string)] or [(argument_test, cfunction_formater)]. See below
        for examples.
    dereference : iterable, optional
        An iterable of symbols that should be dereferenced in the printed code
        expression. These would be values passed by address to the function.
        For example, if ``dereference=[a]``, the resulting code would print
        ``(*a)`` instead of ``a``.
    human : bool, optional
        If True, the result is a single string that may contain some constant
        declarations for the number symbols. If False, the same information is
        returned in a tuple of (symbols_to_declare, not_supported_functions,
        code_text). [default=True].
    contract: bool, optional
        If True, ``Indexed`` instances are assumed to obey tensor contraction
        rules and the corresponding nested loops over indices are generated.
        Setting contract=False will not generate loops, instead the user is
        responsible to provide values for the indices in the code.
        [default=True].

    Examples
    ========

    >>> ccode((2*x)**Rational(7, 2))
    '8*sqrt(2)*pow(x, 7.0L/2.0L)'
    >>> ccode(sin(x), assign_to='s')
    's = sin(x);'

    Simple custom printing can be defined for certain types by passing a
    dictionary of {"type" : "function"} to the ``user_functions`` kwarg.
    Alternatively, the dictionary value can be a list of tuples i.e.
    [(argument_test, cfunction_string)].

    >>> custom_functions = {'ceiling': 'CEIL',
    ...                     'Abs': [(lambda x: not x.is_integer, 'fabs'),
    ...                             (lambda x: x.is_integer, 'ABS')],
    ...                     'func': 'f'}
    >>> func = Function('func')
    >>> ccode(func(abs(x) + ceiling(x)), user_functions=custom_functions)
    'f(fabs(x) + CEIL(x))'

    or if the C-function takes a subset of the original arguments:

    >>> ccode(2**x + 3**x,
    ...       user_functions={'Pow': [(lambda b, e: b == 2,
    ...                                lambda b, e: f'exp2({e})'),
    ...                               (lambda b, e: b != 2, 'pow')]})
    'exp2(x) + pow(3, x)'

    ``Piecewise`` expressions are converted into conditionals. If an
    ``assign_to`` variable is provided an if statement is created, otherwise
    the ternary operator is used. Note that if the ``Piecewise`` lacks a
    default term, represented by ``(expr, True)`` then an error will be thrown.
    This is to prevent generating an expression that may not evaluate to
    anything.

    >>> expr = Piecewise((x + 1, x > 0), (x, True))
    >>> print(ccode(expr, y))
    if (x > 0) {
    y = x + 1;
    }
    else {
    y = x;
    }

    Support for loops is provided through ``Indexed`` types. With
    ``contract=True`` these expressions will be turned into loops, whereas
    ``contract=False`` will just print the assignment expression that should be
    looped over:

    >>> len_y = 5
    >>> y = IndexedBase('y', shape=[len_y])
    >>> t = IndexedBase('t', shape=[len_y])
    >>> Dy = IndexedBase('Dy', shape=[len_y - 1])
    >>> i = Idx('i', len_y-1)
    >>> e = Eq(Dy[i], (y[i+1]-y[i])/(t[i+1]-t[i]))
    >>> ccode(e.rhs, assign_to=e.lhs, contract=False)
    'Dy[i] = (y[i + 1] - y[i])/(t[i + 1] - t[i]);'

    Matrices are also supported, but a ``MatrixSymbol`` of the same dimensions
    must be provided to ``assign_to``. Note that any expression that can be
    generated normally can also exist inside a Matrix:

    >>> mat = Matrix([x**2, Piecewise((x + 1, x > 0), (x, True)), sin(x)])
    >>> A = MatrixSymbol('A', 3, 1)
    >>> print(ccode(mat, A))
    A[0] = pow(x, 2);
    if (x > 0) {
       A[1] = x + 1;
    }
    else {
       A[1] = x;
    }
    A[2] = sin(x);

    """
    return CCodePrinter(settings).doprint(expr, assign_to)
