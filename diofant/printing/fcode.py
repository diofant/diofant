"""
Fortran code printer

The FCodePrinter converts single diofant expressions into single Fortran
expressions, using the functions defined in the Fortran 77 standard where
possible. Some useful pointers to Fortran can be found on wikipedia:

https://en.wikipedia.org/wiki/Fortran

Most of the code below is based on the "Professional Programmer's Guide to
Fortran77" by Clive G. Page:

http://www.star.le.ac.uk/~cgp/prof77.html

Fortran is a case-insensitive language. This might cause trouble because
Diofant is case sensitive. The implementation below does not care and leaves
the responsibility for generating properly cased Fortran code to the user.
"""

import string

from ..core import Add, Function, I, N
from ..logic import true
from .codeprinter import Assignment, CodePrinter
from .precedence import precedence


known_functions = {
    'sin': 'sin',
    'cos': 'cos',
    'tan': 'tan',
    'asin': 'asin',
    'acos': 'acos',
    'atan': 'atan',
    'atan2': 'atan2',
    'sinh': 'sinh',
    'cosh': 'cosh',
    'tanh': 'tanh',
    'log': 'log',
    'exp': 'exp',
    'erf': 'erf',
    'Abs': 'Abs',
    'sign': 'sign',
    'conjugate': 'conjg'
}


class FCodePrinter(CodePrinter):
    """A printer to convert diofant expressions to strings of Fortran code."""

    printmethod = '_fcode'
    language = 'Fortran'

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'source_format': 'fixed',
        'contract': True,
        'standard': 77
    }

    _operators = {
        'and': '.and.',
        'or': '.or.',
        'xor': '.neqv.',
        'equivalent': '.eqv.',
        'not': '.not. ',
    }

    _relationals = {
        '!=': '/=',
    }

    def __init__(self, settings={}):
        """Initialize self."""
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)
        # leading columns depend on fixed or free format
        if self._settings['source_format'] == 'fixed':
            self._lead_code = '      '
            self._lead_cont = '     @ '
            self._lead_comment = 'C     '
        elif self._settings['source_format'] == 'free':
            self._lead_code = ''
            self._lead_cont = '      '
            self._lead_comment = '! '
        else:
            raise ValueError(f"Unknown source format: {self._settings['source_format']}")
        standards = {66, 77, 90, 95, 2003, 2008}
        if self._settings['standard'] not in standards:
            raise ValueError(f"Unknown Fortran standard: {self._settings['standard']}")

    def _rate_index_position(self, p):
        return -p*5

    def _get_statement(self, codestring):
        return codestring

    def _get_comment(self, text):
        return f'! {text}'

    def _declare_number_const(self, name, value):
        return f'parameter ({name} = {value})'

    def _format_code(self, lines):
        return self._wrap_fortran(self.indent_code(lines))

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for j in range(cols) for i in range(rows))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        for i in indices:
            # fortran arrays start at 1 and end at dimension
            var, start, stop = map(self._print,
                                   [i.label, i.lower + 1, i.upper + 1])
            open_lines.append(f'do {var} = {start}, {stop}')
            close_lines.append('end do')
        return open_lines, close_lines

    def _print_Piecewise(self, expr):
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append(f'if ({self._print(c)}) then')
                elif i == len(expr.args) - 1 and c == true:
                    lines.append('else')
                else:
                    lines.append(f'else if ({self._print(c)}) then')
                lines.append(self._print(e))
            lines.append('end if')
            return '\n'.join(lines)
        if self._settings['standard'] >= 95:
            # Only supported in F95 and newer:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            pattern = 'merge({T}, {F}, {COND})'
            code = self._print(expr.args[-1].expr)
            terms = list(expr.args[:-1])
            while terms:
                e, c = terms.pop()
                expr = self._print(e)
                cond = self._print(c)
                code = pattern.format(T=expr, F=code, COND=cond)
            return code
        # `merge` is not supported prior to F95
        raise NotImplementedError('Using Piecewise as an expression using '
                                  'inline operators is not supported in '
                                  'standards earlier than Fortran95.')

    def _print_MatrixElement(self, expr):
        return f'{expr.parent}({expr.i + 1}, {expr.j + 1})'

    def _print_Add(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        # collect the purely real and purely imaginary parts:
        pure_real = []
        pure_imaginary = []
        mixed = []
        for arg in expr.args:
            if arg.is_number and arg.is_extended_real:
                pure_real.append(arg)
            elif arg.is_number and arg.is_imaginary:
                pure_imaginary.append(arg)
            else:
                mixed.append(arg)
        if len(pure_imaginary) > 0:
            if len(mixed) > 0:
                term = Add(*mixed)
                t = self._print(term)
                if t.startswith('-'):
                    sign = '-'
                    t = t[1:]
                else:
                    sign = '+'

                return f'cmplx({self._print(Add(*pure_real))},{self._print(-I * Add(*pure_imaginary))}) {sign} {t}'
            return f'cmplx({self._print(Add(*pure_real))},{self._print(-I * Add(*pure_imaginary))})'
        return CodePrinter._print_Add(self, expr)

    def _print_Function(self, expr):
        # All constant function args are evaluated as floats
        prec = self._settings['precision']
        args = [N(a, prec) for a in expr.args]
        eval_expr = expr.func(*args)
        if not isinstance(eval_expr, Function):
            return self._print(eval_expr)
        return CodePrinter._print_Function(self, expr.func(*args))

    def _print_ImaginaryUnit(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        return 'cmplx(0,1)'

    def _print_int(self, expr):
        return str(expr)

    def _print_Mul(self, expr):
        # purpose: print complex numbers nicely in Fortran.
        if expr.is_number and expr.is_imaginary:
            return f'cmplx(0,{self._print(-I * expr)})'
        return CodePrinter._print_Mul(self, expr)

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp == -1:
            return f'1.0/{self.parenthesize(expr.base, PREC)}'
        if expr.exp == 0.5:
            if expr.base.is_integer:
                # Fortan intrinsic sqrt() does not accept integer argument
                if expr.base.is_Number:
                    return f'sqrt({self._print(expr.base)}.0d0)'
                return f'sqrt(dble({self._print(expr.base)}))'
            return f'sqrt({self._print(expr.base)})'
        return CodePrinter._print_Pow(self, expr)

    def _print_Rational(self, expr):
        p, q = int(expr.numerator), int(expr.denominator)
        return f'{p:d}.0d0/{q:d}.0d0'

    def _print_Float(self, expr):
        printed = CodePrinter._print_Float(self, expr)
        e = printed.find('e')
        if e > -1:
            return f'{printed[:e]}d{printed[e + 1:]}'
        return f'{printed}d0'

    def _print_Indexed(self, expr):
        inds = [self._print(i) for i in expr.indices]
        return f"{self._print(expr.base.label)}({', '.join(inds)})"

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _pad_leading_columns(self, lines):
        result = []
        for line in lines:
            if line.startswith('!'):
                result.append(self._lead_comment + line[1:].lstrip())
            else:
                result.append(self._lead_code + line)
        return result

    def _wrap_fortran(self, lines):
        r"""Wrap long Fortran lines

        Argument: lines  --  a list of lines (without \n character)

        A comment line is split at white space. Code lines are split with a more
        complex rule to give nice results.

        """
        # routine to find split point in a code line
        my_alnum = set('_+-.' + string.digits + string.ascii_letters)
        my_white = set(' \t()')

        def split_pos_code(line, endpos):
            if len(line) <= endpos:
                return len(line)
            pos = endpos

            def split(pos):
                return (line[pos] in my_alnum and line[pos - 1] not in my_alnum) or \
                       (line[pos] not in my_alnum and line[pos - 1] in my_alnum) or \
                       (line[pos] in my_white and line[pos - 1] not in my_white) or \
                       (line[pos] not in my_white and line[pos - 1] in my_white)

            while not split(pos):
                pos -= 1
                if pos == 0:
                    return endpos
            return pos
        # split line by line and add the splitted lines to result
        result = []
        if self._settings['source_format'] == 'free':
            trailing = ' &'
        else:
            trailing = ''
        for line in lines:
            if line.startswith(self._lead_comment):
                # comment line
                if len(line) > 72:
                    pos = line.rfind(' ', 6, 72)
                    if pos == -1:
                        pos = 72
                    hunk = line[:pos]
                    line = line[pos:].lstrip()
                    result.append(hunk)
                    while len(line) > 0:
                        pos = line.rfind(' ', 0, 66)
                        if pos == -1 or len(line) < 66:
                            pos = 66
                        hunk = line[:pos]
                        line = line[pos:].lstrip()
                        result.append(f'{self._lead_comment}{hunk}')
                else:
                    result.append(line)
            elif line.startswith(self._lead_code):
                # code line
                pos = split_pos_code(line, 72)
                hunk = line[:pos].rstrip()
                line = line[pos:].lstrip()
                if line:
                    hunk += trailing
                result.append(hunk)
                while len(line) > 0:
                    pos = split_pos_code(line, 65)
                    hunk = line[:pos].rstrip()
                    line = line[pos:].lstrip()
                    if line:
                        hunk += trailing
                    result.append(f'{self._lead_cont}{hunk}')
            else:
                result.append(line)
        return result

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines."""
        if isinstance(code, str):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        free = self._settings['source_format'] == 'free'
        code = [line.lstrip(' \t') for line in code]

        inc_keyword = ('do ', 'if(', 'if ', 'do\n', 'else')
        dec_keyword = ('end do', 'enddo', 'end if', 'endif', 'else')

        increase = [int(any(map(line.startswith, inc_keyword)))
                    for line in code]
        decrease = [int(any(map(line.startswith, dec_keyword)))
                    for line in code]
        continuation = [int(any(map(line.endswith, ['&', '&\n'])))
                        for line in code]

        level = 0
        cont_padding = 0
        tabwidth = 3
        new_code = []
        for i, line in enumerate(code):
            if line in ('', '\n'):
                new_code.append(line)
                continue
            level -= decrease[i]

            if free:
                padding = ' '*(level*tabwidth + cont_padding)
            else:
                padding = ' '*level*tabwidth

            line = f'{padding}{line}'
            if not free:
                line = self._pad_leading_columns([line])[0]

            new_code.append(line)

            if continuation[i]:
                cont_padding = 2*tabwidth
            else:
                cont_padding = 0
            level += increase[i]

        if not free:
            return self._wrap_fortran(new_code)
        return new_code


def fcode(expr, assign_to=None, **settings):
    """Converts an expr to a string of fortran code

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
        A dictionary where keys are ``FunctionClass`` instances and values are
        their string representations. Alternatively, the dictionary value can
        be a list of tuples i.e. [(argument_test, cfunction_string)]. See below
        for examples.
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
    source_format : optional
        The source format can be either 'fixed' or 'free'. [default='fixed']
    standard : integer, optional
        The Fortran standard to be followed. This is specified as an integer.
        Acceptable standards are 66, 77, 90, 95, 2003, and 2008. Default is 77.
        Note that currently the only distinction internally is between
        standards before 95, and those 95 and after. This may change later as
        more features are added.

    Examples
    ========

    >>> fcode((2*x)**Rational(7, 2))
    '      8*sqrt(2.0d0)*x**(7.0d0/2.0d0)'
    >>> fcode(sin(x), assign_to='s')
    '      s = sin(x)'

    Custom printing can be defined for certain types by passing a dictionary of
    "type" : "function" to the ``user_functions`` kwarg. Alternatively, the
    dictionary value can be a list of tuples i.e. [(argument_test,
    cfunction_string)].

    >>> custom_functions = {'ceiling': 'CEIL',
    ...                     'floor': [(lambda x: not x.is_integer, 'FLOOR1'),
    ...                               (lambda x: x.is_integer, 'FLOOR2')]}
    >>> fcode(floor(x) + ceiling(x), user_functions=custom_functions)
    '      CEIL(x) + FLOOR1(x)'

    ``Piecewise`` expressions are converted into conditionals. If an
    ``assign_to`` variable is provided an if statement is created, otherwise
    the ternary operator is used. Note that if the ``Piecewise`` lacks a
    default term, represented by ``(expr, True)`` then an error will be thrown.
    This is to prevent generating an expression that may not evaluate to
    anything.

    >>> expr = Piecewise((x + 1, x > 0), (x, True))
    >>> print(fcode(expr, y))
          if (x > 0) then
             y = x + 1
          else
             y = x
          end if

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
    >>> fcode(e.rhs, assign_to=e.lhs, contract=False)
    '      Dy(i) = (y(i + 1) - y(i))/(t(i + 1) - t(i))'

    Matrices are also supported, but a ``MatrixSymbol`` of the same dimensions
    must be provided to ``assign_to``. Note that any expression that can be
    generated normally can also exist inside a Matrix:

    >>> mat = Matrix([x**2, Piecewise((x + 1, x > 0), (x, True)), sin(x)])
    >>> A = MatrixSymbol('A', 3, 1)
    >>> print(fcode(mat, A))
          A(1, 1) = x**2
             if (x > 0) then
          A(2, 1) = x + 1
             else
          A(2, 1) = x
             end if
          A(3, 1) = sin(x)

    """
    return FCodePrinter(settings).doprint(expr, assign_to)
