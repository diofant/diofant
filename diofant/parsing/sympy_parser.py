"""Transform a string with Python-like source code into Diofant expression."""

import ast
import unicodedata
from io import BytesIO
from keyword import iskeyword
from tokenize import (ENDMARKER, NAME, NEWLINE, NUMBER, OP, STRING, TokenError,
                      tokenize, untokenize)

from ..core import Basic, Symbol


def _token_splittable(token):
    """
    Predicate for whether a token name can be split into multiple tokens.

    A token is splittable if it does not contain an underscore character and
    it is not the name of a Greek letter. This is used to implicitly convert
    expressions like 'xyz' into 'x*y*z'.

    """
    if '_' in token:
        return False
    else:
        try:
            return not unicodedata.lookup('GREEK SMALL LETTER ' + token)
        except KeyError:
            pass
    if len(token) > 1:
        return True
    return False


def _token_callable(token, local_dict, global_dict, nextToken=None):
    """
    Predicate for whether a token name represents a callable function.

    Essentially wraps ``callable``, but looks up the token name in the
    locals and globals.

    """
    func = local_dict.get(token[1])
    if not func:
        func = global_dict.get(token[1])
    return callable(func) and not isinstance(func, Symbol)


class AppliedFunction:
    """
    A group of tokens representing a function and its arguments.

    `exponent` is for handling the shorthand sin^2, ln^2, etc.

    """

    def __init__(self, function, args, exponent=[]):
        self.function = function
        self.args = args
        self.exponent = exponent
        self.items = ['function', 'args', 'exponent']

    def expand(self):
        """Return a list of tokens representing the function."""
        result = []
        result.append(self.function)
        result.extend(self.args)
        return result

    def __getitem__(self, index):
        return getattr(self, self.items[index])


class ParenthesisGroup(list):
    """List of tokens representing an expression in parentheses."""


def _flatten(result):
    result2 = []
    for tok in result:
        if isinstance(tok, AppliedFunction):
            result2.extend(tok.expand())
        else:
            result2.append(tok)
    return result2


def _group_parentheses(recursor):
    def _inner(tokens, local_dict, global_dict):
        """Group tokens between parentheses with ParenthesisGroup.

        Also processes those tokens recursively.

        """
        result = []
        stacks = []
        stacklevel = 0
        for token in tokens:
            if token[0] == OP:
                if token[1] == '(':
                    stacks.append(ParenthesisGroup([]))
                    stacklevel += 1
                elif token[1] == ')':
                    stacks[-1].append(token)
                    stack = stacks.pop()

                    if len(stacks) > 0:
                        # We don't recurse here since the upper-level stack
                        # would reprocess these tokens
                        stacks[-1].extend(stack)
                    else:
                        # Recurse here to handle nested parentheses
                        # Strip off the outer parentheses to avoid an infinite loop
                        inner = stack[1:-1]
                        inner = recursor(inner,
                                         local_dict,
                                         global_dict)
                        parenGroup = [stack[0]] + inner + [stack[-1]]
                        result.append(ParenthesisGroup(parenGroup))
                    stacklevel -= 1
                    continue
            if stacklevel:
                stacks[-1].append(token)
            else:
                result.append(token)
        if stacklevel:
            raise TokenError('Mismatched parentheses')
        return result
    return _inner


def _apply_functions(tokens, local_dict, global_dict):
    """Convert a NAME token + ParenthesisGroup into an AppliedFunction.

    Note that ParenthesisGroups, if not applied to any function, are
    converted back into lists of tokens.

    """
    result = []
    symbol = None
    for tok in tokens:
        if tok[0] == NAME:
            symbol = tok
            result.append(tok)
        elif isinstance(tok, ParenthesisGroup):
            if symbol and _token_callable(symbol, local_dict, global_dict):
                result[-1] = AppliedFunction(symbol, tok)
                symbol = None
            else:
                result.extend(tok)
        else:
            symbol = None
            result.append(tok)
    return result


def _implicit_multiplication(tokens, local_dict, global_dict):
    """Implicitly adds '*' tokens.

    Cases:

    - Two AppliedFunctions next to each other ("sin(x)cos(x)")

    - AppliedFunction next to an open parenthesis ("sin x (cos x + 1)")

    - A close parenthesis next to an AppliedFunction ("(x+2)sin x")

    - A close parenthesis next to an open parenthesis ("(x+2)(x+3)")

    - AppliedFunction next to an implicitly applied function ("sin(x)cos x")

    """
    result = []
    for tok, nextTok in zip(tokens, tokens[1:]):
        result.append(tok)
        if (isinstance(tok, AppliedFunction) and
                isinstance(nextTok, AppliedFunction)):
            result.append((OP, '*'))
        elif (isinstance(tok, AppliedFunction) and
              nextTok[0] == OP and nextTok[1] == '('):
            # Applied function followed by an open parenthesis
            if tok.function[1] == 'Function':
                result[-1].function = (result[-1].function[0], 'Symbol')
            result.append((OP, '*'))
        elif (tok[0] == OP and tok[1] == ')' and
              isinstance(nextTok, AppliedFunction)):
            # Close parenthesis followed by an applied function
            result.append((OP, '*'))
        elif (tok[0] == OP and tok[1] == ')' and
              nextTok[0] == NAME):
            # Close parenthesis followed by an implicitly applied function
            result.append((OP, '*'))
        elif (tok[0] == nextTok[0] == OP
              and tok[1] == ')' and nextTok[1] == '('):
            # Close parenthesis followed by an open parenthesis
            result.append((OP, '*'))
        elif (isinstance(tok, AppliedFunction) and nextTok[0] == NAME):
            # Applied function followed by implicitly applied function
            result.append((OP, '*'))
        elif (tok[0] == NAME and
              not _token_callable(tok, local_dict, global_dict) and
              nextTok[0] == OP and nextTok[1] == '('):
            # Constant followed by parenthesis
            result.append((OP, '*'))
        elif (tok[0] == NAME and
              not _token_callable(tok, local_dict, global_dict) and
              nextTok[0] == NAME and
              not _token_callable(nextTok, local_dict, global_dict)):
            # Constant followed by constant
            result.append((OP, '*'))
        elif (tok[0] == NAME and
              not _token_callable(tok, local_dict, global_dict) and
              (isinstance(nextTok, AppliedFunction) or nextTok[0] == NAME)):
            # Constant followed by (implicitly applied) function
            result.append((OP, '*'))
    result.append(tokens[-1])
    return result


def _implicit_application(tokens, local_dict, global_dict):
    """Adds parentheses as needed after functions."""
    result = []
    appendParen = 0  # number of closing parentheses to add

    # number of tokens to delay before adding a ')' (to
    # capture **, ^, etc.)
    skip = 0

    # skipping tokens before inserting parentheses to
    # work with function exponentiation
    exponentSkip = False

    for tok, nextTok in zip(tokens, tokens[1:]):
        result.append(tok)
        if (tok[0] == NAME and
            nextTok[0] != OP and
                nextTok[0] not in (ENDMARKER, NEWLINE)):
            if _token_callable(tok, local_dict, global_dict, nextTok):
                result.append((OP, '('))
                appendParen += 1
            else:
                raise NotImplementedError
        # name followed by exponent - function exponentiation
        elif (tok[0] == NAME and nextTok[0] == OP and nextTok[1] == '**'):
            if _token_callable(tok, local_dict, global_dict):
                exponentSkip = True
        elif exponentSkip:
            # if the last token added was an applied function (i.e. the
            # power of the function exponent) OR a multiplication (as
            # implicit multiplication would have added an extraneous
            # multiplication)
            if (isinstance(tok, AppliedFunction)
                    or (tok[0] == OP and tok[1] == '*')):
                # don't add anything if the next token is a multiplication
                # or if there's already a parenthesis (if parenthesis, still
                # stop skipping tokens)
                if not (nextTok[0] == OP and nextTok[1] == '*'):
                    if not(nextTok[0] == OP and nextTok[1] == '('):
                        result.append((OP, '('))
                        appendParen += 1
                    exponentSkip = False
        elif appendParen:
            if nextTok[0] == OP and nextTok[1] in ('^', '**', '*'):
                skip = 1
                continue
            if skip:
                skip -= 1
                continue
            result.append((OP, ')'))
            appendParen -= 1

    result.append(tokens[-1])

    if appendParen:
        result.extend([(OP, ')')] * appendParen)
    return result


def function_exponentiation(tokens, local_dict, global_dict):
    """Allows functions to be exponentiated, e.g. ``cos**2(x)``.

    Examples
    ========

    >>> transformations = standard_transformations + (function_exponentiation,)
    >>> parse_expr('sin**4(x)', transformations=transformations)
    sin(x)**4

    """
    result = []
    exponent = []
    consuming_exponent = False
    level = 0
    for tok, nextTok in zip(tokens, tokens[1:]):
        if tok[0] == NAME and nextTok[0] == OP and nextTok[1] == '**':
            if _token_callable(tok, local_dict, global_dict):
                consuming_exponent = True
        elif consuming_exponent:
            if tok[0] == NAME and tok[1] == 'Function':
                tok = (NAME, 'Symbol')
            # only want to stop after hitting )
            if tok[0] == nextTok[0] == OP and tok[1] == ')' and nextTok[1] == '(':
                consuming_exponent = False
            # if implicit multiplication was used, we may have )*( instead
            if tok[0] == nextTok[0] == OP and tok[1] == '*' and nextTok[1] == '(':
                consuming_exponent = False
            else:
                exponent.append(tok)
            continue
        elif exponent and not consuming_exponent:
            if tok[0] == OP:
                if tok[1] == '(':
                    level += 1
                elif tok[1] == ')':
                    level -= 1
            if level == 0:
                result.append(tok)
                result.extend(exponent)
                exponent = []
                continue
        result.append(tok)
    result.append(tokens[-1])
    if exponent:
        raise NotImplementedError
    return result


def split_symbols_custom(predicate):
    """Creates a transformation that splits symbol names.

    ``predicate`` should return True if the symbol name is to be split.

    For instance, to retain the default behavior but avoid splitting certain
    symbol names, a predicate like this would work:


    >>> def can_split(symbol):
    ...     if symbol not in ('list', 'of', 'unsplittable', 'names'):
    ...         return _token_splittable(symbol)
    ...     return False
    ...
    >>> transformation = split_symbols_custom(can_split)
    >>> parse_expr('unsplittable', transformations=standard_transformations +
    ...            (transformation, implicit_multiplication))
    unsplittable

    """
    def _split_symbols(tokens, local_dict, global_dict):
        result = []
        split = False
        split_previous = False
        for tok in tokens:
            if split_previous:
                # throw out closing parenthesis of Symbol that was split
                split_previous = False
                continue
            split_previous = False
            if tok[0] == NAME and tok[1] in ['Symbol', 'Function']:
                split = True
            elif split and tok[0] == NAME:
                symbol = tok[1][1:-1]
                if predicate(symbol):
                    tok_type = result[-2][1]  # Symbol or Function
                    del result[-2:]  # Get rid of the call to Symbol
                    for char in symbol[:-1]:
                        if char in local_dict or char in global_dict:
                            result.extend([(NAME, f'{char}')])
                        else:
                            result.extend([(NAME, 'Symbol'), (OP, '('),
                                           (NAME, f"'{char}'"), (OP, ')')])
                    char = symbol[-1]
                    if char in local_dict or char in global_dict:
                        result.extend([(NAME, f'{char}')])
                    else:
                        result.extend([(NAME, tok_type), (OP, '('),
                                       (NAME, f"'{char}'"), (OP, ')')])

                    # Set split_previous=True so will skip
                    # the closing parenthesis of the original Symbol
                    split = False
                    split_previous = True
                    continue
                else:
                    split = False
            result.append(tok)
        return result
    return _split_symbols


#: Splits symbol names for implicit multiplication.
#:
#: Intended to let expressions like ``xyz`` be parsed as ``x*y*z``. Does not
#: split Greek character names, so ``theta`` will *not* become
#: ``t*h*e*t*a``. Generally this should be used with
#: ``implicit_multiplication``.
split_symbols = split_symbols_custom(_token_splittable)


def implicit_multiplication(result, local_dict, global_dict):
    """Makes the multiplication operator optional in most cases.

    Use this before :func:`~diofant.parsing.sympy_parser.implicit_application`, otherwise expressions like
    ``sin 2x`` will be parsed as ``x * sin(2)`` rather than ``sin(2*x)``.

    Examples
    ========

    >>> transformations = standard_transformations + (implicit_multiplication,)
    >>> parse_expr('3 x y', transformations=transformations)
    3*x*y

    """
    # These are interdependent steps, so we don't expose them separately
    for step in (_group_parentheses(implicit_multiplication),
                 _apply_functions,
                 _implicit_multiplication):
        result = step(result, local_dict, global_dict)

    result = _flatten(result)
    return result


def implicit_application(result, local_dict, global_dict):
    """Makes parentheses optional in some cases for function calls.

    Use this after :func:`~diofant.parsing.sympy_parser.implicit_multiplication`, otherwise expressions
    like ``sin 2x`` will be parsed as ``x * sin(2)`` rather than
    ``sin(2*x)``.

    Examples
    ========

    >>> transformations = standard_transformations + (implicit_application,)
    >>> parse_expr('cot z + csc z', transformations=transformations)
    cot(z) + csc(z)

    """
    for step in (_group_parentheses(implicit_application),
                 _apply_functions,
                 _implicit_application,):
        result = step(result, local_dict, global_dict)

    result = _flatten(result)
    return result


def implicit_multiplication_application(result, local_dict, global_dict):
    """Allows a slightly relaxed syntax.

    - Parentheses for single-argument method calls are optional.

    - Multiplication is implicit.

    - Symbol names can be split (i.e. spaces are not needed between
      symbols).

    - Functions can be exponentiated.

    Examples
    ========

    >>> parse_expr('10sin**2 x**2 + 3xyz + tan theta',
    ...            transformations=(standard_transformations +
    ...                             (implicit_multiplication_application,)))
    3*x*y*z + 10*sin(x**2)**2 + tan(theta)

    """
    for step in (split_symbols, implicit_multiplication,
                 implicit_application, function_exponentiation):
        result = step(result, local_dict, global_dict)

    return result


def auto_symbol(tokens, local_dict, global_dict):
    """Inserts calls to ``Symbol``/``Function`` for undefined variables."""
    result = []
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == NAME:
            name = tokVal

            if (name in ['True', 'False', 'None']
                or iskeyword(name)
                # Don't convert attribute access
                or (prevTok[0] == OP and prevTok[1] == '.')
                # Don't convert keyword arguments
                or (prevTok[0] == OP and prevTok[1] in ('(', ',')
                    and nextTokNum == OP and nextTokVal == '=')):
                result.append((NAME, name))
                continue
            elif name in local_dict:
                if isinstance(local_dict[name], Symbol) and nextTokVal == '(':
                    result.extend([(NAME, 'Function'),
                                   (OP, '('),
                                   (NAME, repr(str(local_dict[name]))),
                                   (OP, ')')])
                else:
                    result.append((NAME, name))
                continue
            elif name in global_dict:
                obj = global_dict[name]
                assert isinstance(obj, (Basic, type)) or callable(obj)
                result.append((NAME, name))
                continue

            result.extend([
                (NAME, 'Symbol' if nextTokVal != '(' else 'Function'),
                (OP, '('),
                (NAME, repr(str(name))),
                (OP, ')'),
            ])
        else:
            result.append((tokNum, tokVal))

        prevTok = (tokNum, tokVal)

    return result


def lambda_notation(tokens, local_dict, global_dict):
    """Substitutes "lambda" with its Sympy equivalent Lambda().
    However, the conversion doesn't take place if only "lambda"
    is passed because that is a syntax error.

    """
    result = []
    flag = False
    if tokens[0][1] == 'utf-8' and tokens[1] == (NAME, 'lambda'):
        tokens = tokens[1:]
    toknum, tokval = tokens[0]
    tokLen = len(tokens)
    if toknum == NAME and tokval == 'lambda':
        if tokLen == 2 or (tokLen == 3 and tokens[1][0] == NEWLINE):
            result.extend(tokens)
        elif tokLen > 2:
            result.extend([
                (NAME, 'Lambda'),
                (OP, '('),
                (OP, '('),
                (OP, ')'),
                (OP, ')'),
            ])
            for tokNum, tokVal in tokens[1:]:
                if tokNum == OP and tokVal == ':':
                    tokVal = ','
                    flag = True
                if not flag and tokNum == OP and tokVal in ['*', '**']:
                    raise NotImplementedError('Starred arguments in lambda not supported')
                if flag:
                    result.insert(-1, (tokNum, tokVal))
                else:
                    result.insert(-2, (tokNum, tokVal))
        else:
            raise NotImplementedError
    else:
        result.extend(tokens)

    return result


def convert_xor(tokens, local_dict, global_dict):
    """Treats XOR, ``^``, as exponentiation, ``**``."""
    result = []
    for toknum, tokval in tokens:
        if toknum == OP:
            if tokval == '^':
                result.append((OP, '**'))
            else:
                result.append((toknum, tokval))
        else:
            result.append((toknum, tokval))

    return result


def auto_number(tokens, local_dict, global_dict):
    """Converts numeric literals to use Diofant equivalents.

    Complex numbers use ``I``; integer literals use ``Integer``, float
    literals use ``Float``, and repeating decimals use ``Rational``.

    """
    result = []

    for toknum, tokval in tokens:
        if toknum == NUMBER:
            number = tokval
            postfix = []

            if number.endswith('j') or number.endswith('J'):
                number = number[:-1]
                postfix = [(OP, '*'), (NAME, 'I')]

            if '.' in number or (('e' in number or 'E' in number) and
                                 not (number.startswith('0x') or number.startswith('0X'))):
                seq = [(NAME, 'Float'), (OP, '('),
                       (NUMBER, repr(str(number))), (OP, ')')]
            else:
                seq = [(NAME, 'Integer'), (OP, '('), (
                    NUMBER, number), (OP, ')')]

            result.extend(seq + postfix)
        else:
            result.append((toknum, tokval))

    return result


def rationalize(tokens, local_dict, global_dict):
    """Converts floats into ``Rational``. Run AFTER ``auto_number``."""
    result = []
    passed_float = False
    for toknum, tokval in tokens:
        if toknum == NAME:
            if tokval == 'Float':
                passed_float = True
                tokval = 'Rational'
            result.append((toknum, tokval))
        elif passed_float and toknum == NUMBER:
            passed_float = False
            result.append((STRING, tokval))
        else:
            result.append((toknum, tokval))

    return result


#: Standard transformations for :func:`~diofant.parsing.sympy_parser.parse_expr`.
#: Inserts calls to :class:`~diofant.core.symbol.Symbol`, :class:`~diofant.core.numbers.Integer`, and other Diofant datatypes.
standard_transformations = (lambda_notation, auto_symbol, auto_number)


def stringify_expr(s, local_dict, global_dict, transformations):
    """
    Converts the string ``s`` to Python code, in ``local_dict``

    Generally, ``parse_expr`` should be used.

    """
    tokens = []
    input_code = BytesIO(s.encode('utf-8').strip())
    for toknum, tokval, _, _, _ in tokenize(input_code.readline):
        tokens.append((toknum, tokval))

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    return untokenize(tokens)


def eval_expr(code, local_dict, global_dict):
    """
    Evaluate Python code generated by ``stringify_expr``.

    Generally, ``parse_expr`` should be used.

    """
    expr = eval(
        code, global_dict, local_dict)  # take local objects in preference

    return expr


def parse_expr(s, local_dict=None, transformations=standard_transformations,
               global_dict=None, evaluate=True):
    r"""Converts the string ``s`` to a Diofant expression, in ``local_dict``

    Parameters
    ==========

    s : str
        The string to parse.

    local_dict : dict, optional
        A dictionary of local variables to use when parsing.

    global_dict : dict, optional
        A dictionary of global variables. By default, this is initialized
        with ``from diofant import *``; provide this parameter to override
        this behavior (for instance, to parse ``"Q & S"``).

    transformations : tuple, optional
        A tuple of transformation functions used to modify the tokens of the
        parsed expression before evaluation. The default transformations
        convert numeric literals into their Diofant equivalents, convert
        undefined variables into Diofant symbols.

    evaluate : bool, optional
        When False, the order of the arguments will remain as they were in the
        string and automatic simplification that would normally occur is
        suppressed. (see examples)

    Examples
    ========

    >>> parse_expr('1/2')
    1/2
    >>> type(_)
    <class 'diofant.core.numbers.Half'>
    >>> transformations = (standard_transformations +
    ...                    (implicit_multiplication_application,))
    >>> parse_expr('2x', transformations=transformations)
    2*x

    When evaluate=False, some automatic simplifications will not occur:

    >>> parse_expr('2**3'), parse_expr('2**3', evaluate=False)
    (8, 2**3)

    In addition the order of the arguments will not be made canonical.
    This feature allows one to tell exactly how the expression was entered:

    >>> a = parse_expr('1 + x', evaluate=False)
    >>> b = parse_expr('x + 1', evaluate=0)
    >>> a == b
    False
    >>> a.args
    (1, x)
    >>> b.args
    (x, 1)

    See Also
    ========

    diofant.parsing.sympy_parser.stringify_expr
    diofant.parsing.sympy_parser.eval_expr
    diofant.parsing.sympy_parser.standard_transformations
    diofant.parsing.sympy_parser.implicit_multiplication_application

    """
    if local_dict is None:
        local_dict = {}

    if global_dict is None:
        global_dict = {}
        exec('from diofant import *', global_dict)

    code = stringify_expr(s, local_dict, global_dict, transformations)

    if not evaluate:
        code = compile(evaluateFalse(code), '<string>', 'eval')

    return eval_expr(code, local_dict, global_dict)


def evaluateFalse(s):
    """
    Replaces operators with the Diofant equivalent and sets evaluate=False.

    """
    node = ast.parse(s)
    node = EvaluateFalseTransformer().visit(node)
    # node is a Module, we want an Expression
    node = ast.Expression(node.body[0].value)

    return ast.fix_missing_locations(node)


class EvaluateFalseTransformer(ast.NodeTransformer):
    """Transformer class to hold evaluation."""

    operators = {
        ast.Add: 'Add',
        ast.Mult: 'Mul',
        ast.Pow: 'Pow',
        ast.Sub: 'Add',
        ast.Div: 'Mul',
        ast.BitOr: 'Or',
        ast.BitAnd: 'And',
        ast.BitXor: 'Not',
    }

    def flatten(self, args, func):
        result = []
        for arg in args:
            if isinstance(arg, ast.Call) and arg.func.id == func:
                result.extend(self.flatten(arg.args, func))
            else:
                result.append(arg)
        return result

    def visit_BinOp(self, node):
        if node.op.__class__ in self.operators:
            sympy_class = self.operators[node.op.__class__]
            right = self.visit(node.right)

            if isinstance(node.op, ast.Sub):
                right = ast.UnaryOp(op=ast.USub(), operand=right)
            elif isinstance(node.op, ast.Div):
                right = ast.Call(
                    func=ast.Name(id='Pow', ctx=ast.Load()),
                    args=[right, ast.UnaryOp(op=ast.USub(), operand=ast.Num(1))],
                    keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False, kind=None))],
                    starargs=None,
                    kwargs=None
                )

            new_node = ast.Call(
                func=ast.Name(id=sympy_class, ctx=ast.Load()),
                args=[self.visit(node.left), right],
                keywords=[ast.keyword(arg='evaluate', value=ast.Constant(value=False, kind=None))],
                starargs=None,
                kwargs=None
            )

            if sympy_class in ('Add', 'Mul'):
                # Denest Add or Mul as appropriate
                new_node.args = self.flatten(new_node.args, sympy_class)

            return new_node
        else:  # pragma: no cover
            return node
