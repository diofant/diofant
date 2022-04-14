import pytest

from diofant import (Abs, And, Float, Function, I, Integer, Limit, Mul, Pow,
                     Rational, Symbol, exp, sin)
from diofant.parsing.sympy_parser import (TokenError, convert_xor,
                                          function_exponentiation,
                                          implicit_multiplication, parse_expr,
                                          rationalize, split_symbols,
                                          standard_transformations)


__all__ = ()


def test_diofant_parser():
    x = Symbol('x')
    inputs = {
        '2*x': 2 * x,
        '3.00': Float(3),
        '22/7': Rational(22, 7),
        '2+3j': 2 + 3*I,
        'exp(x)': exp(x),
        '-(2)': -Integer(2),
        '[-1, -2, 3]': [Integer(-1), Integer(-2), Integer(3)],
        'Symbol("x").free_symbols': x.free_symbols,
        'Float(Integer(3).evalf(3))': 3.00,
        'factorint(12, visual=True)': Mul(
            Pow(2, 2, evaluate=False),
            Pow(3, 1, evaluate=False),
            evaluate=False),
        'Limit(sin(x), x, 0, dir=1)': Limit(sin(x), x, 0, dir=1),

    }
    for text, result in inputs.items():
        assert parse_expr(text) == result


def test_rationalize():
    inputs = {
        '0.123': Rational(123, 1000)
    }
    transformations = standard_transformations + (rationalize,)
    for text, result in inputs.items():
        assert parse_expr(text, transformations=transformations) == result


def test_local_dict():
    local_dict = {
        'my_function': lambda x: x + 2
    }
    inputs = {
        'my_function(2)': Integer(4)
    }
    for text, result in inputs.items():
        assert parse_expr(text, local_dict=local_dict) == result


def test_local_dict_symbol_to_fcn():
    x = Symbol('x')
    d = {'foo': Function('bar')}
    assert parse_expr('foo(x)', local_dict=d) == d['foo'](x)
    # XXX: bit odd, but would be error if parser left the Symbol
    d = {'foo': Symbol('baz')}
    assert parse_expr('foo(x)', local_dict=d) == Function('baz')(x)


def test_global_dict():
    global_dict = {
        'Symbol': Symbol
    }
    inputs = {
        'Q & S': And(Symbol('Q'), Symbol('S'))
    }
    for text, result in inputs.items():
        assert parse_expr(text, global_dict=global_dict) == result


def test_sympyissue_2515():
    pytest.raises(TokenError, lambda: parse_expr('(()'))
    pytest.raises(TokenError, lambda: parse_expr('"""'))


def test_sympyissue_7663():
    e = '2*(x+1)'
    assert parse_expr(e, evaluate=0) == parse_expr(e, evaluate=False)


def test_split_symbols():
    transformations = standard_transformations + \
        (split_symbols, implicit_multiplication)
    x = Symbol('x')
    y = Symbol('y')
    xy = Symbol('xy')

    assert parse_expr('xy') == xy
    assert parse_expr('xy', transformations=transformations) == x*y


def test_split_symbols_function():
    transformations = standard_transformations + \
        (split_symbols, implicit_multiplication)
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    f = Function('f')

    assert parse_expr('ay(x+1)', transformations=transformations) == a*y*(x + 1)
    assert parse_expr('af(x+1)', transformations=transformations,
                      local_dict={'f': f}) == a*f(x + 1)


def test_functional_exponent():
    t = standard_transformations + (convert_xor, function_exponentiation)
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    yfcn = Function('y')
    assert parse_expr('sin^2(x)', transformations=t) == (sin(x))**2
    assert parse_expr('sin^y(x)', transformations=t) == (sin(x))**y
    assert parse_expr('exp^y(x)', transformations=t) == (exp(x))**y
    assert parse_expr('E^y(x)', transformations=t) == exp(yfcn(x))
    assert parse_expr('a^y(x)', transformations=t) == a**(yfcn(x))


def test_match_parentheses_implicit_multiplication():
    transformations = standard_transformations + (implicit_multiplication,)
    pytest.raises(TokenError, lambda: parse_expr('(1,2),(3,4]', transformations=transformations))


def test_sympyissue_22020():
    x = parse_expr('log((2*V/3-V)/C)/-(R+r)*C')
    y = parse_expr('log((2*V/3-V)/C)/-(R+r)*2')

    assert x.equals(y) is False


def test_builtins():
    cases = [('abs(x)', 'Abs(x)'),
             ('max(x, y)', 'Max(x, y)'),
             ('min(x, y)', 'Min(x, y)'),
             ('pow(x, y)', 'Pow(x, y)')]
    for built_in_func_call, sympy_func_call in cases:
        assert parse_expr(built_in_func_call) == parse_expr(sympy_func_call)
    assert parse_expr('pow(38, -1)') == Rational(1, 38)
    # issue sympy/sympy#22322
    assert parse_expr('abs(-42)', evaluate=False) == Abs(-42, evaluate=False)
