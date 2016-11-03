import pytest

from diofant.core import Symbol, Function, Float, Rational, Integer, I, Mul, Pow
from diofant.functions import exp, sin
from diofant.logic import And
from diofant.series import Limit

from diofant.parsing.sympy_parser import (parse_expr, standard_transformations,
                                          rationalize, TokenError, split_symbols,
                                          implicit_multiplication)


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
        "Float(Integer(3).n(n=3))": 3.00,
        'factorint(12, visual=True)': Mul(
            Pow(2, 2, evaluate=False),
            Pow(3, 1, evaluate=False),
            evaluate=False),
        'Limit(sin(x), x, 0, dir="-")': Limit(sin(x), x, 0, dir='-'),

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
    x = Symbol('x')
    e = '2*(x+1)'
    assert parse_expr(e, evaluate=0) == parse_expr(e, evaluate=False)


def test_split_symbols():
    transformations = standard_transformations + \
        (split_symbols, implicit_multiplication,)
    x = Symbol('x')
    y = Symbol('y')
    xy = Symbol('xy')

    assert parse_expr("xy") == xy
    assert parse_expr("xy", transformations=transformations) == x*y


def test_split_symbols_function():
    transformations = standard_transformations + \
        (split_symbols, implicit_multiplication,)
    x = Symbol('x')
    y = Symbol('y')
    a = Symbol('a')
    f = Function('f')

    assert parse_expr("ay(x+1)", transformations=transformations) == a*y*(x + 1)
    assert parse_expr("af(x+1)", transformations=transformations,
                      local_dict={'f': f}) == a*f(x + 1)


def test_match_parentheses_implicit_multiplication():
    transformations = standard_transformations + \
        (implicit_multiplication,)
    pytest.raises(TokenError, lambda: parse_expr('(1,2),(3,4]', transformations=transformations))
