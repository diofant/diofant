from diofant import (ITE, And, Basic, Derivative, Eq, Equivalent, I, Implies,
                     Integer, Integral, MatrixSymbol, Nand, Nor, Not, Or,
                     Rational, Symbol, Tuple, Xor, cos, count_ops, exp, pi,
                     powsimp, sin, symbols)
from diofant.abc import a, x, y, z


__all__ = ()


def test_count_ops_non_visual():
    def count(val):
        return count_ops(val, visual=False)
    assert count(x) == 0
    assert count(x) is not Integer(0)
    assert count(x + y) == 1
    assert count(x + y) is not Integer(1)
    assert count(x + y*x + 2*y) == 4
    assert count({x + y: x}) == 1
    assert count({x + y: 2 + x}) is not Integer(1)
    assert count(Or(x, y)) == 1
    assert count(And(x, y)) == 1
    assert count(Not(x)) == 1
    assert count(Nor(x, y)) == 2
    assert count(Nand(x, y)) == 2
    assert count(Xor(x, y)) == 1
    assert count(Implies(x, y)) == 1
    assert count(Equivalent(x, y)) == 1
    assert count(ITE(x, y, z)) == 1
    assert count(ITE(True, x, y)) == 0


def test_count_ops_visual():
    ADD, MUL, POW, SIN, COS, AND, D, G = symbols(
        'Add Mul Pow sin cos And Derivative Integral'.upper())
    DIV, SUB, NEG = symbols('DIV SUB NEG')
    NOT, OR, AND, XOR, IMPLIES, EQUIVALENT, _ITE, BASIC, TUPLE = symbols(
        'Not Or And Xor Implies Equivalent ITE Basic Tuple'.upper())

    def count(val):
        return count_ops(val, visual=True)

    assert count(7) is Integer(0)
    assert count(-1) == NEG
    assert count(-2) == NEG
    assert count(Rational(2, 3)) == DIV
    assert count(pi/3) == DIV
    assert count(-pi/3) == DIV + NEG
    assert count(I - 1) == SUB
    assert count(1 - I) == SUB
    assert count(1 - 2*I) == SUB + MUL

    assert count(x) is Integer(0)
    assert count(-x) == NEG
    assert count(-2*x/3) == NEG + DIV + MUL
    assert count(1/x) == DIV
    assert count(1/(x*y)) == DIV + MUL
    assert count(-1/x) == NEG + DIV
    assert count(-2/x) == NEG + DIV
    assert count(x/y) == DIV
    assert count(-x/y) == NEG + DIV

    assert count(x**2) == POW
    assert count(-x**2) == POW + NEG
    assert count(-2*x**2) == POW + MUL + NEG

    assert count(x + pi/3) == ADD + DIV
    assert count(x + Rational(1, 3)) == ADD + DIV
    assert count(x + y) == ADD
    assert count(x - y) == SUB
    assert count(y - x) == SUB
    assert count(-1/(x - y)) == DIV + NEG + SUB
    assert count(-1/(y - x)) == DIV + NEG + SUB
    assert count(1 + x**y) == ADD + POW
    assert count(1 + x + y) == 2*ADD
    assert count(1 + x + y + z) == 3*ADD
    assert count(1 + x**y + 2*x*y + y**2) == 3*ADD + 2*POW + 2*MUL
    assert count(2*z + y + x + 1) == 3*ADD + MUL
    assert count(2*z + y**17 + x + 1) == 3*ADD + MUL + POW
    assert count(2*z + y**17 + x + sin(x)) == 3*ADD + POW + MUL + SIN
    assert count(2*z + y**17 + x + sin(x**2)) == 3*ADD + MUL + 2*POW + SIN
    assert count(2*z + y**17 + x + sin(
        x**2) + exp(cos(x))) == 4*ADD + MUL + 3*POW + COS + SIN

    assert count(Derivative(x, x)) == D
    assert count(Integral(x, x) + 2*x/(1 + x)) == G + DIV + MUL + 2*ADD
    assert count(Basic()) is Integer(0)

    assert count({x + 1: sin(x)}) == ADD + SIN
    assert count([x + 1, sin(x) + y, None]) == ADD + SIN + ADD
    assert count({x + 1: sin(x), y: cos(x) + 1}) == SIN + COS + 2*ADD
    assert count({}) is Integer(0)
    assert count([x + 1, sin(x)*y, None]) == SIN + ADD + MUL
    assert count([]) is Integer(0)

    assert count(Basic()) == 0
    assert count(Basic(Basic(), Basic(x, x + y))) == ADD + 2*BASIC
    assert count(Basic(x, x + y)) == ADD + BASIC
    assert count(Or(x, y)) == OR
    assert count(And(x, y)) == AND
    assert count(And(x**y, z)) == AND + POW
    assert count(Or(x, Or(y, And(z, a)))) == AND + OR
    assert count(Nor(x, y)) == NOT + OR
    assert count(Nand(x, y)) == NOT + AND
    assert count(Xor(x, y)) == XOR
    assert count(Implies(x, y)) == IMPLIES
    assert count(Equivalent(x, y)) == EQUIVALENT
    assert count(ITE(x, y, z)) == _ITE
    assert count([Or(x, y), And(x, y), Basic(x + y)]) == ADD + AND + BASIC + OR

    assert count(Basic(Tuple(x))) == BASIC + TUPLE
    # It checks that TUPLE is counted as an operation.

    assert count(Eq(x + y, 2)) == ADD


def test_sympyissue_9324():
    def count(val):
        return count_ops(val, visual=False)

    M = MatrixSymbol('M', 10, 10)
    assert count(M[0, 0]) == 0
    assert count(2 * M[0, 0] + M[5, 7]) == 2
    P = MatrixSymbol('P', 3, 3)
    Q = MatrixSymbol('Q', 3, 3)
    assert count(P + Q) == 9
    expr = powsimp(M, deep=True)
    assert expr == M
    assert expr.name == 'M'
    m = Symbol('m', integer=True)
    n = Symbol('n', integer=True)
    M = MatrixSymbol('M', m + n, m * m)
    assert count(M[0, 1]) == 0
