"""A module providing information about the necessity of brackets"""

from ..core.function import _coeff_isneg


#: Default precedence values for some basic types.
PRECEDENCE = {
    "Lambda": 1,
    "Xor": 10,
    "Or": 20,
    "And": 30,
    "Relational": 35,
    "Add": 40,
    "Mul": 50,
    "Mod": 50,
    "Pow": 60,
    "Func": 70,
    "Not": 100,
    "Atom": 1000,
    "BitwiseOr": 36,
    "BitwiseAnd": 38,
}

#: A dictionary assigning precedence values to certain classes. These values
#: are treated like they were inherited, so not every single class has to
#: be named here.
PRECEDENCE_VALUES = {
    "Equivalent": PRECEDENCE["Xor"],
    "Xor": PRECEDENCE["Xor"],
    "Implies": PRECEDENCE["Xor"],
    "Or": PRECEDENCE["Or"],
    "And": PRECEDENCE["And"],
    "Add": PRECEDENCE["Add"],
    "Mod": PRECEDENCE["Mod"],
    "Pow": PRECEDENCE["Pow"],
    "Relational": PRECEDENCE["Relational"],
    "Sub": PRECEDENCE["Add"],
    "Not": PRECEDENCE["Not"],
    "factorial": PRECEDENCE["Func"],
    "factorial2": PRECEDENCE["Func"],
    "Function": PRECEDENCE["Func"],
    "NegativeInfinity": PRECEDENCE["Add"],
    "MatAdd": PRECEDENCE["Add"],
    "MatMul": PRECEDENCE["Mul"],
    "MatPow": PRECEDENCE["Pow"],
    "HadamardProduct": PRECEDENCE["Mul"],
    "Fraction": PRECEDENCE["Atom"],
}

# Sometimes it's not enough to assign a fixed precedence value to a
# class. Then a function can be inserted in this dictionary that takes
# an instance of this class as argument and returns the appropriate
# precedence value.

# Precedence functions


def precedence_Mul(item):
    if _coeff_isneg(item):
        return PRECEDENCE["Add"]
    return PRECEDENCE["Mul"]


def precedence_Rational(item):
    if item.numerator < 0:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Mul"]


def precedence_Integer(item):
    if item.numerator < 0:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Atom"]


def precedence_Float(item):
    if item < 0:
        return PRECEDENCE["Add"]
    return PRECEDENCE["Atom"]


def precedence_PolyElement(item):
    if item.is_generator:
        return PRECEDENCE["Atom"]
    elif item.is_ground:
        return precedence(item.coeff(1))
    elif item.is_term:
        return PRECEDENCE["Mul"]
    else:
        return PRECEDENCE["Add"]


def precedence_FracElement(item):
    if item.denominator == 1:
        return precedence_PolyElement(item.numerator)
    else:
        return PRECEDENCE["Mul"]


#: Sometimes it's not enough to assign a fixed precedence value to a class. Then
#: a function can be inserted in this dictionary that takes an instance of this
#: class as argument and returns the appropriate precedence value.
PRECEDENCE_FUNCTIONS = {
    "Integer": precedence_Integer,
    "Mul": precedence_Mul,
    "Rational": precedence_Rational,
    "Float": precedence_Float,
    "PolyElement": precedence_PolyElement,
    "FracElement": precedence_FracElement,
}


def precedence(item):
    """
    Returns the precedence of a given object.
    """
    mro = item.__class__.__mro__
    for i in mro:
        n = i.__name__
        if n in PRECEDENCE_FUNCTIONS:
            return PRECEDENCE_FUNCTIONS[n](item)
        elif n in PRECEDENCE_VALUES:
            return PRECEDENCE_VALUES[n]
    return PRECEDENCE["Atom"]
