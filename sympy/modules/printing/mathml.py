from sympy import Basic

def mathml(expr):
    """Returns the mathml representation of expr"""
    expr = Basic.sympify(expr)
    s = expr.__mathml__().toprettyxml()
    return s


def print_mathml(expr):
    """
    Print's a pretty representation of the mathml code for expr
    
    >>> from sympy import *
    >>> from sympy.modules.printing.mathml import print_mathml
    >>> x = Symbol('x')
    >>> print_mathml(x+1) #doctest: +NORMALIZE_WHITESPACE
    <apply>
        <plus/>
        <cn>
                1
        </cn>
        <ci>
                x
        </ci>
    </apply>

    """
    print mathml(expr)
