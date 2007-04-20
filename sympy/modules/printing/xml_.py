from xml.dom.minidom import parseString
from sympy.core import Basic
from sympy.modules.mathml import add_mathml_headers

def print_xml(expr):
    """
    Print's a pretty representation of the mathml code for expr
    
    >>> from sympy import *
    >>> from sympy.modules.printing import print_xml
    >>> x = Symbol('x')
    >>> print_xml(x+1) #doctest: +NORMALIZE_WHITESPACE
    <?xml version="1.0" ?>
    <math xmlns:mml="http://www.w3.org/1998/Math/MathML" xmlns:sympy="http://www.w3schools.com/furniture" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.w3.org/1998/Math/MathML         http://www.w3.org/Math/XMLSchema/mathml2/mathml2.xsd">
        <apply>
                <plus/>
                <cn sympy:assumptions="is_real:True;is_commutative:True;">
                         1 
                </cn>
                <ci sympy:assumptions="is_commutative:True;">
                         x 
                </ci>
        </apply>
    </math>
    """
    expr = Basic.sympify(expr)
    s = add_mathml_headers(expr.mathml)
    p = parseString(s)
    print p.toprettyxml()
