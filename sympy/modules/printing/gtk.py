"""Printing using GtkMathView"""

from sympy.core import Basic
from sympy.modules.mathml import c2p
import tempfile
import os

def print_gtk(x):
    """Print to Gtkmathview, a gtk widget capable of rendering MathML.
    Needs libgtkmathview-bin"""
    
    assert isinstance(x, Basic)
    
    tmp = tempfile.mktemp() # create a temp file to store the result
    file = open(tmp, 'wb')
    
    file.write( c2p(x.mathml, simple=True) )
    file.close()
    
    os.system("mathmlviewer " + tmp)
    
