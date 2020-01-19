"""Printing subsystem"""

from .ccode import ccode
from .dot import dotprint
from .fcode import fcode
from .latex import latex
from .mathematica import mathematica_code
from .mathml import mathml
from .octave import octave_code
from .pretty import pprint, pprint_use_unicode, pretty, pretty_print
from .python import python
from .repr import srepr
from .str import StrPrinter, sstr, sstrrepr


del str  # or this hide the str function
del repr  # or this hide the repr function


__all__ = ('ccode', 'dotprint', 'fcode', 'latex', 'mathematica_code',
           'mathml', 'octave_code', 'pprint', 'pprint_use_unicode',
           'pretty', 'pretty_print', 'python', 'srepr', 'sstr',
           'StrPrinter', 'sstrrepr')
