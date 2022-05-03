import re

from ..concrete import Sum, product
from ..core import sympify
from ..functions import cos, sin


class MaximaHelpers:
    """Class for maxima parsing helpers."""

    def maxima_expand(self):
        return self.expand()

    def maxima_float(self):
        return self.evalf()

    def maxima_trigexpand(self):
        return self.expand(trig=True)

    def maxima_sum(self, var, low, high):
        return Sum(self, (var, low, high)).doit()

    def maxima_product(self, var, low, high):
        return product(self, (var, low, high))

    def maxima_csc(self):
        return 1/sin(self)

    def maxima_sec(self):
        return 1/cos(self)


sub_dict = {
    'pi': re.compile('%pi'),
    'E': re.compile('%e'),
    'I': re.compile('%i'),
    '**': re.compile(r'\^'),
    'oo': re.compile(r'\binf\b'),
    '-oo': re.compile(r'\bminf\b'),
    '1': re.compile(r'\bminus\b'),
    'maxima_expand': re.compile(r'\bexpand\b'),
    'maxima_float': re.compile(r'\bfloat\b'),
    'maxima_trigexpand': re.compile(r'\btrigexpand'),
    'maxima_sum': re.compile(r'\bsum\b'),
    'maxima_product': re.compile(r'\bproduct\b'),
    'cancel': re.compile(r'\bratsimp\b'),
    'maxima_csc': re.compile(r'\bcsc\b'),
    'maxima_sec': re.compile(r'\bsec\b')
}

var_name = re.compile(r'^\s*(\w+)\s*:')


def parse_maxima(str, globals=None, name_dict={}):
    str = str.strip()
    str = str.rstrip('; ')

    for k, v in sub_dict.items():
        str = v.sub(k, str)

    assign_var = None
    var_match = var_name.search(str)
    if var_match:
        assign_var = var_match.group(1)
        str = str[var_match.end():].strip()

    dct = MaximaHelpers.__dict__.copy()
    dct.update(name_dict)
    obj = sympify(str, locals=dct)

    if assign_var and globals:
        globals[assign_var] = obj

    return obj
