"""Symbolic primitives + unicode/ASCII abstraction for pretty.py"""

from __future__ import annotations

import typing
import unicodedata

from ...core.alphabets import greeks
from ..conventions import split_super_sub


# first, setup unicodedate environment
def U(name):
    """Unicode character by name or None if not found."""
    try:
        u = unicodedata.lookup(name)
    except KeyError:
        u = None

    return u


# prefix conventions when constructing tables
# L   - LATIN     i
# G   - GREEK     beta
# D   - DIGIT     0
# S   - SYMBOL    +


__all__ = ('greek_unicode', 'sub', 'sup', 'xsym', 'vobj', 'hobj',
           'pretty_symbol', 'annotated')


_use_unicode = False


def pretty_use_unicode(flag=None):
    """Set whether pretty-printer should use unicode by default."""
    global _use_unicode
    if flag is None:
        return _use_unicode

    use_unicode_prev = _use_unicode
    _use_unicode = flag
    return use_unicode_prev


# GREEK


def g(l):
    return U('GREEK SMALL LETTER %s' % l.upper())


def G(l):
    return U('GREEK CAPITAL LETTER %s' % l.upper())


greek_letters = list(greeks)  # make a copy
# deal with Unicode's funny spelling of lambda
greek_letters[greek_letters.index('lambda')] = 'lamda'

# {}  greek letter -> (g,G)
greek_unicode: dict[str, typing.Union[tuple, str]] = {l: (g(l), G(l)) for l in greek_letters}
greek_unicode = {L: g(L) for L in greek_letters}
greek_unicode.update((L[0].upper() + L[1:], G(L)) for L in greek_letters)

# aliases
greek_unicode['lambda'] = greek_unicode['lamda']
greek_unicode['Lambda'] = greek_unicode['Lamda']
greek_unicode['varsigma'] = '\N{GREEK SMALL LETTER FINAL SIGMA}'

digit_2txt = {
    '0':    'ZERO',
    '1':    'ONE',
    '2':    'TWO',
    '3':    'THREE',
    '4':    'FOUR',
    '5':    'FIVE',
    '6':    'SIX',
    '7':    'SEVEN',
    '8':    'EIGHT',
    '9':    'NINE',
}

symb_2txt = {
    '+':    'PLUS SIGN',
    '-':    'MINUS',
    '=':    'EQUALS SIGN',
    '(':    'LEFT PARENTHESIS',
    ')':    'RIGHT PARENTHESIS',
    '[':    'LEFT SQUARE BRACKET',
    ']':    'RIGHT SQUARE BRACKET',
    '{':    'LEFT CURLY BRACKET',
    '}':    'RIGHT CURLY BRACKET',

    # non-std
    '{}':   'CURLY BRACKET',
    'sum':  'SUMMATION',
    'int':  'INTEGRAL',
}

# SUBSCRIPT & SUPERSCRIPT


def LSUB(letter):
    return U('LATIN SUBSCRIPT SMALL LETTER %s' % letter.upper())


def GSUB(letter):
    return U('GREEK SUBSCRIPT SMALL LETTER %s' % letter.upper())


def DSUB(digit):
    return U(f'SUBSCRIPT {digit_2txt[digit]}')


def SSUB(symb):
    return U(f'SUBSCRIPT {symb_2txt[symb]}')


def LSUP(letter):
    return U('SUPERSCRIPT LATIN SMALL LETTER %s' % letter.upper())


def DSUP(digit):
    return U(f'SUPERSCRIPT {digit_2txt[digit]}')


def SSUP(symb):
    return U(f'SUPERSCRIPT {symb_2txt[symb]}')


sub = {}    # symb -> subscript symbol
sup = {}    # symb -> superscript symbol

# latin subscripts
for l in 'aeioruvxhklmnpst':
    sub[l] = LSUB(l)

for l in 'in':
    sup[l] = LSUP(l)

for gl in ['beta', 'gamma', 'rho', 'phi', 'chi']:
    sub[gl] = GSUB(gl)

for d in [str(i) for i in range(10)]:
    sub[d] = DSUB(d)
    sup[d] = DSUP(d)

for s in '+-=()':
    sub[s] = SSUB(s)
    sup[s] = SSUP(s)

# Variable modifiers
# TODO: Is it worth trying to handle faces with, e.g., 'MATHEMATICAL BOLD CAPITAL A'?
# TODO: Make brackets adjust to height of contents
modifier_dict = {
    # Accents
    'mathring': lambda s: s + '\N{COMBINING RING ABOVE}',
    'ddddot': lambda s: s + '\N{COMBINING DIAERESIS}\N{COMBINING DIAERESIS}',
    'dddot': lambda s: s + '\N{COMBINING DIAERESIS}\N{COMBINING DOT ABOVE}',
    'ddot': lambda s: s + '\N{COMBINING DIAERESIS}',
    'dot': lambda s: s + '\N{COMBINING DOT ABOVE}',
    'check': lambda s: s + '\N{COMBINING CARON}',
    'breve': lambda s: s + '\N{COMBINING BREVE}',
    'acute': lambda s: s + '\N{COMBINING ACUTE ACCENT}',
    'grave': lambda s: s + '\N{COMBINING GRAVE ACCENT}',
    'tilde': lambda s: s + '\N{COMBINING TILDE}',
    'hat': lambda s: s + '\N{COMBINING CIRCUMFLEX ACCENT}',
    'bar': lambda s: s + '\N{COMBINING OVERLINE}',
    'vec': lambda s: s + '\N{COMBINING RIGHT ARROW ABOVE}',
    'prime': lambda s: s + '\N{PRIME}',
    'prm': lambda s: s + '\N{PRIME}',
    # # Faces -- these are here for some compatibility with latex printing
    # 'bold': lambda s: s,
    # 'bm': lambda s: s,
    # 'cal': lambda s: s,
    # 'scr': lambda s: s,
    # 'frak': lambda s: s,
    # Brackets
    'norm': lambda s: '\N{DOUBLE VERTICAL LINE}' + s + '\N{DOUBLE VERTICAL LINE}',
    'avg': lambda s: '\N{MATHEMATICAL LEFT ANGLE BRACKET}' + s + '\N{MATHEMATICAL RIGHT ANGLE BRACKET}',
    'abs': lambda s: '\N{VERTICAL LINE}' + s + '\N{VERTICAL LINE}',
    'mag': lambda s: '\N{VERTICAL LINE}' + s + '\N{VERTICAL LINE}',
}

# VERTICAL OBJECTS


def HUP(symb):
    return U(f'{symb_2txt[symb]} UPPER HOOK')


def CUP(symb):
    return U(f'{symb_2txt[symb]} UPPER CORNER')


def MID(symb):
    return U(f'{symb_2txt[symb]} MIDDLE PIECE')


def EXT(symb):
    return U(f'{symb_2txt[symb]} EXTENSION')


def HLO(symb):
    return U(f'{symb_2txt[symb]} LOWER HOOK')


def CLO(symb):
    return U(f'{symb_2txt[symb]} LOWER CORNER')


# {} '('  ->  (extension, start, end, middle) 1-character
_xobj_unicode = {

    # vertical symbols
    #       (( ext, top, bot, mid ), c1)
    '(':    (( EXT('('), HUP('('), HLO('(') ), '('),
    ')':    (( EXT(')'), HUP(')'), HLO(')') ), ')'),
    '[':    (( EXT('['), CUP('['), CLO('[') ), '['),
    ']':    (( EXT(']'), CUP(']'), CLO(']') ), ']'),
    '{':    (( EXT('{}'), HUP('{'), HLO('{'), MID('{') ), '{'),
    '}':    (( EXT('{}'), HUP('}'), HLO('}'), MID('}') ), '}'),
    '|':    '\N{BOX DRAWINGS LIGHT VERTICAL}',

    '<':    (('\N{BOX DRAWINGS LIGHT VERTICAL}',
              '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}',
              '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}'), '<'),

    '>':    (('\N{BOX DRAWINGS LIGHT VERTICAL}',
              '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}',
              '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}'), '>'),

    'lfloor': (( EXT('['), EXT('['), CLO('[') ), '\N{LEFT FLOOR}'),
    'rfloor': (( EXT(']'), EXT(']'), CLO(']') ), '\N{RIGHT FLOOR}'),
    'lceil':  (( EXT('['), CUP('['), EXT('[') ), '\N{LEFT CEILING}'),
    'rceil':  (( EXT(']'), CUP(']'), EXT(']') ), '\N{RIGHT CEILING}'),

    'int':  (( EXT('int'), '\N{TOP HALF INTEGRAL}', '\N{BOTTOM HALF INTEGRAL}' ), '\N{INTEGRAL}'),
    'sum':  (( '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}', '_', '\N{OVERLINE}', '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}'), '\N{N-ARY SUMMATION}'),

    # horizontal objects
    # '-':   '-',
    '-':    '\N{BOX DRAWINGS LIGHT HORIZONTAL}',
    '_':    '\N{LOW LINE}',
    # We used to use this, but LOW LINE looks better for roots, as it's a
    # little lower (i.e., it lines up with the / perfectly.  But perhaps this
    # one would still be wanted for some cases?
    # '_':    '\N{HORIZONTAL SCAN LINE-9}',

    # diagonal objects '\' & '/' ?
    '/':    '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT}',
    '\\':   '\N{BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT}',
}

_xobj_ascii = {
    # vertical symbols
    #       (( ext, top, bot, mid ), c1)
    '(':    (( '|', '/', '\\' ), '('),
    ')':    (( '|', '\\', '/' ), ')'),

    # XXX this looks ugly
    #   '[':    (( '|', '-', '-' ), '['),
    #   ']':    (( '|', '-', '-' ), ']'),
    # XXX not so ugly :(
    '[':    (( '[', '[', '[' ), '['),
    ']':    (( ']', ']', ']' ), ']'),

    '{':    (( '|', '/', '\\', '<' ), '{'),
    '}':    (( '|', '\\', '/', '>' ), '}'),
    '|':    '|',

    '<':    (( '|', '/', '\\' ), '<'),
    '>':    (( '|', '\\', '/' ), '>'),

    'int':  ( ' | ', '  /', '/  ' ),

    # horizontal objects
    '-':    '-',
    '_':    '_',

    # diagonal objects '\' & '/' ?
    '/':    '/',
    '\\':   '\\',
}


def xobj(symb, length):
    """Construct spatial object of given length.

    return: [] of equal-length strings
    """
    if length <= 0:
        raise ValueError('Length should be greater than 0')

    # TODO robustify when no unicodedat available
    if _use_unicode:
        _xobj = _xobj_unicode
    else:
        _xobj = _xobj_ascii

    vinfo = _xobj[symb]

    c1 = top = bot = mid = None

    if not isinstance(vinfo, tuple):        # 1 entry
        ext = vinfo
    else:
        if isinstance(vinfo[0], tuple):     # (vlong), c1
            vlong = vinfo[0]
            c1 = vinfo[1]
        else:                               # (vlong), c1
            vlong = vinfo

        ext = vlong[0]

        try:
            top = vlong[1]
            bot = vlong[2]
            mid = vlong[3]
        except IndexError:
            pass

    if c1 is None:
        c1 = ext
    if top is None:
        top = ext
    if bot is None:
        bot = ext
    if mid is not None:
        if (length % 2) == 0:
            # even height, but we have to print it somehow anyway...
            # XXX is it ok?
            length += 1

    else:
        mid = ext

    if length == 1:
        return c1

    res = []
    next = (length - 2)//2
    nmid = (length - 2) - next*2

    res += [top]
    res += [ext]*next
    res += [mid]*nmid
    res += [ext]*next
    res += [bot]

    return res


def vobj(symb, height):
    """Construct vertical object of a given height.

    See Also
    ========

    xobj
    """
    return '\n'.join( xobj(symb, height) )


def hobj(symb, width):
    """Construct horizontal object of a given width.

    See Also
    ========

    xobj
    """
    return ''.join( xobj(symb, width) )


# RADICAL
# n -> symbol
root = {
    2: '\N{SQUARE ROOT}',   # '\N{RADICAL SYMBOL BOTTOM}'
    3: '\N{CUBE ROOT}',
    4: '\N{FOURTH ROOT}',
}


# RATIONAL
def VF(txt):
    return U(f'VULGAR FRACTION {txt}')


# (p,q) -> symbol
frac = {
    (1, 2): VF('ONE HALF'),
    (1, 3): VF('ONE THIRD'),
    (2, 3): VF('TWO THIRDS'),
    (1, 4): VF('ONE QUARTER'),
    (3, 4): VF('THREE QUARTERS'),
    (1, 5): VF('ONE FIFTH'),
    (2, 5): VF('TWO FIFTHS'),
    (3, 5): VF('THREE FIFTHS'),
    (4, 5): VF('FOUR FIFTHS'),
    (1, 6): VF('ONE SIXTH'),
    (5, 6): VF('FIVE SIXTHS'),
    (1, 8): VF('ONE EIGHTH'),
    (3, 8): VF('THREE EIGHTHS'),
    (5, 8): VF('FIVE EIGHTHS'),
    (7, 8): VF('SEVEN EIGHTHS'),
}


# atom symbols
_xsym = {
    '==':  ('=', '='),
    '<':   ('<', '<'),
    '>':   ('>', '>'),
    '<=':  ('<=', '\N{LESS-THAN OR EQUAL TO}'),
    '>=':  ('>=', '\N{GREATER-THAN OR EQUAL TO}'),
    '!=':  ('!=', '\N{NOT EQUAL TO}'),
    '*':   ('*', '\N{DOT OPERATOR}'),
    '-->': ('-->', '\N{EM DASH}' + '\N{EM DASH}' +
            '\N{BLACK RIGHT-POINTING TRIANGLE}'),
    '==>': ('==>', '\N{BOX DRAWINGS DOUBLE HORIZONTAL}' +
            '\N{BOX DRAWINGS DOUBLE HORIZONTAL}' +
            '\N{BLACK RIGHT-POINTING TRIANGLE}'),
    '.':   ('*', '\N{RING OPERATOR}'),
}


def xsym(sym):
    """Get symbology for a 'character'."""
    op = _xsym[sym]

    if _use_unicode:
        return op[1]
    else:
        return op[0]


# SYMBOLS

atoms_table = {
    # class                    how-to-display
    'Exp1':                    '\N{SCRIPT SMALL E}',
    'Pi':                      '\N{GREEK SMALL LETTER PI}',
    'Infinity':                '\N{INFINITY}',
    'NegativeInfinity':        '\N{INFINITY}' and ('-' + '\N{INFINITY}'),  # XXX what to do here
    # 'ImaginaryUnit':          '\N{GREEK SMALL LETTER IOTA}',
    # 'ImaginaryUnit':          '\N{MATHEMATICAL ITALIC SMALL I}',
    'ImaginaryUnit':           '\N{DOUBLE-STRUCK ITALIC SMALL I}',
    'EmptySet':                '\N{EMPTY SET}',
    'Naturals':                '\N{DOUBLE-STRUCK CAPITAL N}',
    'Naturals0':               ('\N{DOUBLE-STRUCK CAPITAL N}' and
                                ('\N{DOUBLE-STRUCK CAPITAL N}' +
                                 '\N{SUBSCRIPT ZERO}')),
    'Integers':                '\N{DOUBLE-STRUCK CAPITAL Z}',
    'Rationals':               '\N{DOUBLE-STRUCK CAPITAL Q}',
    'Reals':                   '\N{DOUBLE-STRUCK CAPITAL R}',
    'Union':                   '\N{UNION}',
    'SymmetricDifference':     '\N{INCREMENT}',
    'Intersection':            '\N{INTERSECTION}',
    'Ring':                    '\N{RING OPERATOR}'
}


def pretty_atom(atom_name, default=None):
    """Return pretty representation of an atom."""
    if _use_unicode:
        return atoms_table[atom_name]
    else:
        if default is not None:
            return default

        raise KeyError('only unicode')  # send it default printer


def pretty_symbol(symb_name):
    """Return pretty representation of a symbol."""
    # let's split symb_name into symbol + index
    # UC: beta1
    # UC: f_beta

    if not _use_unicode:
        return symb_name

    name, sups, subs = split_super_sub(symb_name)

    def translate(s):
        gG = greek_unicode.get(s)
        if gG is not None:
            return gG
        for key in sorted(modifier_dict, key=lambda k: len(k), reverse=True):
            if s.lower().endswith(key) and len(s) > len(key):
                return modifier_dict[key](translate(s[:-len(key)]))
        return s

    name = translate(name)

    # Let's prettify sups/subs. If it fails at one of them, pretty sups/subs are
    # not used at all.
    def pretty_list(l, mapping):
        result = []
        for s in l:
            pretty = mapping.get(s)
            if pretty is None:
                try:  # match by separate characters
                    pretty = ''.join([mapping[c] for c in s])
                except (TypeError, KeyError):
                    return
            result.append(pretty)
        return result

    pretty_sups = pretty_list(sups, sup)
    if pretty_sups is not None:
        pretty_subs = pretty_list(subs, sub)
    else:
        pretty_subs = None

    # glue the results into one string
    if pretty_subs is None:  # nice formatting of sups/subs did not work
        if subs:
            name += '_'+'_'.join([translate(s) for s in subs])
        if sups:
            name += '__'+'__'.join([translate(s) for s in sups])
        return name
    else:
        sups_result = ' '.join(pretty_sups)
        subs_result = ' '.join(pretty_subs)

    return ''.join([name, sups_result, subs_result])


def annotated(letter):
    """
    Return a stylised drawing of the letter ``letter``, together with
    information on how to put annotations (super- and subscripts to the
    left and to the right) on it.

    See pretty.py functions _print_meijerg, _print_hyper on how to use this
    information.
    """
    ucode_pics = {
        'F': (2, 0, 2, 0, '\N{BOX DRAWINGS LIGHT DOWN AND RIGHT}\N{BOX DRAWINGS LIGHT HORIZONTAL}\n'
                          '\N{BOX DRAWINGS LIGHT VERTICAL AND RIGHT}\N{BOX DRAWINGS LIGHT HORIZONTAL}\n'
                          '\N{BOX DRAWINGS LIGHT UP}'),
        'G': (3, 0, 3, 1,
              '\N{BOX DRAWINGS LIGHT ARC DOWN AND RIGHT}\N{BOX DRAWINGS LIGHT HORIZONTAL}\N{BOX DRAWINGS LIGHT ARC DOWN AND LEFT}\n'
              '\N{BOX DRAWINGS LIGHT VERTICAL}\N{BOX DRAWINGS LIGHT RIGHT}\N{BOX DRAWINGS LIGHT DOWN AND LEFT}\n'
              '\N{BOX DRAWINGS LIGHT ARC UP AND RIGHT}\N{BOX DRAWINGS LIGHT HORIZONTAL}\N{BOX DRAWINGS LIGHT ARC UP AND LEFT}')
    }
    ascii_pics = {
        'F': (3, 0, 3, 0, ' _\n|_\n|\n'),
        'G': (3, 0, 3, 1, ' __\n/__\n\\_|')
    }

    if _use_unicode:
        return ucode_pics[letter]
    else:
        return ascii_pics[letter]
