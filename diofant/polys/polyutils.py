"""Useful utilities for higher level polynomial classes."""

import re

from ..core import Add, Mul, Pow
from ..core.exprtools import decompose_power
from ..utilities import default_sort_key
from .polyerrors import GeneratorsNeededError, PolynomialError
from .polyoptions import build_options


_gens_order = {
    'a': 301, 'b': 302, 'c': 303, 'd': 304,
    'e': 305, 'f': 306, 'g': 307, 'h': 308,
    'i': 309, 'j': 310, 'k': 311, 'l': 312,
    'm': 313, 'n': 314, 'o': 315, 'p': 216,
    'q': 217, 'r': 218, 's': 219, 't': 220,
    'u': 221, 'v': 222, 'w': 223, 'x': 124,
    'y': 125, 'z': 126,
}

_max_order = 1000
_re_gen = re.compile(r'^(.+?)(\d*)$')


def _nsort(roots, separated=False):
    """Sort the numerical roots putting the real roots first, then sorting
    according to real and imaginary parts. If ``separated`` is True, then
    the real and imaginary roots will be returned in two lists, respectively.

    This routine tries to avoid issue sympy/sympy#6137 by separating the roots into real
    and imaginary parts before evaluation. In addition, the sorting will raise
    an error if any computation cannot be done with precision.
    """
    if len(roots) == 1:
        if not separated:
            return list(roots)
        r = list(roots)[0]
        if r.is_extended_real:
            return [[r], []]
        if r.is_real is False and r.is_complex:
            return [[], [r]]
        raise NotImplementedError
    all_numbers = all(r.is_number for r in roots)
    if not all_numbers:
        raise NotImplementedError
    # see issue sympy/sympy#6137:
    # get the real part of the evaluated real and imaginary parts of each root
    key = [[i.evalf(2).as_real_imag()[0] for i in r.as_real_imag()] for r in roots]
    # insert a key to indicate if the root has an imaginary part
    key = [(1 if i else 0, r, -abs(i), i.is_positive) for r, i in key]
    key = sorted(zip(key, roots))
    # return the real and imaginary roots separately if desired
    if separated:
        r = []
        i = []
        for (im, _, _, _), v in key:
            if im:
                i.append(v)
            else:
                r.append(v)
        return r, i
    _, roots = zip(*key)
    return list(roots)


def _sort_gens(gens, **args):
    """Sort generators in a reasonably intelligent way."""
    opt = build_options(args)
    gens_order, wrt = _gens_order.copy(), opt.wrt

    for i, gen in enumerate(opt.sort, start=1):
        gens_order[gen] = i

    def order_key(gen):
        gen = str(gen)

        if wrt is not None:
            try:
                return -len(wrt) + wrt.index(gen), gen, 0
            except ValueError:
                pass

        name, index = _re_gen.match(gen).groups()
        index = int(index) if index else 0

        return gens_order.get(name, _max_order), name, index

    return tuple(sorted(gens, key=order_key))


def _unify_gens(f_gens, g_gens):
    """Unify generators in a reasonably intelligent way."""
    f_gens = list(f_gens)
    g_gens = list(g_gens)

    if f_gens == g_gens:
        return tuple(f_gens)

    gens, common, k = [], [], 0

    for gen in f_gens:
        if gen in g_gens:
            common.append(gen)

    for i, gen in enumerate(g_gens):
        if gen in common:
            g_gens[i], k = common[k], k + 1

    for gen in common:
        i = f_gens.index(gen)

        gens.extend(f_gens[:i])
        f_gens = f_gens[i + 1:]

        i = g_gens.index(gen)

        gens.extend(g_gens[:i])
        g_gens = g_gens[i + 1:]

        gens.append(gen)

    gens.extend(f_gens)
    gens.extend(g_gens)

    return tuple(gens)


def _find_gens(exprs, opt):
    """Find generators in a reasonably intelligent way."""
    if opt.domain is not None:
        def _is_coeff(factor):
            return factor in opt.domain
    elif opt.extension is True:
        def _is_coeff(factor):
            return factor.is_number and factor.is_algebraic
    elif opt.greedy is not False:
        def _is_coeff(factor):
            return factor.is_Number and factor.is_finite is not False
    else:
        def _is_coeff(factor):
            return factor.is_number and factor.is_finite is not False

    gens = set()

    for expr in exprs:
        for term in Add.make_args(expr):
            for factor in Mul.make_args(term):
                try:
                    if factor.is_Add and opt.expand:
                        gens |= set(_find_gens([factor], opt))
                    elif not _is_coeff(factor):
                        base, exp = decompose_power(factor)
                        if exp < 0:
                            base = Pow(base, -1)

                        if opt.expand and exp > 1:
                            gens |= set(_find_gens([base], opt))
                        else:
                            gens.add(base)
                except GeneratorsNeededError:
                    pass

    if not gens:
        raise GeneratorsNeededError(f'specify generators to give {exprs} a meaning')

    return _sort_gens(gens, opt=opt)


def _sort_factors(factors, **args):
    """Sort low-level factors in increasing 'complexity' order."""
    def order_if_multiple_key(factor):
        f, n = factor
        return len(f), n, default_sort_key(f)

    def order_no_multiple_key(f):
        return len(f), default_sort_key(f)

    if args.get('multiple', True):
        return sorted(factors, key=order_if_multiple_key)
    return sorted(factors, key=order_no_multiple_key)


def _parallel_dict_from_expr_if_gens(exprs, opt):
    """Transform expressions into a multinomial form given generators."""
    indices = {g: i for i, g in enumerate(opt.gens)}
    zero_monom = [0]*len(opt.gens)
    polys = []

    for expr in exprs:
        poly = {}

        for term in Add.make_args(expr):
            coeff, monom = [], zero_monom.copy()

            for factor in Mul.make_args(term):
                base, exp = decompose_power(factor)
                if exp < 0:
                    exp, base = -exp, Pow(base, -1)
                try:
                    monom[indices[base]] += exp
                    continue
                except KeyError as exc:
                    if factor.free_symbols & set(opt.gens):
                        raise PolynomialError(f'{factor} contains an element'
                                              ' of the generators set') from exc

                coeff.append(factor)

            monom = tuple(monom)
            poly[monom] = Mul(*coeff) + poly.get(monom, 0)

        polys.append(poly)

    return polys


def parallel_dict_from_expr(exprs, **args):
    """Transform expressions into a multinomial form."""
    reps, opt = _parallel_dict_from_expr(exprs, build_options(args))
    return reps, opt.gens


def _parallel_dict_from_expr(exprs, opt):
    """Transform expressions into a multinomial form."""
    if any(not expr.is_commutative for expr in exprs):
        raise PolynomialError('non-commutative expressions are not supported')

    if opt.expand is not False:
        exprs = [expr.expand() for expr in exprs]

    if not opt.gens:
        opt = opt.clone({'gens': _find_gens(exprs, opt)})

    reps = _parallel_dict_from_expr_if_gens(exprs, opt)

    return reps, opt.clone()
