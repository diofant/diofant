"""Useful utilities for higher level polynomial classes. """

import re
from collections import defaultdict

from ..core import Add, Mul, Pow, S, nan, oo, zoo
from ..core.exprtools import decompose_power
from .polyerrors import GeneratorsError, GeneratorsNeeded, PolynomialError
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
_re_gen = re.compile(r"^(.+?)(\d*)$")


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
        elif r.is_real is False and r.is_complex:
            return [[], [r]]
    if not all(r.is_number for r in roots):  # pragma: no cover
        raise NotImplementedError
    # see issue sympy/sympy#6137:
    # get the real part of the evaluated real and imaginary parts of each root
    key = [[i.n(2).as_real_imag()[0] for i in r.as_real_imag()] for r in roots]
    # make sure the parts were computed with precision
    if any(i._prec == 1 for k in key for i in k):  # pragma: no cover
        raise NotImplementedError("could not compute root with precision")
    # insert a key to indicate if the root has an imaginary part
    key = [(1 if i else 0, r, i) for r, i in key]
    key = sorted(zip(key, roots))
    # return the real and imaginary roots separately if desired
    if separated:
        r = []
        i = []
        for (im, _, _), v in key:
            if im:
                i.append(v)
            else:
                r.append(v)
        return r, i
    _, roots = zip(*key)
    return list(roots)


def _sort_gens(gens, **args):
    """Sort generators in a reasonably intelligent way. """
    opt = build_options(args)

    gens_order, wrt = {}, opt.wrt
    for i, gen in enumerate(opt.sort):
        gens_order[gen] = i + 1

    def order_key(gen):
        gen = str(gen)

        if wrt is not None:
            try:
                return -len(wrt) + wrt.index(gen), gen, 0
            except ValueError:
                pass

        name, index = _re_gen.match(gen).groups()

        if index:
            index = int(index)
        else:
            index = 0

        try:
            return gens_order[name], name, index
        except KeyError:
            pass

        try:
            return _gens_order[name], name, index
        except KeyError:
            pass

        return _max_order, name, index

    try:
        gens = sorted(gens, key=order_key)
    except TypeError:  # pragma: no cover
        pass

    return tuple(gens)


def _unify_gens(f_gens, g_gens):
    """Unify generators in a reasonably intelligent way. """
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


def _analyze_gens(gens):
    """Support for passing generators as `*gens` and `[gens]`. """
    if len(gens) == 1 and hasattr(gens[0], '__iter__'):
        return tuple(gens[0])
    else:
        return tuple(gens)


def _sort_factors(factors, **args):
    """Sort low-level factors in increasing 'complexity' order. """
    def order_if_multiple_key(factor):
        (f, n) = factor
        return len(f), n, f

    def order_no_multiple_key(f):
        return len(f), f

    if args.get('multiple', True):
        return sorted(factors, key=order_if_multiple_key)
    else:
        return sorted(factors, key=order_no_multiple_key)


def _not_a_coeff(expr):
    """Do not treat NaN and infinities as valid polynomial coefficients. """
    return expr in [nan, oo, -oo, zoo]


def _parallel_dict_from_expr_if_gens(exprs, opt):
    """Transform expressions into a multinomial form given generators. """
    k, indices = len(opt.gens), {}

    for i, g in enumerate(opt.gens):
        indices[g] = i

    polys = []

    for expr in exprs:
        poly = {}

        if expr.is_Equality:
            expr = expr.lhs - expr.rhs

        if not expr.is_commutative:
            raise PolynomialError('non-commutative expressions are not supported')

        for term in Add.make_args(expr):
            coeff, monom = [], [0]*k

            for factor in Mul.make_args(term):
                if not _not_a_coeff(factor) and factor.is_Number:
                    coeff.append(factor)
                else:
                    try:
                        base, exp = decompose_power(factor)

                        if exp < 0:
                            exp, base = -exp, Pow(base, -S.One)

                        monom[indices[base]] += exp
                    except KeyError:
                        if not factor.free_symbols.intersection(opt.gens):
                            coeff.append(factor)
                        else:
                            raise PolynomialError("%s contains an element of the generators set" % factor)

            monom = tuple(monom)

            if monom in poly:
                poly[monom] += Mul(*coeff)
            else:
                poly[monom] = Mul(*coeff)

        polys.append(poly)

    return polys, opt.gens


def _parallel_dict_from_expr_no_gens(exprs, opt):
    """Transform expressions into a multinomial form and figure out generators. """
    if opt.domain is not None:
        def _is_coeff(factor):
            return factor in opt.domain
    elif opt.extension is True:
        def _is_coeff(factor):
            return factor.is_algebraic
    elif opt.greedy is not False:
        def _is_coeff(factor):
            return False
    else:
        def _is_coeff(factor):
            return factor.is_number

    gens, reprs = set(), []

    for expr in exprs:
        terms = []

        if expr.is_Equality:
            expr = expr.lhs - expr.rhs

        if not expr.is_commutative:
            raise PolynomialError('non-commutative expressions are not supported')

        for term in Add.make_args(expr):
            coeff, elements = [], defaultdict(int)

            for factor in Mul.make_args(term):
                if not _not_a_coeff(factor) and (factor.is_Number or _is_coeff(factor)):
                    coeff.append(factor)
                else:
                    base, exp = decompose_power(factor)

                    if exp < 0:
                        exp, base = -exp, Pow(base, -S.One)

                    elements[base] += exp
                    gens.add(base)

            terms.append((coeff, elements))

        reprs.append(terms)

    if not gens:
        if len(exprs) == 1:
            arg = exprs[0]
        else:
            arg = (exprs,)

        raise GeneratorsNeeded("specify generators to give %s a meaning" % arg)

    gens = _sort_gens(gens, opt=opt)
    k, indices = len(gens), {}

    for i, g in enumerate(gens):
        indices[g] = i

    polys = []

    for terms in reprs:
        poly = {}

        for coeff, term in terms:
            monom = [0]*k

            for base, exp in term.items():
                monom[indices[base]] = exp

            monom = tuple(monom)

            if monom in poly:
                poly[monom] += Mul(*coeff)
            else:
                poly[monom] = Mul(*coeff)

        polys.append(poly)

    return polys, tuple(gens)


def _dict_from_expr_if_gens(expr, opt):
    """Transform an expression into a multinomial form given generators. """
    (poly,), gens = _parallel_dict_from_expr_if_gens((expr,), opt)
    return poly, gens


def _dict_from_expr_no_gens(expr, opt):
    """Transform an expression into a multinomial form and figure out generators. """
    (poly,), gens = _parallel_dict_from_expr_no_gens((expr,), opt)
    return poly, gens


def parallel_dict_from_expr(exprs, **args):
    """Transform expressions into a multinomial form. """
    reps, opt = _parallel_dict_from_expr(exprs, build_options(args))
    return reps, opt.gens


def _parallel_dict_from_expr(exprs, opt):
    """Transform expressions into a multinomial form. """
    if opt.expand is not False:
        exprs = [ expr.expand() for expr in exprs ]

    if opt.gens:
        reps, gens = _parallel_dict_from_expr_if_gens(exprs, opt)
    else:
        reps, gens = _parallel_dict_from_expr_no_gens(exprs, opt)

    return reps, opt.clone({'gens': gens})


def dict_from_expr(expr, **args):
    """Transform an expression into a multinomial form. """
    rep, opt = _dict_from_expr(expr, build_options(args))
    return rep, opt.gens


def _dict_from_expr(expr, opt):
    """Transform an expression into a multinomial form. """
    if opt.expand is not False:
        expr = expr.expand()

    if opt.gens:
        rep, gens = _dict_from_expr_if_gens(expr, opt)
    else:
        rep, gens = _dict_from_expr_no_gens(expr, opt)

    return rep, opt.clone({'gens': gens})


def expr_from_dict(rep, *gens):
    """Convert a multinomial form into an expression. """
    result = []

    for monom, coeff in rep.items():
        term = [coeff]
        for g, m in zip(gens, monom):
            if m:
                term.append(Pow(g, m))

        result.append(Mul(*term))

    return Add(*result)


parallel_dict_from_basic = parallel_dict_from_expr
dict_from_basic = dict_from_expr
basic_from_dict = expr_from_dict


def _dict_reorder(rep, gens, new_gens):
    """Reorder levels using dict representation. """
    gens = list(gens)

    monoms = rep.keys()
    coeffs = rep.values()

    new_monoms = [ [] for _ in range(len(rep)) ]
    used_indices = set()

    for gen in new_gens:
        try:
            j = gens.index(gen)
            used_indices.add(j)

            for M, new_M in zip(monoms, new_monoms):
                new_M.append(M[j])
        except ValueError:
            for new_M in new_monoms:
                new_M.append(0)

    for i, _ in enumerate(gens):
        if i not in used_indices:
            for monom in monoms:
                if monom[i]:
                    raise GeneratorsError("unable to drop generators")

    return map(tuple, new_monoms), coeffs
