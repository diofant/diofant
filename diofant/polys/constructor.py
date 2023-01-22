"""Tools for constructing domains for expressions."""

from ..core import I
from ..core.sympify import sympify
from ..domains import EX, QQ, RR, ZZ
from ..domains.realfield import RealField
from ..utilities import ordered
from .polyerrors import GeneratorsNeededError
from .polyoptions import build_options
from .polyutils import parallel_dict_from_expr


def _construct_simple(coeffs, opt):
    """Handle simple domains, e.g.: ZZ, QQ, RR and algebraic domains."""
    result, rationals, reals, algebraics = {}, False, False, False

    if opt.extension is False:
        def is_algebraic(coeff):
            return False
    else:
        def is_algebraic(coeff):
            return coeff.is_number and coeff.is_algebraic and not coeff.is_rational

    for coeff in coeffs:
        if coeff.is_Rational:
            if not coeff.is_Integer:
                rationals = True
        elif coeff.is_Float:
            if not algebraics:
                reals = True
            else:
                # there are both reals and algebraics -> EX
                return False
        elif is_algebraic(coeff):
            if not reals:
                algebraics = True
            else:
                # there are both algebraics and reals -> EX
                return False
        else:
            # this is a composite domain, e.g. ZZ[X], EX
            return

    if algebraics:
        domain, result = _construct_algebraic(coeffs, opt)
    else:
        if reals:
            # Use the maximum precision of all coefficients for the RR's
            # precision
            max_prec = max(c._prec for c in coeffs)
            domain = RealField(prec=max_prec)
        else:
            if opt.field or rationals:
                domain = QQ
            else:
                domain = ZZ

        result = []

        for coeff in coeffs:
            result.append(domain.convert(coeff))

    return domain, result


def _construct_algebraic(coeffs, opt):
    """We know that coefficients are algebraic so construct the extension."""
    result, exts = [], set()

    for coeff in coeffs:
        if coeff.is_Rational:
            coeff = (None, 0, QQ.convert(coeff))
        else:
            a = coeff.as_coeff_add()[0]
            coeff -= a

            b = coeff.as_coeff_mul()[0]
            coeff /= b

            exts.add(coeff)

            a = QQ.convert(a)
            b = QQ.convert(b)

            coeff = (coeff, b, a)

        result.append(coeff)

    exts = list(ordered(exts))

    if all(e.is_real for e in exts):
        domain = QQ.algebraic_field(*exts)
    else:
        ground_exts = list(ordered(set().union(*[_.as_real_imag() for _ in exts])))
        domain = QQ.algebraic_field(*ground_exts).algebraic_field(I)

    H = [domain.from_expr(e).rep for e in exts]

    for i, (coeff, a, b) in enumerate(result):
        if coeff is not None:
            coeff = a*domain.dtype(H[exts.index(coeff)]) + b
        else:
            coeff = domain.dtype([b])

        result[i] = coeff

    return domain, result


def _construct_composite(coeffs, opt):
    """Handle composite domains, e.g.: ZZ[X], QQ[X], ZZ(X), QQ(X)."""
    numers, denoms = [], []

    for coeff in coeffs:
        numer, denom = coeff.as_numer_denom()

        numers.append(numer)
        denoms.append(denom)

    try:
        polys, gens = parallel_dict_from_expr(numers + denoms)  # XXX: sorting
    except GeneratorsNeededError:
        return

    if opt.composite is None:
        if any(g.is_number and not g.is_transcendental for g in gens):
            return  # generators are number-like so lets better use EX

        all_symbols = set()

        for gen in gens:
            symbols = gen.free_symbols

            if all_symbols & symbols:
                return  # there could be algebraic relations between generators
            all_symbols |= symbols

    n = len(gens)
    k = len(polys)//2

    numers = polys[:k]
    denoms = polys[k:]

    if opt.field:
        fractions = True
    else:
        fractions, zeros = False, (0,)*n

        for denom in denoms:
            if len(denom) > 1 or zeros not in denom:
                fractions = True
                break

    coeffs = set()

    if not fractions:
        for numer, denom in zip(numers, denoms):
            denom = denom[zeros]

            for monom, coeff in numer.items():
                coeff /= denom
                coeffs.add(coeff)
                numer[monom] = coeff
    else:
        for numer, denom in zip(numers, denoms):
            coeffs.update(list(numer.values()))
            coeffs.update(list(denom.values()))

    rationals, reals = False, False

    for coeff in coeffs:
        if coeff.is_Rational:
            if not coeff.is_Integer:
                rationals = True
        elif coeff.is_Float:
            reals = True
            break
        else:
            raise NotImplementedError

    if reals:
        ground = RR
    elif rationals:
        ground = QQ
    else:
        ground = ZZ

    result = []

    if not fractions:
        domain = ground.inject(*gens)

        for numer in numers:
            for monom, coeff in numer.items():
                numer[monom] = ground.convert(coeff)

            result.append(domain(numer))
    else:
        domain = ground.inject(*gens).field

        for numer, denom in zip(numers, denoms):
            for monom, coeff in numer.items():
                numer[monom] = ground.convert(coeff)

            for monom, coeff in denom.items():
                denom[monom] = ground.convert(coeff)

            result.append(domain((numer, denom)))

    return domain, result


def _construct_expression(coeffs, opt):
    """The last resort case, i.e. use the expression domain."""
    domain, result = EX, []

    for coeff in coeffs:
        result.append(domain.convert(coeff))

    return domain, result


def construct_domain(obj, **args):
    """Construct a minimal domain for the list of coefficients."""
    opt = build_options(args)

    if hasattr(obj, '__iter__'):
        if isinstance(obj, dict):
            if not obj:
                monoms, coeffs = [], []
            else:
                monoms, coeffs = zip(*obj.items())
        else:
            coeffs = obj
    else:
        coeffs = [obj]

    coeffs = list(map(sympify, coeffs))
    result = _construct_simple(coeffs, opt)

    if result is not None:
        if result is not False:
            domain, coeffs = result
        else:
            domain, coeffs = _construct_expression(coeffs, opt)
    else:
        if opt.composite is False:
            result = None
        else:
            result = _construct_composite(coeffs, opt)

        if result is not None:
            domain, coeffs = result
        else:
            domain, coeffs = _construct_expression(coeffs, opt)

    if hasattr(obj, '__iter__'):
        if isinstance(obj, dict):
            return domain, dict(zip(monoms, coeffs))
        return domain, coeffs
    return domain, coeffs[0]
