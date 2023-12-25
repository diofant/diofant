from collections import defaultdict

from ..core import Add, Function, Integer, Mul, Pow, Rational, pi
from ..functions import binomial, factorial, gamma, sin, sqrt
from ..polys import cancel, factor
from ..utilities import default_sort_key, ordered, sift
from ..utilities.iterables import uniq


def combsimp(expr):
    r"""
    Simplify combinatorial expressions.

    This function takes as input an expression containing factorials,
    binomials, Pochhammer symbol and other "combinatorial" functions,
    and tries to minimize the number of those functions and reduce
    the size of their arguments.

    The algorithm works by rewriting all combinatorial functions as
    expressions involving rising factorials (Pochhammer symbols) and
    applies recurrence relations and other transformations applicable
    to rising factorials, to reduce their arguments, possibly letting
    the resulting rising factorial to cancel. Rising factorials with
    the second argument being an integer are expanded into polynomial
    forms and finally all other rising factorial are rewritten in terms
    of more familiar functions. If the initial expression consisted of
    gamma functions alone, the result is expressed in terms of gamma
    functions. If the initial expression consists of gamma function
    with some other combinatorial, the result is expressed in terms of
    gamma functions.

    If the result is expressed using gamma functions, the following three
    additional steps are performed:

    1. Reduce the number of gammas by applying the reflection theorem
       gamma(x)*gamma(1-x) == pi/sin(pi*x).
    2. Reduce the number of gammas by applying the multiplication theorem
       gamma(x)*gamma(x+1/n)*...*gamma(x+(n-1)/n) == C*gamma(n*x).
    3. Reduce the number of prefactors by absorbing them into gammas, where
       possible.

    All transformation rules can be found (or was derived from) here:

    1. http://functions.wolfram.com/GammaBetaErf/Pochhammer/17/01/02/
    2. http://functions.wolfram.com/GammaBetaErf/Pochhammer/27/01/0005/

    Examples
    ========

    >>> combsimp(factorial(n)/factorial(n - 3))
    n*(n - 2)*(n - 1)
    >>> combsimp(binomial(n+1, k+1)/binomial(n, k))
    (n + 1)/(k + 1)

    """
    # as a rule of thumb, if the expression contained gammas initially, it
    # probably makes sense to retain them
    as_gamma = expr.has(gamma)
    as_factorial = expr.has(factorial)
    as_binomial = expr.has(binomial)

    expr = expr.replace(binomial,
                        lambda n, k: _rf((n - k + 1).expand(), k.expand())/_rf(1, k.expand()))
    expr = expr.replace(factorial,
                        lambda n: _rf(1, n.expand()))
    expr = expr.rewrite(gamma)
    expr = expr.replace(gamma,
                        lambda n: _rf(1, (n - 1).expand()))

    if as_gamma:
        expr = expr.replace(_rf,
                            lambda a, b: gamma(a + b)/gamma(a))
    else:
        expr = expr.replace(_rf,
                            lambda a, b: binomial(a + b - 1, b)*gamma(b + 1))

    def rule(n, k):
        coeff, rewrite = Integer(1), False

        cn, _n = n.as_coeff_Add()

        if _n and cn.is_Integer and cn:
            coeff *= _rf(_n + 1, cn)/_rf(_n - k + 1, cn)
            rewrite = True
            n = _n

        if rewrite:
            return coeff*binomial(n, k)

    expr = expr.replace(binomial, rule)

    def rule_gamma(expr, level=0):
        """Simplify products of gamma functions further."""
        if expr.is_Atom:
            return expr

        def gamma_rat(x):
            # helper to simplify ratios of gammas
            was = x.count(gamma)
            xx = x.replace(gamma, lambda n: _rf(1, (n - 1).expand()).replace(_rf, lambda a, b: gamma(a + b)/gamma(a)))
            if xx.count(gamma) < was:
                x = xx
            return x

        def gamma_factor(x):
            # return True if there is a gamma factor in shallow args
            if isinstance(x, gamma):
                return True
            if x.is_Add or x.is_Mul:
                return any(gamma_factor(xi) for xi in x.args)
            if x.is_Pow and (x.exp.is_integer or x.base.is_positive):
                return gamma_factor(x.base)
            return False

        # recursion step
        if level == 0:
            expr = expr.func(*[rule_gamma(x, level + 1) for x in expr.args])
            level += 1

        if not expr.is_Mul:
            return expr

        # non-commutative step
        if level == 1:
            args, nc = expr.args_cnc()
            if not args:
                return expr
            if nc:
                return rule_gamma(Mul._from_args(args), level + 1)*Mul._from_args(nc)
            level += 1

        # pure gamma handling, not factor absorbtion
        if level == 2:
            sifted = sift(expr.args, gamma_factor)
            gamma_ind = Mul(*sifted.pop(False, []))
            d = Mul(*sifted.pop(True, []))
            assert not sifted

            nd, dd = d.as_numer_denom()
            for ipass in range(2):
                args = list(ordered(Mul.make_args(nd)))
                for i, ni in enumerate(args):
                    if ni.is_Add:
                        ni, dd = Add(*[
                            rule_gamma(gamma_rat(a/dd), level + 1) for a in ni.args]
                        ).as_numer_denom()
                        args[i] = ni
                        if not dd.has(gamma):
                            break
                nd = Mul(*args)
                if ipass == 0 and not gamma_factor(nd):
                    break
                nd, dd = dd, nd  # now process in reversed order
            expr = gamma_ind*nd/dd
            if not (expr.is_Mul and (gamma_factor(dd) or gamma_factor(nd))):
                return expr
            level += 1

        # iteration until constant
        if level == 3:
            while True:
                was = expr
                expr = rule_gamma(expr, 4)
                if expr == was:
                    return expr

        numer_gammas = []
        denom_gammas = []
        numer_others = []
        denom_others = []

        def explicate(p):
            if p == 1:
                return None, []
            b, e = p.as_base_exp()
            if e.is_Integer:
                if isinstance(b, gamma):
                    return True, [b.args[0]]*e
                return False, [b]*e
            return False, [p]

        newargs = list(ordered(expr.args))
        while newargs:
            n, d = newargs.pop().as_numer_denom()
            isg, l = explicate(n)
            if isg:
                numer_gammas.extend(l)
            elif isg is False:
                numer_others.extend(l)
            isg, l = explicate(d)
            if isg:
                denom_gammas.extend(l)
            elif isg is False:
                denom_others.extend(l)

        # =========== level 2 work: pure gamma manipulation =========

        # Try to reduce the number of gamma factors by applying the
        # reflection formula gamma(x)*gamma(1-x) = pi/sin(pi*x)
        for gammas, numer, denom in [(
            numer_gammas, numer_others, denom_others),
                (denom_gammas, denom_others, numer_others)]:
            new = []
            while gammas:
                g1 = gammas.pop()
                if g1.is_integer:
                    new.append(g1)
                    continue
                for i, g2 in enumerate(gammas):
                    n = g1 + g2 - 1
                    if not n.is_Integer:
                        continue
                    numer.append(pi)
                    denom.append(sin(pi*g1))
                    gammas.pop(i)
                    if n > 0:
                        for k in range(n):
                            numer.append(1 - g1 + k)
                    elif n < 0:
                        for k in range(-n):
                            denom.append(-g1 - k)
                    break
                else:
                    new.append(g1)
            # /!\ updating IN PLACE
            gammas[:] = new

        # Try to reduce the number of gammas by using the duplication
        # theorem to cancel an upper and lower: gamma(2*s)/gamma(s) =
        # 2**(2*s + 1)/(4*sqrt(pi))*gamma(s + 1/2). Although this could
        # be done with higher argument ratios like gamma(3*x)/gamma(x),
        # this would not reduce the number of gammas as in this case.
        for ng, dg, no, do in [(numer_gammas, denom_gammas, numer_others,
                                denom_others),
                               (denom_gammas, numer_gammas, denom_others,
                                numer_others)]:

            while True:
                for x in ng:
                    for y in dg:
                        n = x - 2*y
                        if n.is_Integer:
                            break
                    else:
                        continue
                    break
                else:
                    break
                ng.remove(x)
                dg.remove(y)
                if n > 0:
                    for k in range(n):
                        no.append(2*y + k)
                elif n < 0:
                    for k in range(-n):
                        do.append(2*y - 1 - k)
                ng.append(y + Rational(1, 2))
                no.append(2**(2*y - 1))
                do.append(sqrt(pi))

        # Try to reduce the number of gamma factors by applying the
        # multiplication theorem (used when n gammas with args differing
        # by 1/n mod 1 are encountered).
        #
        # run of 2 with args differing by 1/2
        #
        # >>> combsimp(gamma(x)*gamma(x+Rational(1, 2)))
        # 2*sqrt(2)*2**(-2*x - 1/2)*sqrt(pi)*gamma(2*x)
        #
        # run of 3 args differing by 1/3 (mod 1)
        #
        # >>> combsimp(gamma(x)*gamma(x+Rational(1, 3))*gamma(x+Rational(2, 3)))
        # 6*3**(-3*x - 1/2)*pi*gamma(3*x)
        # >>> combsimp(gamma(x)*gamma(x+Rational(1, 3))*gamma(x+Rational(5, 3)))
        # 2*3**(-3*x - 1/2)*pi*(3*x + 2)*gamma(3*x)
        #
        def _run(coeffs):
            # find runs in coeffs such that the difference in terms (mod 1)
            # of t1, t2, ..., tn is 1/n
            u = list(uniq(coeffs))
            for i in range(len(u)):
                dj = ([((u[j] - u[i]) % 1, j) for j in range(i + 1, len(u))])
                for one, j in dj:
                    if one.numerator == 1 and one.denominator != 1:
                        n = Integer(one.denominator)
                        got = [i]
                        get = list(range(1, n))
                        for d, j in dj:
                            m = n*d
                            if m.is_Integer and m in get:
                                get.remove(m)
                                got.append(j)
                                if not get:
                                    break
                        else:
                            continue
                        for i, j in enumerate(got):
                            c = u[j]
                            coeffs.remove(c)
                            got[i] = c
                        return Integer(one.denominator), got[0], got[1:]

        def _mult_thm(gammas, numer, denom):
            # pull off and analyze the leading coefficient from each gamma arg
            # looking for runs in those Rationals

            # expr -> coeff + resid -> rats[resid] = coeff
            rats = defaultdict(list)
            for g in gammas:
                c, resid = g.as_coeff_Add()
                rats[resid].append(c)

            # look for runs in Rationals for each resid
            keys = sorted(rats, key=default_sort_key)
            for resid in keys:
                coeffs = sorted(rats[resid])
                new = []
                while True:
                    run = _run(coeffs)
                    if run is None:
                        break

                    # process the sequence that was found:
                    # 1) convert all the gamma functions to have the right
                    #    argument (could be off by an integer)
                    # 2) append the factors corresponding to the theorem
                    # 3) append the new gamma function

                    n, ui, other = run

                    # (1)
                    for u in other:
                        con = resid + u - 1
                        for k in range(int(u - ui)):
                            numer.append(con - k)

                    con = n*(resid + ui)  # for (2) and (3)

                    # (2)
                    numer.append((2*pi)**(Integer(n - 1)/2) *
                                 n**(Rational(1, 2) - con))
                    # (3)
                    new.append(con)

                # restore resid to coeffs
                rats[resid] = [resid + c for c in coeffs] + new

            # rebuild the gamma arguments
            g = []
            for resid in keys:
                g += rats[resid]
            # /!\ updating IN PLACE
            gammas[:] = g

        for l, numer, denom in [(numer_gammas, numer_others, denom_others),
                                (denom_gammas, denom_others, numer_others)]:
            _mult_thm(l, numer, denom)

        # =========== level >= 2 work: factor absorbtion =========

        if level >= 2:
            # Try to absorb factors into the gammas: x*gamma(x) -> gamma(x + 1)
            # and gamma(x)/(x - 1) -> gamma(x - 1)
            # This code (in particular repeated calls to find_fuzzy) can be very
            # slow.
            def find_fuzzy(l, x):
                if not l:
                    return
                S1, T1 = compute_ST(x)
                for y in l:
                    S2, T2 = inv[y]
                    if T1 != T2 or (not S1.intersection(S2) and
                                    (S1 != set() or S2 != set())):
                        continue
                    # XXX we want some simplification (e.g. cancel or
                    # simplify) but no matter what it's slow.
                    a = len(cancel(x/y).free_symbols)
                    b = len(x.free_symbols)
                    c = len(y.free_symbols)
                    # TODO is there a better heuristic?
                    if a == 0 and (b > 0 or c > 0):
                        return y

            # We thus try to avoid expensive calls by building the following
            # "invariants": For every factor or gamma function argument
            #   - the set of free symbols S
            #   - the set of functional components T
            # We will only try to absorb if T1==T2 and (S1 intersect S2 != emptyset
            # or S1 == S2 == emptyset)
            inv = {}

            def compute_ST(expr):
                if expr in inv:
                    return inv[expr]
                return (expr.free_symbols,
                        expr.atoms(Function) | {e.exp for e in expr.atoms(Pow)})

            def update_ST(expr):
                inv[expr] = compute_ST(expr)
            for expr in numer_gammas + denom_gammas + numer_others + denom_others:
                update_ST(expr)

            for gammas, numer, denom in [(
                numer_gammas, numer_others, denom_others),
                    (denom_gammas, denom_others, numer_others)]:
                new = []
                while gammas:
                    g = gammas.pop()
                    cont = True
                    while cont:
                        cont = False
                        y = find_fuzzy(numer, g)
                        if y is not None:
                            numer.remove(y)
                            if y != g:
                                numer.append(y/g)
                                update_ST(y/g)
                            g += 1
                            cont = True
                        y = find_fuzzy(denom, g - 1)
                        if y is not None:
                            denom.remove(y)
                            if y != g - 1:
                                numer.append((g - 1)/y)
                                update_ST((g - 1)/y)
                            g -= 1
                            cont = True
                    new.append(g)
                # /!\ updating IN PLACE
                gammas[:] = new

        # =========== rebuild expr ==================================

        return Mul(*[gamma(g) for g in numer_gammas]) \
            / Mul(*[gamma(g) for g in denom_gammas]) \
            * Mul(*numer_others) / Mul(*denom_others)

    # (for some reason we cannot use Basic.replace in this case)
    was = factor(expr)
    expr = rule_gamma(was)
    if expr != was:
        expr = factor(expr)

    if not as_gamma:
        if as_factorial:
            expr = expr.rewrite(factorial)
        elif as_binomial:
            expr = expr.rewrite(binomial)

    return expr


class _rf(Function):
    @classmethod
    def eval(cls, a, b):
        if b.is_Integer:
            if not b:
                return Integer(1)

            n, result = int(b), Integer(1)

            if n > 0:
                for i in range(n):
                    result *= a + i

                return result
            for i in range(1, -n + 1):
                result *= a - i

            return 1/result
        if b.is_Add:
            c, _b = b.as_coeff_Add()

            if c.is_Integer:
                if c > 0:
                    return _rf(a, _b)*_rf(a + _b, c)
                if c < 0:
                    return _rf(a, _b)/_rf(a + _b + c, -c)

        if a.is_Add:
            c, _a = a.as_coeff_Add()

            if c.is_Integer:
                if c > 0:
                    return _rf(_a, b)*_rf(_a + b, c)/_rf(_a, c)
                if c < 0:
                    return _rf(_a, b)*_rf(_a + c, -c)/_rf(_a + b + c, -c)
