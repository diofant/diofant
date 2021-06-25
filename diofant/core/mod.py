from .function import Function
from .numbers import nan


class Mod(Function):
    """Represents a modulo operation on symbolic expressions.

    Receives two arguments, dividend p and divisor q.

    The convention used is the same as Python's: the remainder always has the
    same sign as the divisor.

    Examples
    ========

    >>> x**2 % y
    x**2%y
    >>> _.subs({x: 5, y: 6})
    1

    """

    @classmethod
    def eval(cls, p, q):
        from ..polys.polytools import gcd
        from .add import Add
        from .exprtools import gcd_terms
        from .mul import Mul
        from .numbers import Integer

        def doit(p, q):
            """Try to return p % q if both are numbers or +/-p is known
            to be less than or equal q.

            """
            if p.is_infinite or q.is_infinite:
                return nan
            if (p == q or p == -q or
                    p.is_Pow and p.exp.is_Integer and p.base == q or
                    p.is_integer and q == 1):
                return Integer(0)

            if q.is_Number:
                if p.is_Number:
                    return p % q
                if q == 2:
                    if p.is_even:
                        return Integer(0)
                    elif p.is_odd:
                        return Integer(1)

            # by ratio
            r = p/q
            try:
                d = int(r)
            except TypeError:
                pass
            else:
                rv = p - d*q
                if (rv*q).is_nonnegative:
                    return rv
                elif (rv*q).is_nonpositive:
                    return rv + q

            # by difference
            d = p - q
            if d.is_negative:
                if q.is_negative:
                    return d
                elif q.is_positive:
                    return p

        rv = doit(p, q)
        if rv is not None:
            return rv

        # denest
        if isinstance(p, cls):
            # easy
            qinner = p.args[1]
            if qinner == q:
                return p
            # XXX other possibilities?

        # extract gcd; any further simplification should be done by the user
        G = gcd(p, q)
        if G != 1:
            p, q = [
                gcd_terms(i/G, clear=False, fraction=False) for i in (p, q)]
        pwas, qwas = p, q

        # simplify terms
        # (x + y + 2) % x -> Mod(y + 2, x)
        if p.is_Add:
            args = []
            for i in p.args:
                a = cls(i, q)
                if a.count(cls) > i.count(cls):
                    args.append(i)
                else:
                    args.append(a)
            if args != list(p.args):
                p = Add(*args)

        else:
            # handle coefficients if they are not Rational
            # since those are not handled by factor_terms
            # e.g. Mod(.6*x, .3*y) -> 0.3*Mod(2*x, y)
            cp, p = p.as_coeff_Mul()
            cq, q = q.as_coeff_Mul()
            ok = False
            if not cp.is_Rational or not cq.is_Rational:
                r = cp % cq
                if r == 0:
                    G *= cq
                    p *= int(cp/cq)
                    ok = True
            if not ok:
                p = cp*p
                q = cq*q

        # simple -1 extraction
        if p.could_extract_minus_sign() and q.could_extract_minus_sign():
            G, p, q = [-i for i in (G, p, q)]

        # check again to see if p and q can now be handled as numbers
        rv = doit(p, q)
        if rv is not None:
            return rv*G

        # put 1.0 from G on inside
        if G.is_Float and G == 1:
            p *= G
            return cls(p, q, evaluate=False)
        elif G.is_Mul and G.args[0].is_Float and G.args[0] == 1:
            p = G.args[0]*p
            G = Mul._from_args(G.args[1:])
        return G*cls(p, q, evaluate=(p, q) != (pwas, qwas))

    def _eval_is_integer(self):
        p, q = self.args
        if p.is_integer and q.is_integer and q.is_nonzero:
            return True

    def _eval_is_nonnegative(self):
        p, q = self.args
        if p.is_real and q.is_real and q.is_positive:
            return True

    def _eval_is_nonpositive(self):
        p, q = self.args
        if p.is_real and q.is_real and q.is_negative:
            return True
