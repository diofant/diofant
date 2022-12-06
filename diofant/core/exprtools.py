"""Tools for manipulating of large commutative expressions."""

import numbers
from collections import defaultdict

from ..utilities import default_sort_key, ordered, variations
from ..utilities.iterables import common_prefix, common_suffix
from .add import Add
from .basic import Basic, preorder_traversal
from .compatibility import is_sequence, iterable
from .containers import Dict, Tuple
from .coreerrors import NonCommutativeExpression
from .expr import Expr
from .mul import Mul, _keep_coeff
from .numbers import I, Integer, Number, Rational, oo
from .power import Pow
from .symbol import Dummy
from .sympify import sympify


def _isnumber(i):
    return isinstance(i, (numbers.Integral, float)) or i.is_Number


def decompose_power(expr):
    """
    Decompose power into symbolic base and integer exponent.

    This is strictly only valid if the exponent from which
    the integer is extracted is itself an integer or the
    base is positive. These conditions are assumed and not
    checked here.

    Examples
    ========

    >>> decompose_power(x)
    (x, 1)
    >>> decompose_power(x**2)
    (x, 2)
    >>> decompose_power(x**(2*y))
    (x**y, 2)
    >>> decompose_power(x**(2*y/3))
    (x**(y/3), 2)

    """
    base, exp = expr.as_base_exp()

    if exp.is_Number:
        if exp.is_Rational:
            if not exp.is_Integer:
                base = Pow(base, Rational(1, exp.denominator))

            exp = exp.numerator
        else:
            base, exp = expr, 1
    else:
        exp, tail = exp.as_coeff_Mul(rational=True)

        if exp == -1:
            base, exp = Pow(base, tail), -1
        elif exp != 1:
            tail = _keep_coeff(Rational(1, exp.denominator), tail)
            base, exp = Pow(base, tail), exp.numerator
        else:
            base, exp = expr, 1

    return base, exp


class Factors:
    """Efficient representation of ``f_1*f_2*...*f_n``."""

    def __init__(self, factors=None):
        """Initialize Factors from dict or expr.

        Examples
        ========

        >>> e = 2*x**3
        >>> Factors(e)
        Factors({2: 1, x: 3})
        >>> Factors(e.as_powers_dict())
        Factors({2: 1, x: 3})
        >>> f = _
        >>> f.factors  # underlying dictionary
        {2: 1, x: 3}
        >>> f.gens  # base of each factor
        frozenset({2, x})
        >>> Factors(0)
        Factors({0: 1})
        >>> Factors(I)
        Factors({I: 1})

        Notes
        =====

        Although a dictionary can be passed, only minimal checking is
        performed: powers of -1 and I are made canonical.

        """
        if isinstance(factors, (numbers.Integral, float)):
            factors = sympify(factors)

        if isinstance(factors, Factors):
            factors = factors.factors.copy()
        elif factors in (None, 1):
            factors = {}
        elif factors == 0:
            factors = {Integer(0): Integer(1)}
        elif isinstance(factors, Number):
            n = factors
            factors = {}
            if n < 0:
                factors[Integer(-1)] = Integer(1)
                n = -n
            if n.is_Float or n.is_Integer or n is oo:
                factors[n] = Integer(1)
            elif n.is_Rational and n != 1:
                # since we're processing Numbers, the denominator is
                # stored with a negative exponent; all other factors
                # are left .
                if n.numerator != 1:
                    factors[Integer(n.numerator)] = Integer(1)
                factors[Integer(n.denominator)] = Integer(-1)
            else:  # pragma: no cover
                raise ValueError(f'Expected Float|Rational|Integer, not {n}')
        elif isinstance(factors, Basic) and not factors.args:
            factors = {factors: Integer(1)}
        elif isinstance(factors, Expr):
            c, nc = factors.args_cnc()
            i = c.count(I)
            for _ in range(i):
                c.remove(I)
            factors = dict(Mul._from_args(c).as_powers_dict())
            if i:
                factors[I] = Integer(1)*i
            if nc:
                factors[Mul(*nc, evaluate=False)] = Integer(1)
        else:
            factors = factors.copy()  # /!\ should be dict-like

            # tidy up -/+1 and I exponents if Rational

            handle = []
            for k in factors:
                if k in (I, -1, 1):
                    handle.append(k)
            if handle:
                i1 = Integer(1)
                for k in handle:
                    if not _isnumber(factors[k]):
                        continue
                    i1 *= k**factors.pop(k)
                if i1 != 1 or i1.is_Float:
                    for a in i1.args if i1.is_Mul else [i1]:  # at worst, -1.0*I*(-1)**e
                        if a == -1 and not a.is_Float:
                            factors[a] = Integer(1)
                        elif a is I:
                            factors[I] = Integer(1)
                        elif a.is_Pow:
                            if -1 not in factors:
                                factors[Integer(-1)] = Integer(0)
                            factors[Integer(-1)] += a.exp
                        elif a == 1:
                            factors[a] = Integer(1)
                        elif a == -1:
                            factors[-a] = Integer(1)
                            factors[Integer(-1)] = Integer(1)
                        else:  # pragma: no cover
                            raise RuntimeError(f'unexpected factor in i1: {a}')

        self.factors = factors
        self.gens = frozenset(factors)

    def __hash__(self):
        keys = tuple(ordered(self.factors))
        values = tuple(self.factors[k] for k in keys)
        return hash((keys, values))

    def __repr__(self):
        return f"Factors({{{', '.join([f'{k!s}: {v!s}' for k, v in ordered(self.factors.items())])}}})"

    @property
    def is_zero(self):
        """
        >>> Factors(0).is_zero
        True

        """
        f = self.factors
        return len(f) == 1 and 0 in f

    @property
    def is_one(self):
        """
        >>> Factors(1).is_one
        True

        """
        return not self.factors

    def as_expr(self):
        """Return the underlying expression.

        Examples
        ========

        >>> Factors((x*y**2).as_powers_dict()).as_expr()
        x*y**2

        """
        args = []
        for factor, exp in self.factors.items():
            if exp != 1:
                b, e = factor.as_base_exp()
                if isinstance(exp, int):
                    e = _keep_coeff(Integer(exp), e)
                elif isinstance(exp, Rational):
                    e = _keep_coeff(exp, e)
                else:
                    e *= exp
                args.append(b**e)
            else:
                args.append(factor)
        return Mul(*args)

    def mul(self, other):
        """Return Factors of ``self * other``.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.mul(b)
        Factors({x: 2, y: 3, z: -1})
        >>> a*b
        Factors({x: 2, y: 3, z: -1})

        """
        if not isinstance(other, Factors):
            other = Factors(other)
        if any(f.is_zero for f in (self, other)):
            return Factors(Integer(0))
        factors = dict(self.factors)

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = factors[factor] + exp

                if not exp:
                    del factors[factor]
                    continue

            factors[factor] = exp

        return Factors(factors)

    def normal(self, other):
        """Return ``self`` and ``other`` with ``gcd`` removed from each.
        The only differences between this and method ``div`` is that this
        is 1) optimized for the case when there are few factors in common and
        2) this does not raise an error if ``other`` is zero.

        See Also
        ========
        div

        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return Factors(), Factors(Integer(0))
            if self.is_zero:
                return Factors(Integer(0)), Factors()

        self_factors = dict(self.factors)
        other_factors = dict(other.factors)

        for factor, self_exp in self.factors.items():
            try:
                other_exp = other.factors[factor]
            except KeyError:
                continue

            exp = self_exp - other_exp

            if not exp:
                del self_factors[factor]
                del other_factors[factor]
            elif _isnumber(exp):
                if exp > 0:
                    self_factors[factor] = exp
                    del other_factors[factor]
                else:
                    del self_factors[factor]
                    other_factors[factor] = -exp
            else:
                r = self_exp.extract_additively(other_exp)
                if r is not None:
                    assert r
                    self_factors[factor] = r
                    del other_factors[factor]
                else:
                    sc, sa = self_exp.as_coeff_Add()
                    if sc:
                        oc, oa = other_exp.as_coeff_Add()
                        diff = sc - oc
                        if diff > 0:
                            self_factors[factor] -= oc
                            other_exp = oa
                        elif diff < 0:
                            self_factors[factor] -= sc
                            other_factors[factor] -= sc
                            other_exp = oa - diff
                        else:
                            self_factors[factor] = sa
                            other_exp = oa
                    if other_exp:
                        other_factors[factor] = other_exp
                    else:
                        del other_factors[factor]

        return Factors(self_factors), Factors(other_factors)

    def div(self, other):
        """Return ``self`` and ``other`` with ``gcd`` removed from each.
        This is optimized for the case when there are many factors in common.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> a.div(a)
        (Factors({}), Factors({}))
        >>> a.div(x*z)
        (Factors({y: 2}), Factors({z: 1}))

        The ``/`` operator only gives ``quo``:

        >>> a/x
        Factors({y: 2})

        Factors treats its factors as though they are all in the numerator, so
        if you violate this assumption the results will be correct but will
        not strictly correspond to the numerator and denominator of the ratio:

        >>> a.div(x/z)
        (Factors({y: 2}), Factors({z: -1}))

        Factors is also naive about bases: it does not attempt any denesting
        of Rational-base terms, for example the following does not become
        2**(2*x)/2.

        >>> Factors(2**(2*x + 2)).div(Integer(8))
        (Factors({2: 2*x + 2}), Factors({8: 1}))

        factor_terms can clean up such Rational-bases powers:

        >>> n, d = Factors(2**(2*x + 2)).div(Integer(8))
        >>> n.as_expr()/d.as_expr()
        2**(2*x + 2)/8
        >>> factor_terms(_)
        2**(2*x)/2

        """
        quo, rem = dict(self.factors), {}

        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                raise ZeroDivisionError
            if self.is_zero:
                return Factors(Integer(0)), Factors()

        for factor, exp in other.factors.items():
            if factor in quo:
                d = quo[factor] - exp
                if _isnumber(d):
                    if d <= 0:
                        del quo[factor]

                    if d >= 0:
                        if d:
                            quo[factor] = d
                    else:
                        exp = -d
                        rem[factor] = exp
                else:
                    r = quo[factor].extract_additively(exp)
                    if r is not None:
                        assert r
                        quo[factor] = r
                    else:
                        other_exp = exp
                        sc, sa = quo[factor].as_coeff_Add()
                        if sc:
                            oc, oa = other_exp.as_coeff_Add()
                            diff = sc - oc
                            if diff > 0:
                                quo[factor] -= oc
                                other_exp = oa
                            elif diff < 0:
                                quo[factor] -= sc
                                other_exp = oa - diff
                            else:
                                quo[factor] = sa
                                other_exp = oa
                        if other_exp:
                            rem[factor] = other_exp
                        else:
                            assert factor not in rem
            else:
                rem[factor] = exp

        return Factors(quo), Factors(rem)

    def quo(self, other):
        """Return numerator Factor of ``self / other``.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.quo(b)  # same as a/b
        Factors({y: 1})

        """
        return self.div(other)[0]

    def rem(self, other):
        """Return denominator Factors of ``self / other``.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.rem(b)
        Factors({z: -1})
        >>> a.rem(a)
        Factors({})

        """
        return self.div(other)[1]

    def pow(self, other):
        """Return self raised to a non-negative integer power.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> a**2
        Factors({x: 2, y: 4})

        """
        if isinstance(other, Factors):
            other = other.as_expr()
            if other.is_Integer:
                other = int(other)
        if isinstance(other, numbers.Integral) and other >= 0:
            factors = {}

            if other:
                for factor, exp in self.factors.items():
                    factors[factor] = exp*other

            return Factors(factors)
        else:
            raise ValueError(f'expected non-negative integer, got {other}')

    def gcd(self, other):
        """Return Factors of ``gcd(self, other)``. The keys are
        the intersection of factors with the minimum exponent for
        each factor.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.gcd(b)
        Factors({x: 1, y: 1})

        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return Factors(self.factors)

        factors = {}

        for factor, exp in self.factors.items():
            factor, exp = sympify(factor), sympify(exp)
            if factor in other.factors:
                lt = (exp - other.factors[factor]).is_negative
                if lt:
                    factors[factor] = exp
                elif lt is False:
                    factors[factor] = other.factors[factor]

        return Factors(factors)

    def lcm(self, other):
        """Return Factors of ``lcm(self, other)`` which are
        the union of factors with the maximum exponent for
        each factor.

        Examples
        ========

        >>> a = Factors((x*y**2).as_powers_dict())
        >>> b = Factors((x*y/z).as_powers_dict())
        >>> a.lcm(b)
        Factors({x: 1, y: 2, z: -1})

        """
        if not isinstance(other, Factors):
            other = Factors(other)
            if any(f.is_zero for f in (self, other)):
                return Factors(Integer(0))

        factors = dict(self.factors)

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = max(exp, factors[factor])

            factors[factor] = exp

        return Factors(factors)

    def __mul__(self, other):
        return self.mul(other)

    def __divmod__(self, other):
        return self.div(other)

    def __truediv__(self, other):
        return self.quo(other)

    def __mod__(self, other):
        return self.rem(other)

    def __pow__(self, other):
        return self.pow(other)

    def __eq__(self, other):
        if not isinstance(other, Factors):
            other = Factors(other)
        return self.factors == other.factors


class Term:
    """Efficient representation of ``coeff*(numer/denom)``."""

    def __init__(self, term, numer=None, denom=None):
        if numer is None and denom is None:
            if not term.is_commutative:
                raise NonCommutativeExpression(
                    'commutative expression expected')

            coeff, factors = term.as_coeff_mul()
            numer, denom = defaultdict(int), defaultdict(int)

            for factor in factors:
                base, exp = decompose_power(factor)

                if base.is_Add:
                    cont, base = base.primitive()
                    coeff *= cont**exp

                if exp > 0:
                    numer[base] += exp
                else:
                    denom[base] += -exp

            numer = Factors(numer)
            denom = Factors(denom)
        else:
            coeff = term

            if numer is None:
                numer = Factors()

            if denom is None:
                denom = Factors()

        self.coeff = coeff
        self.numer = numer
        self.denom = denom

    def as_expr(self):
        return self.coeff*(self.numer.as_expr()/self.denom.as_expr())

    def mul(self, other):
        coeff = self.coeff*other.coeff
        numer = self.numer.mul(other.numer)
        denom = self.denom.mul(other.denom)

        numer, denom = numer.normal(denom)

        return Term(coeff, numer, denom)

    def inv(self):
        return Term(1/self.coeff, self.denom, self.numer)

    def quo(self, other):
        return self.mul(other.inv())

    def pow(self, other):
        if other < 0:
            return self.inv().pow(-other)
        else:
            return Term(self.coeff ** other,
                        self.numer.pow(other),
                        self.denom.pow(other))

    def gcd(self, other):
        return Term(self.coeff.gcd(other.coeff),
                    self.numer.gcd(other.numer),
                    self.denom.gcd(other.denom))

    def lcm(self, other):
        return Term(self.coeff.lcm(other.coeff),
                    self.numer.lcm(other.numer),
                    self.denom.lcm(other.denom))

    def __mul__(self, other):
        if isinstance(other, Term):
            return self.mul(other)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Term):
            return self.quo(other)
        else:
            return NotImplemented

    def __pow__(self, other):
        if isinstance(other, numbers.Integral):
            return self.pow(other)
        else:
            return NotImplemented

    def __eq__(self, other):
        return (self.coeff == other.coeff and
                self.numer == other.numer and
                self.denom == other.denom)


def _gcd_terms(terms, isprimitive=False, fraction=True):
    """Helper function for :func:`gcd_terms`.

    If ``isprimitive`` is True then the call to primitive
    for an Add will be skipped. This is useful when the
    content has already been extrated.

    If ``fraction`` is True then the expression will appear over a common
    denominator, the lcm of all term denominators.

    """
    if isinstance(terms, Basic) and not isinstance(terms, Tuple):
        terms = Add.make_args(terms)

    terms = list(map(Term, (t for t in terms if t)))

    # there is some simplification that may happen if we leave this
    # here rather than duplicate it before the mapping of Term onto
    # the terms
    if len(terms) == 0:
        return Integer(0), Integer(0), Integer(1)

    if len(terms) == 1:
        cont = terms[0].coeff
        numer = terms[0].numer.as_expr()
        denom = terms[0].denom.as_expr()

    else:
        cont = terms[0]
        for term in terms[1:]:
            cont = cont.gcd(term)

        for i, term in enumerate(terms):
            terms[i] = term.quo(cont)

        if fraction:
            denom = terms[0].denom

            for term in terms[1:]:
                denom = denom.lcm(term.denom)

            numers = []
            for term in terms:
                numer = term.numer.mul(denom.quo(term.denom))
                numers.append(term.coeff*numer.as_expr())
        else:
            numers = [t.as_expr() for t in terms]
            denom = Term(Integer(1)).numer

        cont = cont.as_expr()
        numer = Add(*numers)
        denom = denom.as_expr()

    if not isprimitive and numer.is_Add:
        _cont, numer = numer.primitive()
        cont *= _cont

    return cont, numer, denom


def gcd_terms(terms, isprimitive=False, clear=True, fraction=True):
    """Compute the GCD of ``terms`` and put them together.

    ``terms`` can be an expression or a non-Basic sequence of expressions
    which will be handled as though they are terms from a sum.

    If ``isprimitive`` is True the _gcd_terms will not run the primitive
    method on the terms.

    ``clear`` controls the removal of integers from the denominator of an Add
    expression. When True (default), all numerical denominator will be cleared;
    when False the denominators will be cleared only if all terms had numerical
    denominators other than 1.

    ``fraction``, when True (default), will put the expression over a common
    denominator.

    Examples
    ========

    >>> gcd_terms((x + 1)**2*y + (x + 1)*y**2)
    y*(x + 1)*(x + y + 1)
    >>> gcd_terms(x/2 + 1)
    (x + 2)/2
    >>> gcd_terms(x/2 + 1, clear=False)
    x/2 + 1
    >>> gcd_terms(x/2 + y/2, clear=False)
    (x + y)/2
    >>> gcd_terms(x/2 + 1/x)
    (x**2 + 2)/(2*x)
    >>> gcd_terms(x/2 + 1/x, fraction=False)
    (x + 2/x)/2
    >>> gcd_terms(x/2 + 1/x, fraction=False, clear=False)
    x/2 + 1/x

    >>> gcd_terms(x/2/y + 1/x/y)
    (x**2 + 2)/(2*x*y)
    >>> gcd_terms(x/2/y + 1/x/y, fraction=False, clear=False)
    (x + 2/x)/(2*y)

    The ``clear`` flag was ignored in this case because the returned
    expression was a rational expression, not a simple sum.

    See Also
    ========
    factor_terms, diofant.polys.polytools.terms_gcd

    """
    def mask(terms):
        """Replace nc portions of each term with a unique Dummy symbols
        and return the replacements to restore them.

        """
        args = [(a, []) if a.is_commutative else a.args_cnc() for a in terms]
        reps = []
        for i, (c, nc) in enumerate(args):
            if nc:
                nc = Mul._from_args(nc)
                d = Dummy()
                reps.append((d, nc))
                c.append(d)
                args[i] = Mul._from_args(c)
            else:
                args[i] = c
        return args, dict(reps)

    isadd = isinstance(terms, Add)
    addlike = isadd or not isinstance(terms, Basic) and \
        is_sequence(terms, include=set) and \
        not isinstance(terms, Dict)

    if addlike:
        if isadd:  # i.e. an Add
            terms = list(terms.args)
        else:
            terms = sympify(terms)
        terms, reps = mask(terms)
        cont, numer, denom = _gcd_terms(terms, isprimitive, fraction)
        numer = numer.xreplace(reps)
        coeff, factors = cont.as_coeff_Mul()
        return _keep_coeff(coeff, factors*numer/denom, clear=clear)

    if not isinstance(terms, Basic):
        return terms

    if terms.is_Atom:
        return terms

    if terms.is_Mul:
        c, args = terms.as_coeff_mul()
        return _keep_coeff(c, Mul(*[gcd_terms(i, isprimitive, clear, fraction)
                                    for i in args]), clear=clear)

    def handle(a):
        # don't treat internal args like terms of an Add
        if not isinstance(a, Expr):
            return a.func(*[handle(i) for i in a.args])
        return gcd_terms(a, isprimitive, clear, fraction)

    if isinstance(terms, Dict):
        return Dict(*[(k, handle(v)) for k, v in terms.args])
    return terms.func(*[handle(i) for i in terms.args])


def factor_terms(expr, radical=False, clear=False, fraction=False, sign=True):
    """Remove common factors from terms in all arguments without
    changing the underlying structure of the expr. No expansion or
    simplification (and no processing of non-commutatives) is performed.

    If radical=True then a radical common to all terms will be factored
    out of any Add sub-expressions of the expr.

    If clear=False (default) then coefficients will not be separated
    from a single Add if they can be distributed to leave one or more
    terms with integer coefficients.

    If fraction=True (default is False) then a common denominator will be
    constructed for the expression.

    If sign=True (default) then even if the only factor in common is a -1,
    it will be factored out of the expression.

    Examples
    ========

    >>> factor_terms(x + x*(2 + 4*y)**3)
    x*(8*(2*y + 1)**3 + 1)
    >>> A = Symbol('A', commutative=False)
    >>> factor_terms(x*A + x*A + x*y*A)
    x*(y*A + 2*A)

    When ``clear`` is False, a rational will only be factored out of an
    Add expression if all terms of the Add have coefficients that are
    fractions:

    >>> factor_terms(x/2 + 1, clear=False)
    x/2 + 1
    >>> factor_terms(x/2 + 1, clear=True)
    (x + 2)/2

    This only applies when there is a single Add that the coefficient
    multiplies:

    >>> factor_terms(x*y/2 + y, clear=True)
    y*(x + 2)/2
    >>> factor_terms(x*y/2 + y, clear=False) == _
    True

    If a -1 is all that can be factored out, to *not* factor it out, the
    flag ``sign`` must be False:

    >>> factor_terms(-x - y)
    -(x + y)
    >>> factor_terms(-x - y, sign=False)
    -x - y
    >>> factor_terms(-2*x - 2*y, sign=False)
    -2*(x + y)

    See Also
    ========
    gcd_terms, diofant.polys.polytools.terms_gcd

    """
    def do(expr):
        is_iterable = iterable(expr)

        if not isinstance(expr, Basic) and is_iterable:
            return type(expr)([do(i) for i in expr])

        if expr.is_Atom:
            return expr

        if expr.is_Pow or expr.is_Function or \
                is_iterable or not hasattr(expr, 'args_cnc'):
            args = expr.args
            newargs = tuple(do(i) for i in args)
            if newargs == args:
                return expr
            return expr.func(*newargs)

        cont, p = expr.as_content_primitive(radical=radical)
        if p.is_Add:
            list_args = [do(a) for a in Add.make_args(p)]
            # get a common negative (if there) which gcd_terms does not remove
            if all(a.as_coeff_Mul()[0] < 0 for a in list_args):
                cont = -cont
                list_args = [-a for a in list_args]
            # watch out for exp(-(x+2)) which gcd_terms will change to exp(-x-2)
            special = {}
            for i, a in enumerate(list_args):
                _, e = a.as_base_exp()
                if e.is_Mul and e != Mul(*e.args):
                    list_args[i] = Dummy()
                    special[list_args[i]] = a
            # rebuild p not worrying about the order which gcd_terms will fix
            p = Add._from_args(list_args)
            p = gcd_terms(p,
                          isprimitive=True,
                          clear=clear,
                          fraction=fraction).xreplace(special)
        elif p.args:
            p = p.func(
                *[do(a) for a in p.args])
        rv = _keep_coeff(cont, p, clear=clear, sign=sign)
        return rv
    expr = sympify(expr)
    return do(expr)


def _mask_nc(eq, name=None):
    """
    Return ``eq`` with non-commutative objects replaced with Dummy
    symbols. A dictionary that can be used to restore the original
    values is returned: if it is None, the expression is noncommutative
    and cannot be made commutative. The third value returned is a list
    of any non-commutative symbols that appear in the returned equation.

    ``name``, if given, is the name that will be used with numered Dummy
    variables that will replace the non-commutative objects and is mainly
    used for doctesting purposes.

    Notes
    =====

    All non-commutative objects other than Symbols are replaced with
    a non-commutative Symbol. Identical objects will be identified
    by identical symbols.

    If there is only 1 non-commutative object in an expression it will
    be replaced with a commutative symbol. Otherwise, the non-commutative
    entities are retained and the calling routine should handle
    replacements in this case since some care must be taken to keep
    track of the ordering of symbols when they occur within Muls.

    Examples
    ========

    >>> A, B, C = symbols('A B C', commutative=False)

    One nc-symbol:

    >>> _mask_nc(A**2 - x**2, 'd')
    (-x**2 + _d0**2, {_d0: A}, [])

    Multiple nc-symbols:

    >>> _mask_nc(A**2 - B**2, 'd')
    (A**2 - B**2, None, [A, B])

    If there is an object that:

        - doesn't contain nc-symbols
        - but has arguments which derive from Expr
        - and doesn't define an _eval_is_commutative routine

    then it will give False (or None?) for the is_commutative test. Such
    objects are also removed by this routine:

    >>> eq = (1 + Mul(Expr(), Expr(), evaluate=False))
    >>> eq.is_commutative is None
    True
    >>> _mask_nc(eq, 'd')
    (_d0**2 + 1, {_d0: Expr()}, [])

    """
    name = name or 'mask'
    # Make Dummy() append sequential numbers to the name

    def numbered_names():
        i = 0
        while True:
            yield name + str(i)
            i += 1

    names = numbered_names()

    def Dummy(*args, **kwargs):
        from .symbol import Dummy
        return Dummy(next(names), *args, **kwargs)

    expr = eq
    if expr.is_commutative:
        return eq, {}, []

    # identify nc-objects; symbols and other
    rep = []
    nc_obj = set()
    nc_syms = set()
    pot = preorder_traversal(expr, keys=default_sort_key)
    for a in pot:
        if any(a == r[0] for r in rep):
            pot.skip()
        elif not a.is_commutative:
            if a.is_Symbol:
                nc_syms.add(a)
            elif not (a.is_Add or a.is_Mul or a.is_Pow):
                if all(s.is_commutative for s in a.free_symbols):
                    rep.append((a, Dummy()))
                else:
                    nc_obj.add(a)
                pot.skip()

    # If there is only one nc symbol, it can be factored regularly
    # but polys is going to complain, so replace it with a Dummy.
    if len(nc_syms) == 1 and not nc_obj:
        rep.append((nc_syms.pop(), Dummy()))

    # Any remaining nc-objects will be replaced with an nc-Dummy and
    # identified as an nc-Symbol to watch out for
    nc_obj = sorted(nc_obj, key=default_sort_key)
    for n in nc_obj:
        nc = Dummy(commutative=False)
        rep.append((n, nc))
        nc_syms.add(nc)
    expr = expr.subs(rep)

    nc_syms = list(nc_syms)
    nc_syms.sort(key=default_sort_key)
    return expr, {v: k for k, v in rep} or None, nc_syms


def factor_nc(expr):
    """Return the factored form of ``expr`` while handling non-commutative
    expressions.

    Examples
    ========

    >>> A = Symbol('A', commutative=False)
    >>> B = Symbol('B', commutative=False)
    >>> factor_nc((x**2 + 2*A*x + A**2).expand())
    (x + A)**2
    >>> factor_nc(((x + A)*(x + B)).expand())
    (x + A)*(x + B)

    """
    from ..polys import factor, gcd
    from ..simplify.simplify import powsimp

    def _pemexpand(expr):
        """Expand with the minimal set of hints necessary to check the result."""
        return expr.expand(deep=True, mul=True, power_exp=True,
                           power_base=False, basic=False, multinomial=True, log=False)

    expr = sympify(expr)
    if not isinstance(expr, Expr) or not expr.args:
        return expr
    if not expr.is_Add:
        return expr.func(*[factor_nc(a) for a in expr.args])

    expr, rep, nc_symbols = _mask_nc(expr)
    if rep:
        return factor(expr).subs(rep)
    else:
        args = [a.args_cnc() for a in Add.make_args(expr)]
        c = g = l = r = Integer(1)
        hit = False
        # find any commutative gcd term
        for i, a in enumerate(args):
            if i == 0:
                c = Mul._from_args(a[0])
            elif a[0]:
                c = gcd(c, Mul._from_args(a[0]))
            else:
                c = Integer(1)
        if c != 1:
            hit = True
            c, g = c.as_coeff_Mul()
            if g != 1:
                for i, (cc, _) in enumerate(args):
                    cc = list(Mul.make_args(Mul._from_args(list(cc))/g))
                    args[i][0] = cc
            for i, (cc, _) in enumerate(args):
                cc[0] = cc[0]/c
                args[i][0] = cc
        # find any noncommutative common prefix
        for i, a in enumerate(args):
            if i == 0:
                n = a[1][:]
            else:
                n = common_prefix(n, a[1])
            if not n:
                # is there a power that can be extracted?
                if not args[0][1]:
                    break
                b, e = args[0][1][0].as_base_exp()
                ok = False
                if e.is_Integer:
                    for t in args:
                        if not t[1]:
                            break
                        bt, et = t[1][0].as_base_exp()
                        if et.is_Integer and bt == b:
                            e = min(e, et)
                        else:
                            break
                    else:
                        ok = hit = True
                        l = b**e
                        il = b**-e
                        for i, a in enumerate(args):
                            args[i][1][0] = il*args[i][1][0]
                        break
                break
        else:
            hit = True
            lenn = len(n)
            l = Mul(*n)
            for i, a in enumerate(args):
                args[i][1] = args[i][1][lenn:]
        # find any noncommutative common suffix
        for i, a in enumerate(args):
            if i == 0:
                n = a[1][:]
            else:
                n = common_suffix(n, a[1])
            if not n:
                # is there a power that can be extracted?
                if not args[0][1]:
                    break
                b, e = args[0][1][-1].as_base_exp()
                ok = False
                if e.is_Integer:
                    for t in args:
                        if not t[1]:
                            break
                        bt, et = t[1][-1].as_base_exp()
                        if et.is_Integer and bt == b:
                            e = min(e, et)
                        else:
                            break
                    else:
                        ok = hit = True
                        r = b**e
                        il = b**-e
                        for i, a in enumerate(args):
                            args[i][1][-1] = args[i][1][-1]*il
                        break
                break
        else:
            hit = True
            lenn = len(n)
            r = Mul(*n)
            for i, a in enumerate(args):
                args[i][1] = a[1][:len(a[1]) - lenn]
        if hit:
            mid = Add(*[Mul(*cc)*Mul(*nc) for cc, nc in args])
        else:
            mid = expr

        # sort the symbols so the Dummys would appear in the same
        # order as the original symbols, otherwise you may introduce
        # a factor of -1, e.g. A**2 - B**2) -- {A:y, B:x} --> y**2 - x**2
        # and the former factors into two terms, (A - B)*(A + B) while the
        # latter factors into 3 terms, (-1)*(x - y)*(x + y)
        rep1 = [(n, Dummy()) for n in sorted(nc_symbols, key=default_sort_key)]
        unrep1 = [(v, k) for k, v in rep1]
        unrep1.reverse()
        new_mid, r2, _ = _mask_nc(mid.subs(rep1))
        new_mid = powsimp(factor(new_mid))

        new_mid = new_mid.subs(r2).subs(unrep1)

        if new_mid.is_Pow:
            return _keep_coeff(c, g*l*new_mid*r)

        if new_mid.is_Mul:
            # XXX TODO there should be a way to inspect what order the terms
            # must be in and just select the plausible ordering without
            # checking permutations
            cfac = []
            ncfac = []
            for f in new_mid.args:
                if f.is_commutative:
                    cfac.append(f)
                else:
                    b, e = f.as_base_exp()
                    if e.is_Integer:
                        ncfac.extend([b]*e)
                    else:
                        ncfac.append(f)
            pre_mid = g*Mul(*cfac)*l
            target = _pemexpand(expr/c)
            for s in variations(ncfac, len(ncfac)):
                ok = pre_mid*Mul(*s)*r
                if _pemexpand(ok) == target:
                    return _keep_coeff(c, ok)

        # mid was an Add that didn't factor successfully
        return _keep_coeff(c, g*l*mid*r)
