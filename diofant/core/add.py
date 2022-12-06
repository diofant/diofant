import collections
import functools
import math

from ..utilities import default_sort_key
from .cache import cacheit
from .compatibility import is_sequence
from .logic import _fuzzy_group
from .numbers import Integer, nan, oo, zoo
from .operations import AssocOp


class Add(AssocOp):
    """Symbolic addition class."""

    is_Add = True

    identity = Integer(0)

    @classmethod
    def flatten(cls, seq):
        """
        Takes the sequence "seq" of nested Adds and returns a flatten list.

        Returns: (commutative_part, noncommutative_part, order_symbols)

        Applies associativity, all terms are commutable with respect to
        addition.

        See also
        ========

        diofant.core.mul.Mul.flatten

        """
        from ..series.order import Order
        from .mul import Mul

        rv = None
        if len(seq) == 2:
            a, b = seq
            if b.is_Rational:
                a, b = b, a
            if a.is_Rational:
                if b.is_Mul:
                    rv = [a, b], [], None
            if rv:
                if all(s.is_commutative for s in rv[0]):
                    return rv
                return [], rv[0], None

        # term -> coeff
        # e.g. x**2 -> 5   for ... + 5*x**2 + ...
        terms = {}

        # coefficient (Number or zoo) to always be in slot 0, e.g. 3 + ...
        coeff = Integer(0)

        order_factors = []

        for o in seq:

            # O(x)
            if o.is_Order:
                for o1 in order_factors:
                    if o1.contains(o):
                        o = None
                        break
                if o is None:
                    continue
                order_factors = [o] + [
                    o1 for o1 in order_factors if not o.contains(o1)]
                continue

            # 3 or NaN
            if o.is_Number:
                if (o is nan or coeff is zoo and
                        o.is_finite is False):
                    # we know for sure the result will be nan
                    return [nan], [], None
                if coeff.is_Number:
                    coeff += o
                    if coeff is nan:
                        # we know for sure the result will be nan
                        return [nan], [], None
                continue

            if o is zoo:
                if coeff.is_finite is False:
                    # we know for sure the result will be nan
                    return [nan], [], None
                coeff = zoo
                continue

            # Add([...])
            if o.is_Add:
                # NB: here we assume Add is always commutative
                seq.extend(o.args)  # TODO zerocopy?
                continue

            if o.has(Order):
                c, s = Integer(1), o

            # Mul([...])
            elif o.is_Mul:
                c, s = o.as_coeff_Mul()

            # check for unevaluated Pow, e.g. 2**3 or 2**(-1/2)
            elif o.is_Pow:
                b, e = o.as_base_exp()
                if b.is_Number and (e.is_Integer or
                                    (e.is_Rational and e.is_negative)):
                    seq.append(b**e)
                    continue
                c, s = Integer(1), o

            else:
                # everything else
                c = Integer(1)
                s = o

            # now we have:
            # o = c*s, where
            #
            # c is a Number
            # s is an expression with number factor extracted
            # let's collect terms with the same s, so e.g.
            # 2*x**2 + 3*x**2  ->  5*x**2
            if s in terms:
                terms[s] += c
                if terms[s] is nan:
                    # we know for sure the result will be nan
                    return [nan], [], None
            else:
                terms[s] = c

        # now let's construct new args:
        # [2*x**2, x**3, 7*x**4, pi, ...]
        newseq = []
        noncommutative = False
        for s, c in terms.items():
            # 0*s
            if c == 0:
                continue
            # 1*s
            if c == 1 and not c.is_Float:
                newseq.append(s)
            # c*s
            else:
                if s.is_Mul:
                    # Mul, already keeps its arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = s._new_rawargs(*((c,) + s.args))
                    newseq.append(cs)
                elif s.is_Add:
                    # we just re-create the unevaluated Mul
                    newseq.append(Mul(c, s, evaluate=False))
                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(Mul(c, s))

            noncommutative = noncommutative or not s.is_commutative

        # oo, -oo
        if coeff is oo:
            newseq = [f for f in newseq if not
                      (f.is_nonnegative or f.is_extended_real and f.is_finite)]

        elif coeff == -oo:
            newseq = [f for f in newseq if not
                      (f.is_nonpositive or f.is_extended_real and f.is_finite)]

        if coeff is zoo:
            # zoo might be
            #   infinite_real + finite_im
            #   finite_real + infinite_im
            #   infinite_real + infinite_im
            # addition of a finite real or imaginary number won't be able to
            # change the zoo nature; adding an infinite qualtity would result
            # in a NaN condition if it had sign opposite of the infinite
            # portion of zoo, e.g., infinite_real - infinite_real.
            newseq = [c for c in newseq if not (c.is_finite and
                                                (c.is_extended_real is not None
                                                 or c.is_imaginary is not None))]

        # process O(x)
        if order_factors:
            newseq2 = []
            for t in newseq:
                for o in order_factors:
                    # x + O(x) -> O(x)
                    if o.contains(t):
                        t = None
                        break
                # x + O(x**2) -> x + O(x**2)
                if t is not None:
                    newseq2.append(t)
            newseq = newseq2 + order_factors
            # 1 + O(1) -> O(1)
            for o in order_factors:
                if o.contains(coeff):
                    coeff = Integer(0)
                    break

        # order args canonically
        newseq.sort(key=default_sort_key)

        # current code expects coeff to be first
        if coeff != 0:
            newseq.insert(0, coeff)

        # we are done
        if noncommutative:
            return [], newseq, None
        else:
            return newseq, [], None

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 4, 1, cls.__name__

    def as_coefficients_dict(self):
        """Return a dictionary mapping terms to their Rational coefficient.

        Since the dictionary is a defaultdict, inquiries about terms which
        were not present will return a coefficient of 0. If an expression is
        not an Add it is considered to have a single term.

        Examples
        ========

        >>> (3*x + x*y + 4).as_coefficients_dict()
        {1: 4, x: 3, x*y: 1}
        >>> _[y]
        0
        >>> (3*y*x).as_coefficients_dict()
        {x*y: 3}

        """
        d = collections.defaultdict(list)
        for ai in self.args:
            c, m = ai.as_coeff_Mul()
            d[m].append(c)
        for k, v in d.items():
            if len(v) == 1:
                d[k] = v[0]
            else:
                d[k] = Add(*v)
        di = collections.defaultdict(int)
        di.update(d)
        return di

    @cacheit
    def as_coeff_add(self, *deps):
        """
        Returns a tuple (coeff, args) where self is treated as an Add and coeff
        is the Number term and args is a tuple of all other terms.

        Examples
        ========

        >>> (7 + 3*x).as_coeff_add()
        (7, (3*x,))
        >>> (7*x).as_coeff_add()
        (0, (7*x,))

        """
        if deps:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(*deps):
                    l2.append(f)
                else:
                    l1.append(f)
            return self._new_rawargs(*l1), tuple(l2)
        coeff, notrat = self.args[0].as_coeff_add()
        if coeff != 0:
            return coeff, notrat + self.args[1:]
        return Integer(0), self.args

    def as_coeff_Add(self, rational=False):
        """Efficiently extract the coefficient of a summation."""
        coeff, args = self.args[0], self.args[1:]

        if coeff.is_Number and not rational or coeff.is_Rational:
            return coeff, self._new_rawargs(*args)
        return Integer(0), self

    # Note, we intentionally do not implement Add.as_coeff_mul().  Rather, we
    # let Expr.as_coeff_mul() just always return (Integer(1), self) for an Add.  See
    # issue sympy/sympy#5524.

    @cacheit
    def _eval_derivative(self, s):
        return self.func(*[a.diff(s) for a in self.args])

    def _eval_nseries(self, x, n, logx):
        terms = [t.nseries(x, n=n, logx=logx) for t in self.args]
        return self.func(*terms)

    def _matches_simple(self, expr, repl_dict):
        # handle (w+3)._matches('x+5') -> {w: x+2}
        coeff, terms = self.as_coeff_add()
        if len(terms) == 1:
            return terms[0]._matches(expr - coeff, repl_dict)

    def _matches(self, expr, repl_dict={}):
        """Helper method for match().

        See Also
        ========

        diofant.core.basic.Basic.matches

        """
        return AssocOp._matches_commutative(self, expr, repl_dict)

    @staticmethod
    def _combine_inverse(lhs, rhs):
        """
        Returns lhs - rhs, but treats arguments like symbols, so things like
        oo - oo return 0, instead of a nan.

        """
        from ..simplify import signsimp
        from . import I, oo
        if lhs == oo and rhs == oo or lhs == oo*I and rhs == oo*I:
            return Integer(0)
        return signsimp(lhs - rhs)

    @cacheit
    def as_two_terms(self):
        """Return head and tail of self.

        This is the most efficient way to get the head and tail of an
        expression.

        - if you want only the head, use self.args[0];
        - if you want to process the arguments of the tail then use
          self.as_coef_add() which gives the head and a tuple containing
          the arguments of the tail when treated as an Add.
        - if you want the coefficient when self is treated as a Mul
          then use self.as_coeff_mul()[0]

        >>> (3*x*y).as_two_terms()
        (3, x*y)

        """
        return self.args[0], self._new_rawargs(*self.args[1:])

    def _eval_as_numer_denom(self):
        """Expression -> a/b -> a, b.

        See Also
        ========

        diofant.core.expr.Expr.as_numer_denom

        """
        from .mul import Mul, _keep_coeff

        # clear rational denominator
        content, expr = self.primitive()
        if not isinstance(expr, self.func):
            return Mul(content, expr, evaluate=False).as_numer_denom()
        ncon, dcon = content.as_numer_denom()

        # collect numerators and denominators of the terms
        nd = collections.defaultdict(list)
        for f in expr.args:
            ni, di = f.as_numer_denom()
            nd[di].append(ni)

        # check for quick exit
        if len(nd) == 1:
            d, n = nd.popitem()
            return self.func(
                *[_keep_coeff(ncon, ni) for ni in n]), _keep_coeff(dcon, d)

        # sum up the terms having a common denominator
        for d, n in nd.items():
            if len(n) == 1:
                nd[d] = n[0]
            else:
                nd[d] = self.func(*n)

        # assemble single numerator and denominator
        denoms, numers = map(list, zip(*nd.items()))
        n, d = self.func(*[Mul(*(denoms[:i] + [numers[i]] + denoms[i + 1:]))
                           for i in range(len(numers))]), Mul(*denoms)

        return _keep_coeff(ncon, n), _keep_coeff(dcon, d)

    def _eval_is_polynomial(self, syms):
        return all(term._eval_is_polynomial(syms) for term in self.args)

    def _eval_is_rational_function(self, syms):
        return all(term._eval_is_rational_function(syms) for term in self.args)

    def _eval_is_algebraic_expr(self, syms):
        return all(term._eval_is_algebraic_expr(syms) for term in self.args)

    # assumption methods

    def _eval_is_commutative(self):
        return _fuzzy_group((a.is_commutative for a in self.args),
                            quick_exit=True)

    def _eval_is_real(self):
        return _fuzzy_group((a.is_real for a in self.args), quick_exit=True)

    def _eval_is_extended_real(self):
        r = _fuzzy_group((a.is_extended_real for a in self.args), quick_exit=True)
        if r is not True:
            return r
        else:
            nfin = [_ for _ in self.args if not _.is_finite]
            if len(nfin) <= 1:
                return True
            elif (all(_.is_nonnegative for _ in nfin) or
                  all(_.is_nonpositive for _ in nfin)):
                return True

    def _eval_is_complex(self):
        return _fuzzy_group((a.is_complex for a in self.args), quick_exit=True)

    def _eval_is_finite(self):
        return _fuzzy_group((a.is_finite for a in self.args), quick_exit=True)

    def _eval_is_integer(self):
        return _fuzzy_group((a.is_integer for a in self.args), quick_exit=True)

    def _eval_is_rational(self):
        return _fuzzy_group((a.is_rational for a in self.args), quick_exit=True)

    def _eval_is_algebraic(self):
        return _fuzzy_group((a.is_algebraic for a in self.args), quick_exit=True)

    def _eval_is_imaginary(self):
        return _fuzzy_group((a.is_imaginary for a in self.args), quick_exit=True)

    def _eval_is_odd(self):
        l = [f for f in self.args if not f.is_even]
        if not l:
            return False
        if l[0].is_odd:
            return self._new_rawargs(*l[1:]).is_even

    def _eval_is_irrational(self):
        for t in self.args:
            a = t.is_irrational
            if a:
                if all(x.is_rational for x in self.args if x != t):
                    return True
                return
            if a is None:
                return
        return False

    def _eval_is_positive(self):
        if self.is_number:
            n = super()._eval_is_positive()
            if n is not None:
                return n

        if any(a.is_infinite for a in self.args):
            args = [a for a in self.args if not a.is_finite]
        else:
            args = self.args

        nonpos = nonneg = 0
        for a in args:
            if a.is_positive:
                continue
            if a.is_nonnegative:
                nonneg += 1
                if a.is_zero:
                    nonpos += 1
            elif a.is_nonpositive:
                nonpos += 1
            else:
                break
        else:
            if not nonpos and nonneg < len(args):
                return True
            elif nonpos == len(args):
                return False

    def _eval_is_negative(self):
        if self.is_number:
            n = super()._eval_is_negative()
            if n is not None:
                return n

        if any(a.is_infinite for a in self.args):
            args = [a for a in self.args if not a.is_finite]
        else:
            args = self.args

        nonneg = nonpos = 0
        for a in args:
            if a.is_negative:
                continue
            if a.is_nonpositive:
                nonpos += 1
                if a.is_zero:
                    nonneg += 1
            elif a.is_nonnegative:
                nonneg += 1
            else:
                break
        else:
            if not nonneg and nonpos < len(args):
                return True
            elif nonneg == len(args):
                return False

    def _eval_subs(self, old, new):
        if not old.is_Add:
            return

        coeff_self, terms_self = self.as_coeff_Add()
        coeff_old, terms_old = old.as_coeff_Add()

        if coeff_self.is_Rational and coeff_old.is_Rational:
            if terms_self == terms_old:   # (2 + a).subs({+3 + a: y}) -> -1 + y
                return self.func(new, coeff_self, -coeff_old)
            if terms_self == -terms_old:  # (2 + a).subs({-3 - a: y}) -> -1 - y
                return self.func(-new, coeff_self, coeff_old)

        if coeff_self.is_Rational and coeff_old.is_Rational \
                or coeff_self == coeff_old:
            args_old, args_self = self.func.make_args(
                terms_old), self.func.make_args(terms_self)
            if len(args_old) < len(args_self):  # (a+b+c).subs({b+c: x}) -> a+x
                self_set = set(args_self)
                old_set = set(args_old)

                if old_set < self_set:
                    ret_set = self_set - old_set
                    return self.func(new, coeff_self, -coeff_old,
                                     *[s._subs(old, new) for s in ret_set])

                args_old = self.func.make_args(
                    -terms_old)     # (a+b+c+d).subs({-b-c: x}) -> a-x+d
                old_set = set(args_old)
                if old_set < self_set:
                    ret_set = self_set - old_set
                    return self.func(-new, coeff_self, coeff_old,
                                     *[s._subs(old, new) for s in ret_set])

    def removeO(self):
        """Removes the additive O(..) symbol.

        See Also
        ========

        diofant.core.expr.Expr.removeO

        """
        args = [a for a in self.args if not a.is_Order]
        return self._new_rawargs(*args)

    def getO(self):
        """Returns the additive O(..) symbol.

        See Also
        ========

        diofant.core.expr.Expr.getO

        """
        args = [a for a in self.args if a.is_Order]
        if args:
            return self._new_rawargs(*args)

    @cacheit
    def extract_leading_order(self, symbols):
        """Returns the leading term and its order.

        Examples
        ========

        >>> (x + 1 + 1/x**5).extract_leading_order(x)
        ((x**(-5), O(x**(-5))),)
        >>> (1 + x).extract_leading_order(x)
        ((1, O(1)),)
        >>> (x + x**2).extract_leading_order(x)
        ((x, O(x)),)

        """
        from ..series import Order
        lst = []
        symbols = list(symbols if is_sequence(symbols) else [symbols])
        point = [0]*len(symbols)
        seq = [(f, Order(f, *zip(symbols, point))) for f in self.args]
        for ef, of in seq:
            for e, o in lst:
                if o.contains(of) and o != of:
                    of = None
                    break
            if of is None:
                continue
            new_lst = [(ef, of)]
            for e, o in lst:
                if of.contains(o) and o != of:
                    continue
                new_lst.append((e, o))
            lst = new_lst
        return tuple(lst)

    def as_real_imag(self, deep=True, **hints):
        """
        Returns a tuple representing a complex number.

        Examples
        ========

        >>> (7 + 9*I).as_real_imag()
        (7, 9)
        >>> ((1 + I)/(1 - I)).as_real_imag()
        (0, 1)
        >>> ((1 + 2*I)*(1 + 3*I)).as_real_imag()
        (-5, 5)

        """
        sargs = self.args
        re_part, im_part = [], []
        for term in sargs:
            re, im = term.as_real_imag(deep=deep)
            re_part.append(re)
            im_part.append(im)
        return self.func(*re_part), self.func(*im_part)

    def _eval_as_leading_term(self, x):
        from ..series import Order
        from . import factor_terms

        by_O = functools.cmp_to_key(lambda f, g: 1 if Order(g, x).contains(f) is not False else -1)
        expr = Integer(0)

        for t in sorted((_.as_leading_term(x) for _ in self.args), key=by_O):
            expr += t
            if not expr:
                # simple leading term analysis gave us 0 but we have to send
                # back a term, so compute the leading term (via series)
                return self.compute_leading_term(x)

        expr = expr.removeO()

        if not expr.is_Add:
            return expr
        else:
            plain = expr.func(*[s for s, _ in expr.extract_leading_order(x)])
            rv = factor_terms(plain, fraction=False)
            rv_simplify = rv.simplify()
            # if it simplifies to an x-free expression, return that;
            # tests don't fail if we don't but it seems nicer to do this
            if x not in rv_simplify.free_symbols:
                if rv_simplify.is_zero and plain:
                    return (expr - plain)._eval_as_leading_term(x)
                return rv_simplify
            return rv

    def _eval_adjoint(self):
        return self.func(*[t.adjoint() for t in self.args])

    def _eval_conjugate(self):
        return self.func(*[t.conjugate() for t in self.args])

    def _eval_transpose(self):
        return self.func(*[t.transpose() for t in self.args])

    def __neg__(self):
        return self.func(*[-t for t in self.args])

    def primitive(self):
        """
        Return ``(R, self/R)`` where ``R`` is the Rational GCD of ``self``.

        ``R`` is collected only from the leading coefficient of each term.

        Examples
        ========

        >>> (2*x + 4*y).primitive()
        (2, x + 2*y)

        >>> (2*x/3 + 4*y/9).primitive()
        (2/9, 3*x + 2*y)

        >>> (2*x/3 + 4.2*y).primitive()
        (1/3, 2*x + 12.6*y)

        No subprocessing of term factors is performed:

        >>> ((2 + 2*x)*x + 2).primitive()
        (1, x*(2*x + 2) + 2)

        Recursive subprocessing can be done with the as_content_primitive()
        method:

        >>> ((2 + 2*x)*x + 2).as_content_primitive()
        (2, x*(x + 1) + 1)

        See Also
        ========

        diofant.polys.polytools.primitive

        """
        from .mul import _keep_coeff
        from .numbers import Rational

        terms = []
        inf = False
        for a in self.args:
            c, m = a.as_coeff_Mul()
            if not c.is_Rational:
                c = Integer(1)
                m = a
            inf = inf or m is zoo
            terms.append((c.numerator, c.denominator, m))

        if not inf:
            ngcd = functools.reduce(math.gcd, [t[0] for t in terms], 0)
            dlcm = functools.reduce(math.lcm, [t[1] for t in terms], 1)
        else:
            ngcd = functools.reduce(math.gcd, [t[0] for t in terms if t[1]], 0)
            dlcm = functools.reduce(math.lcm, [t[1] for t in terms if t[1]], 1)

        if ngcd == dlcm == 1:
            return Integer(1), self
        for i, (p, q, term) in enumerate(terms):
            terms[i] = _keep_coeff(Rational((p//ngcd)*(dlcm//q)), term)

        # we don't need a complete re-flattening since no new terms will join
        # so we just use the same sort as is used in Add.flatten. When the
        # coefficient changes, the ordering of terms may change, e.g.
        #     (3*x, 6*y) -> (2*y, x)
        #
        # We do need to make sure that term[0] stays in position 0, however.
        #
        if terms[0].is_Number or terms[0] is zoo:
            c = terms.pop(0)
        else:
            c = None
        terms.sort(key=default_sort_key)
        if c:
            terms.insert(0, c)
        return Rational(ngcd, dlcm), self._new_rawargs(*terms)

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self. If radical is True (default is False) then
        common radicals will be removed and included as a factor of the
        primitive expression.

        Examples
        ========

        >>> (3 + 3*sqrt(2)).as_content_primitive()
        (3, 1 + sqrt(2))

        Radical content can also be factored out of the primitive:

        >>> (2*sqrt(2) + 4*sqrt(10)).as_content_primitive(radical=True)
        (2, sqrt(2)*(1 + 2*sqrt(5)))

        See Also
        ========

        diofant.core.expr.Expr.as_content_primitive

        """
        from ..functions import root
        from .mul import Mul, _keep_coeff

        con, prim = self.func(*[_keep_coeff(*a.as_content_primitive(
            radical=radical)) for a in self.args]).primitive()
        if radical and prim.is_Add:
            # look for common radicals that can be removed
            args = prim.args
            rads = []
            common_q = None
            for m in args:
                term_rads = collections.defaultdict(list)
                for ai in Mul.make_args(m):
                    if ai.is_Pow:
                        b, e = ai.as_base_exp()
                        if e.is_Rational and b.is_Integer:
                            term_rads[e.denominator].append(abs(int(b))**e.numerator)
                if not term_rads:
                    break
                if common_q is None:
                    common_q = set(term_rads)
                else:
                    common_q = common_q & set(term_rads)
                    if not common_q:
                        break
                rads.append(term_rads)
            else:
                # process rads
                # keep only those in common_q
                for r in rads:
                    for q in list(r):
                        if q not in common_q:
                            r.pop(q)
                    for q in r:
                        r[q] = math.prod(r[q])
                # find the gcd of bases for each q
                G = []
                for q in common_q:
                    g = functools.reduce(math.gcd, [r[q] for r in rads], 0)
                    if g != 1:
                        G.append(root(g, q))
                if G:
                    G = Mul(*G)
                    args = [ai/G for ai in args]
                    prim = G*prim.func(*args)

        return con, prim

    @property
    def _sorted_args(self):
        return tuple(sorted(self.args, key=default_sort_key))
