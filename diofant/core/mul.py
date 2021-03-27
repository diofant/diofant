from collections import defaultdict

from ..utilities import default_sort_key
from .basic import Basic
from .cache import cacheit
from .logic import _fuzzy_group, fuzzy_and
from .operations import AssocOp
from .singleton import S
from .sympify import sympify


# internal marker to indicate:
#   "there are still non-commutative objects -- don't forget to process them"


class NC_Marker:
    """Helper class to mark non-commutative exponents."""

    is_Order = False
    is_Mul = False
    is_Number = False
    is_Poly = False

    is_commutative = False


def _unevaluated_Mul(*args):
    """Return a well-formed unevaluated Mul: Numbers are collected and
    put in slot 0, any arguments that are Muls will be flattened, and args
    are sorted. Use this when args have changed but you still want to return
    an unevaluated Mul.

    Examples
    ========

    >>> a = _unevaluated_Mul(*[Float(3.0), x, Integer(2)])
    >>> a.args[0]
    6.00000000000000
    >>> a.args[1]
    x

    Two unevaluated Muls with the same arguments will
    always compare as equal during testing:

    >>> m = _unevaluated_Mul(sqrt(2), sqrt(3))
    >>> m == _unevaluated_Mul(sqrt(3), sqrt(2))
    True
    >>> u = Mul(sqrt(3), sqrt(2), evaluate=False)
    >>> m == _unevaluated_Mul(u)
    True
    >>> m == Mul(*m.args)
    False

    """
    args = list(args)
    newargs = []
    ncargs = []
    co = S.One
    while args:
        a = args.pop()
        if a.is_Mul:
            c, nc = a.args_cnc()
            args.extend(c)
            if nc:
                ncargs.append(Mul._from_args(nc))
        elif a.is_Number:
            co *= a
        else:
            newargs.append(a)
    newargs.sort(key=default_sort_key)
    if co is not S.One:
        newargs.insert(0, co)
    if ncargs:
        newargs.append(Mul._from_args(ncargs))
    return Mul._from_args(newargs)


class Mul(AssocOp):
    """Symbolic multiplication class."""

    is_Mul = True

    identity = S.One

    @classmethod
    def flatten(cls, seq):
        """Return commutative, noncommutative and order arguments by
        combining related terms.

        Notes
        =====
            * In an expression like ``a*b*c``, python process this through diofant
              as ``Mul(Mul(a, b), c)``. This can have undesirable consequences.

              -  Sometimes terms are not combined as one would like:
                 {c.f. https://github.com/sympy/sympy/issues/4596}

                >>> 2*(x + 1)  # this is the 2-arg Mul behavior
                2*x + 2
                >>> y*(x + 1)*2
                2*y*(x + 1)
                >>> 2*(x + 1)*y  # 2-arg result will be obtained first
                y*(2*x + 2)
                >>> Mul(2, x + 1, y)  # all 3 args simultaneously processed
                2*y*(x + 1)
                >>> 2*((x + 1)*y)  # parentheses can control this behavior
                2*y*(x + 1)

                Powers with compound bases may not find a single base to
                combine with unless all arguments are processed at once.
                Post-processing may be necessary in such cases.
                {c.f. https://github.com/sympy/sympy/issues/5728}

                >>> a = sqrt(x*sqrt(y))
                >>> a**3
                (x*sqrt(y))**(3/2)
                >>> Mul(a, a, a)
                (x*sqrt(y))**(3/2)
                >>> a*a*a
                x*sqrt(y)*sqrt(x*sqrt(y))
                >>> _.subs({a.base: z}).subs({z: a.base})
                (x*sqrt(y))**(3/2)

              -  If more than two terms are being multiplied then all the
                 previous terms will be re-processed for each new argument.
                 So if each of ``a``, ``b`` and ``c`` were :class:`Mul`
                 expression, then ``a*b*c`` (or building up the product
                 with ``*=``) will process all the arguments of ``a`` and
                 ``b`` twice: once when ``a*b`` is computed and again when
                 ``c`` is multiplied.

                 Using ``Mul(a, b, c)`` will process all arguments once.

            * The results of Mul are cached according to arguments, so flatten
              will only be called once for ``Mul(a, b, c)``. If you can
              structure a calculation so the arguments are most likely to be
              repeats then this can save time in computing the answer. For
              example, say you had a Mul, M, that you wished to divide by ``d[i]``
              and multiply by ``n[i]`` and you suspect there are many repeats
              in ``n``. It would be better to compute ``M*n[i]/d[i]`` rather
              than ``M/d[i]*n[i]`` since every time n[i] is a repeat, the
              product, ``M*n[i]`` will be returned without flattening -- the
              cached value will be returned. If you divide by the ``d[i]``
              first (and those are more unique than the ``n[i]``) then that will
              create a new Mul, ``M/d[i]`` the args of which will be traversed
              again when it is multiplied by ``n[i]``.

              {c.f. https://github.com/sympy/sympy/issues/5706}

              This consideration is moot if the cache is turned off.

            The validity of the above notes depends on the implementation
            details of Mul and flatten which may change at any time. Therefore,
            you should only consider them when your code is highly performance
            sensitive.

        """
        from ..series.order import Order

        rv = None
        if len(seq) == 2:
            a, b = seq
            if b.is_Rational:
                a, b = b, a
            assert a is not S.One
            if not a.is_zero and a.is_Rational:
                r, b = b.as_coeff_Mul()
                if b.is_Add:
                    if r is not S.One:  # 2-arg hack
                        if a*r is S.One:
                            rv = [b], [], None
                        else:
                            # leave the Mul as a Mul
                            rv = [cls(a*r, b, evaluate=False)], [], None
                    elif b.is_commutative:
                        r, b = b.as_coeff_Add()
                        bargs = [_keep_coeff(a, bi) for bi in Add.make_args(b)]
                        bargs.sort(key=default_sort_key)
                        ar = a*r
                        if ar:
                            bargs.insert(0, ar)
                        bargs = [Add._from_args(bargs)]
                        rv = bargs, [], None
            if rv:
                return rv

        # apply associativity, separate commutative part of seq
        c_part = []         # out: commutative factors
        nc_part = []        # out: non-commutative factors

        nc_seq = []

        # standalone term e.g. 3 * ...
        coeff = S.One

        #                           n
        # (base,exp)e.g. (x,n) for x
        c_powers = []

        #                                           y
        # (num-base, exp) e.g.  (3, y)  for  ... * 3  * ...
        num_exp = []

        neg1e = S.Zero  # exponent on -1 extracted from Number-based Pow and I

        #                                                 1/2
        # (num-base, Rat-exp) e.g.  (3, 1/2)  for  ... * 3     * ...
        pnum_rat = defaultdict(list)

        order_symbols = None

        # --- PART 1 ---
        #
        # "collect powers and coeff":
        #
        # o coeff
        # o c_powers
        # o num_exp
        # o neg1e
        # o pnum_rat
        #
        # NOTE: this is optimized for all-objects-are-commutative case
        for o in seq:
            # O(x)
            if o.is_Order:
                o, order_symbols = o.as_expr_variables(order_symbols)

            # Mul([...])
            if o.is_Mul:
                if o.is_commutative:
                    seq.extend(o.args)    # XXX zerocopy?

                else:
                    # NCMul can have commutative parts as well
                    for q in o.args:
                        if q.is_commutative:
                            seq.append(q)
                        else:
                            nc_seq.append(q)

                    # append non-commutative marker, so we don't forget to
                    # process scheduled non-commutative objects
                    seq.append(NC_Marker)

                continue

            # 3
            elif o.is_Number:
                if o is nan or coeff is zoo and o is S.Zero:
                    # we know for sure the result will be nan
                    return [nan], [], None
                if coeff.is_Number:  # it could be zoo
                    coeff *= o
                    if coeff is nan:
                        # we know for sure the result will be nan
                        return [nan], [], None
                o  # XXX "peephole" optimization, http://bugs.python.org/issue2506
                continue

            elif o is zoo:
                if not coeff:
                    # 0 * zoo = NaN
                    return [nan], [], None
                coeff = zoo
                continue

            elif o is I:
                neg1e += S.Half
                continue

            elif o.is_commutative:
                #      e
                # o = b
                b, e = o.as_base_exp()

                if o.has(Order):
                    b, e = o, S.One

                #  y
                # 3
                elif o.is_Pow:
                    if b.is_Number:

                        # get all the factors with numeric base so they can be
                        # combined below, but don't combine negatives unless
                        # the exponent is an integer
                        if e.is_Rational:
                            if e.is_Integer:
                                coeff *= Pow(b, e)  # it is an unevaluated power
                                continue
                            elif e.is_negative:    # also a sign of an unevaluated power
                                seq.append(Pow(b, e))
                                continue
                            elif b.is_negative:
                                neg1e += e
                                b = -b
                            if b is not S.One:
                                pnum_rat[b].append(e)
                            o  # XXX "peephole" optimization, http://bugs.python.org/issue2506
                            continue
                        elif b.is_positive or e.is_integer:
                            num_exp.append((b, e))
                            continue

                    elif b is I and e.is_Rational:
                        neg1e += e/2
                        continue

                c_powers.append((b, e))

            # NON-COMMUTATIVE
            # TODO: Make non-commutative exponents not combine automatically
            else:
                if o is not NC_Marker:
                    nc_seq.append(o)

                # process nc_seq (if any)
                while nc_seq:
                    o = nc_seq.pop(0)
                    if not nc_part:
                        nc_part.append(o)
                        continue

                    #                             b    c       b+c
                    # try to combine last terms: a  * a   ->  a
                    o1 = nc_part.pop()
                    b1, e1 = o1.as_base_exp()
                    b2, e2 = o.as_base_exp()
                    new_exp = e1 + e2
                    # Only allow powers to combine if the new exponent is
                    # not an Add. This allow things like a**2*b**3 == a**5
                    # if a.is_commutative == False, but prohibits
                    # a**x*a**y and x**a*x**b from combining (x,y commute).
                    if b1 == b2 and (not new_exp.is_Add) and not (o.has(Order) or o1.has(Order)):
                        o12 = b1 ** new_exp

                        # now o12 could be a commutative object
                        if o12.is_commutative:
                            seq.append(o12)
                            continue
                        else:
                            nc_seq.insert(0, o12)

                    else:
                        nc_part.append(o1)
                        nc_part.append(o)

        # We do want a combined exponent if it would not be an Add, such as
        #  y    2y     3y
        # x  * x   -> x
        # We determine if two exponents have the same term by using
        # as_coeff_Mul.
        #
        # Unfortunately, this isn't smart enough to consider combining into
        # exponents that might already be adds, so things like:
        #  z - y    y
        # x      * x  will be left alone.  This is because checking every possible
        # combination can slow things down.

        # gather exponents of common bases...
        def _gather(c_powers):
            new_c_powers = []
            common_b = {}  # b:e
            for b, e in c_powers:
                co = e.as_coeff_Mul()
                common_b.setdefault(b, {}).setdefault(co[1], []).append(co[0])
            for b, d in common_b.items():
                for di, li in d.items():
                    d[di] = Add(*li)
            for b, e in common_b.items():
                for t, c in e.items():
                    new_c_powers.append((b, c*t))
            return new_c_powers

        # in c_powers
        c_powers = _gather(c_powers)

        # and in num_exp
        num_exp = _gather(num_exp)

        # --- PART 2 ---
        #
        # o process collected powers  (x**0 -> 1; x**1 -> x; otherwise Pow)
        # o combine collected powers  (2**x * 3**x -> 6**x)
        #   with numeric base

        # ................................
        # now we have:
        # - coeff:
        # - c_powers:    (b, e)
        # - num_exp:     (2, e)
        # - pnum_rat:    {(1/3, [1/3, 2/3, 1/4])}

        #  0             1
        # x  -> 1       x  -> x
        for b, e in c_powers:
            if e is S.One:
                assert not b.is_Number
                c_part.append(b)
            elif e is not S.Zero:
                c_part.append(Pow(b, e))

        #  x    x     x
        # 2  * 3  -> 6
        # exp:Mul(num-bases)     x    x
        # e.g.  x:6  for  ... * 2  * 3  * ...
        inv_exp_dict = defaultdict(list)

        for b, e in num_exp:
            inv_exp_dict[e].append(b)
        for e, b in inv_exp_dict.items():
            inv_exp_dict[e] = cls(*b)
        c_part.extend([Pow(b, e) for e, b in inv_exp_dict.items() if e])

        # b, e -> e' = sum(e), b
        # {(1/5, [1/3]), (1/2, [1/12, 1/4]} -> {(1/3, [1/5, 1/2])}
        comb_e = defaultdict(list)
        for b, e in pnum_rat.items():
            comb_e[Add(*e)].append(b)
        del pnum_rat
        # process them, reducing exponents to values less than 1
        # and updating coeff if necessary else adding them to
        # num_rat for further processing
        num_rat = []
        for e, b in comb_e.items():
            b = cls(*b)
            if e.denominator == 1:
                coeff *= Pow(b, e)
                continue
            if e.numerator > e.denominator:
                e_i, ep = divmod(e.numerator, e.denominator)
                coeff *= Pow(b, e_i)
                e = Rational(ep, e.denominator)
            num_rat.append((b, e))
        del comb_e

        # extract gcd of bases in num_rat
        # 2**(1/3)*6**(1/4) -> 2**(1/3+1/4)*3**(1/4)
        pnew = defaultdict(list)
        i = 0  # steps through num_rat which may grow
        while i < len(num_rat):
            bi, ei = num_rat[i]
            grow = []
            for j in range(i + 1, len(num_rat)):
                bj, ej = num_rat[j]
                g = bi.gcd(bj)
                if g is not S.One:
                    # 4**r1*6**r2 -> 2**(r1+r2)  *  2**r1 *  3**r2
                    # this might have a gcd with something else
                    e = ei + ej
                    if e.denominator == 1:
                        coeff *= Pow(g, e)
                    else:
                        if e.numerator > e.denominator:
                            e_i, ep = divmod(e.numerator, e.denominator)  # change e in place
                            coeff *= Pow(g, e_i)
                            e = Rational(ep, e.denominator)
                        grow.append((g, e))
                    # update the jth item
                    num_rat[j] = (bj/g, ej)
                    # update bi that we are checking with
                    bi = bi/g
                    if bi is S.One:
                        break
            if bi is not S.One:
                obj = Pow(bi, ei)
                if obj.is_Number:
                    coeff *= obj
                else:
                    # changes like sqrt(12) -> 2*sqrt(3)
                    for obj in Mul.make_args(obj):
                        if obj.is_Number:
                            coeff *= obj
                        else:
                            assert obj.is_Pow
                            bi, ei = obj.base, obj.exp
                            pnew[ei].append(bi)

            num_rat.extend(grow)
            i += 1

        # combine bases of the new powers
        for e, b in pnew.items():
            pnew[e] = cls(*b)

        # handle -1 and I
        if neg1e:
            # treat I as (-1)**(1/2) and compute -1's total exponent
            p, q = neg1e.as_numer_denom()
            # if the integer part is odd, extract -1
            n, p = divmod(p, q)
            if n % 2:
                coeff = -coeff
            # if it's a multiple of 1/2 extract I
            if q == 2:
                c_part.append(I)
            elif p:
                # see if there is any positive base this power of
                # -1 can join
                neg1e = Rational(p, q)
                for e, b in pnew.items():
                    if e == neg1e and b.is_positive:
                        pnew[e] = -b
                        break
                else:
                    # keep it separate; we've already evaluated it as
                    # much as possible so evaluate=False
                    c_part.append(Pow(S.NegativeOne, neg1e, evaluate=False))

        # add all the pnew powers
        c_part.extend([Pow(b, e) for e, b in pnew.items()])

        # oo, -oo
        if coeff in (oo, -oo):
            def _handle_for_oo(c_part, coeff_sign):
                new_c_part = []
                for t in c_part:
                    if t.is_positive:
                        continue
                    if t.is_negative:
                        coeff_sign *= -1
                        continue
                    new_c_part.append(t)
                return new_c_part, coeff_sign
            c_part, coeff_sign = _handle_for_oo(c_part, 1)
            nc_part, coeff_sign = _handle_for_oo(nc_part, coeff_sign)
            coeff *= coeff_sign

        # zoo
        if coeff is zoo:
            # zoo might be
            #   infinite_real + bounded_im
            #   bounded_real + infinite_im
            #   infinite_real + infinite_im
            # and non-zero real or imaginary will not change that status.
            c_part = [c for c in c_part if not (c.is_nonzero and
                                                c.is_extended_real is not None)]
            nc_part = [c for c in nc_part if not (c.is_nonzero and
                                                  c.is_extended_real is not None)]

        # 0
        elif coeff is S.Zero:
            # we know for sure the result will be 0 except the multiplicand
            # is infinity
            if any(c.is_finite is False for c in c_part):
                return [nan], [], order_symbols
            return [coeff], [], order_symbols

        # check for straggling Numbers that were produced
        _new = []
        for i in c_part:
            if i.is_Number:
                coeff *= i
            else:
                _new.append(i)
        c_part = _new

        # order commutative part canonically
        c_part.sort(key=default_sort_key)

        # current code expects coeff to be always in slot-0
        if coeff is not S.One:
            c_part.insert(0, coeff)

        # we are done
        if (not nc_part and len(c_part) == 2 and c_part[0].is_Number and
                c_part[1].is_Add):
            # 2*(1+a) -> 2 + 2 * a
            coeff = c_part[0]
            c_part = [Add(*[coeff*f for f in c_part[1].args])]

        return c_part, nc_part, order_symbols

    def _eval_power(self, e):

        # don't break up NC terms: (A*B)**3 != A**3*B**3, it is A*B*A*B*A*B
        cargs, nc = self.args_cnc(split_1=False)

        if e.is_Integer:
            return Mul(*[Pow(b, e, evaluate=False) for b in cargs]) * \
                Pow(Mul._from_args(nc), e, evaluate=False)

        p = Pow(self, e, evaluate=False)

        if e.is_Rational or e.is_Float:
            return p._eval_expand_power_base()

        return p

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 4, 0, cls.__name__

    def _eval_evalf(self, prec):
        c, m = self.as_coeff_Mul()
        if c is S.NegativeOne:
            if m.is_Mul:
                rv = -AssocOp._eval_evalf(m, prec)
            else:
                mnew = m._eval_evalf(prec)
                if mnew is not None:
                    m = mnew
                rv = -m
        else:
            rv = AssocOp._eval_evalf(self, prec)
        return rv

    @cacheit
    def as_two_terms(self):
        """Return head and tail of self.

        This is the most efficient way to get the head and tail of an
        expression.

        - if you want only the head, use self.args[0];
        - if you want to process the arguments of the tail then use
          self.as_coef_mul() which gives the head and a tuple containing
          the arguments of the tail when treated as a Mul.
        - if you want the coefficient when self is treated as an Add
          then use self.as_coeff_add()[0]

        >>> (3*x*y).as_two_terms()
        (3, x*y)

        """
        args = self.args

        if len(args) == 2:
            return args
        else:
            return args[0], self._new_rawargs(*args[1:])

    @cacheit
    def as_coeff_mul(self, *deps, **kwargs):
        """Return the tuple (c, args) where self is written as a Mul.

        See Also
        ========

        diofant.core.expr.Expr.as_coeff_mul

        """
        rational = kwargs.pop('rational', True)
        if deps:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(*deps):
                    l2.append(f)
                else:
                    l1.append(f)
            return self._new_rawargs(*l1), tuple(l2)
        args = self.args
        if args[0].is_Number:
            if not rational or args[0].is_Rational:
                return args[0], args[1:]
            elif args[0].is_negative:
                return S.NegativeOne, (-args[0],) + args[1:]
        return S.One, args

    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product."""
        coeff, args = self.args[0], self.args[1:]

        if coeff.is_Number:
            if not rational or coeff.is_Rational:
                if len(args) == 1:
                    return coeff, args[0]
                else:
                    return coeff, self._new_rawargs(*args)
            elif coeff.is_negative:
                return S.NegativeOne, self._new_rawargs(*((-coeff,) + args))
        return S.One, self

    def as_real_imag(self, deep=True, **hints):
        """Returns real and imaginary parts of self

        See Also
        ========

        diofant.core.expr.Expr.as_real_imag

        """
        from .function import expand_mul
        from ..functions import Abs, im, re
        other = []
        coeffr = []
        coeffi = []
        addterms = S.One
        for a in self.args:
            if a.is_extended_real:
                coeffr.append(a)
            elif a.is_imaginary:
                coeffi.append(a)
            elif a.is_commutative:
                # search for complex conjugate pairs:
                for i, x in enumerate(other):
                    if x == a.conjugate():
                        coeffr.append(Abs(x)**2)
                        del other[i]
                        break
                else:
                    if a.is_Add:
                        addterms *= a
                    else:
                        other.append(a)
            else:
                other.append(a)
        m = self.func(*other)
        if hints.get('ignore') == m:
            return
        if len(coeffi) % 2:
            imco = im(coeffi.pop(0))
            # all other pairs make a real factor; they will be
            # put into reco below
        else:
            imco = S.Zero
        reco = self.func(*(coeffr + coeffi))
        r, i = (reco*re(m), reco*im(m))
        if addterms == 1:
            if m == 1:
                if imco is S.Zero:
                    return reco, S.Zero
                else:
                    return S.Zero, reco*imco
            if imco is S.Zero:
                return r, i
            return -imco*i, imco*r
        addre, addim = expand_mul(addterms, deep=False).as_real_imag()
        if imco is S.Zero:
            return r*addre - i*addim, i*addre + r*addim
        else:
            r, i = -imco*i, imco*r
            return r*addre - i*addim, r*addim + i*addre

    @staticmethod
    def _expandsums(sums):
        """
        Helper function for _eval_expand_mul.

        sums must be a list of instances of Basic.

        """
        L = len(sums)
        if L == 1:
            return sums[0].args
        terms = []
        left = Mul._expandsums(sums[:L//2])
        right = Mul._expandsums(sums[L//2:])

        terms = [Mul(a, b) for a in left for b in right]
        added = Add(*terms)
        return Add.make_args(added)  # it may have collapsed down to one term

    def _eval_expand_mul(self, **hints):
        from ..simplify import fraction

        # Handle things like 1/(x*(x + 1)), which are automatically converted
        # to 1/x*1/(x + 1)
        expr = self
        n, d = fraction(expr)
        if d.is_Mul:
            n, d = [i._eval_expand_mul(**hints) if i.is_Mul else i
                    for i in (n, d)]
            expr = n/d
            if not expr.is_Mul:
                return expr

        plain, sums, rewrite = [], [], False
        for factor in expr.args:
            if factor.is_Add:
                sums.append(factor)
                rewrite = True
            else:
                if factor.is_commutative:
                    plain.append(factor)
                else:
                    sums.append(Basic(factor))  # Wrapper

        if not rewrite:
            return expr
        else:
            plain = self.func(*plain)
            terms = self.func._expandsums(sums)
            args = []
            for term in terms:
                t = self.func(plain, term)
                if t.is_Mul and any(a.is_Add for a in t.args):
                    t = t._eval_expand_mul()
                args.append(t)
            return Add(*args)

    @cacheit
    def _eval_derivative(self, s):
        args = list(self.args)
        terms = []
        for i in range(len(args)):
            d = args[i].diff(s)
            if d:
                terms.append(self.func(*(args[:i] + [d] + args[i + 1:])))
        return Add(*terms)

    def _matches_simple(self, expr, repl_dict):
        # handle (w*3)._matches('x*5') -> {w: x*5/3}
        coeff, terms = self.as_coeff_Mul()
        terms = Mul.make_args(terms)
        if len(terms) == 1:
            newexpr = self.__class__._combine_inverse(expr, coeff)
            return terms[0]._matches(newexpr, repl_dict)

    def _matches(self, expr, repl_dict={}):
        """Helper method for match().

        See Also
        ========

        diofant.core.basic.Basic.matches

        """
        expr = sympify(expr)
        if self.is_commutative and expr.is_commutative:
            return AssocOp._matches_commutative(self, expr, repl_dict)
        elif self.is_commutative is not expr.is_commutative:
            return
        c1, nc1 = self.args_cnc()
        c2, nc2 = expr.args_cnc()
        repl_dict = repl_dict.copy()
        if c1:
            if not c2:
                c2 = [1]
            a = self.func(*c1)
            if isinstance(a, AssocOp):
                repl_dict = a._matches_commutative(self.func(*c2), repl_dict)
            else:
                repl_dict = a._matches(self.func(*c2), repl_dict)
        if repl_dict:
            a = self.func(*nc1)
            if not isinstance(a, self.func):
                repl_dict = a._matches(self.func(*nc2), repl_dict)
            else:
                raise NotImplementedError
        return repl_dict or None

    @staticmethod
    def _combine_inverse(lhs, rhs):
        """
        Returns lhs/rhs, but treats arguments like symbols, so things like
        oo/oo return 1, instead of a nan.

        """
        if lhs == rhs:
            return S.One

        def check(l, r):
            if l.is_Float and r.is_comparable:
                # if both objects are added to 0 they will share the same "normalization"
                # and are more likely to compare the same. Since Add(foo, 0) will not allow
                # the 0 to pass, we use __add__ directly.
                return l.__add__(0) == r.evalf(strict=False).__add__(0)
            return False
        if check(lhs, rhs) or check(rhs, lhs):
            return S.One
        if lhs.is_Mul and rhs.is_Mul:
            a = list(lhs.args)
            b = [1]
            for x in rhs.args:
                if x in a:
                    a.remove(x)
                elif -x in a:
                    a.remove(-x)
                    b.append(-1)
                else:
                    b.append(x)
            return lhs.func(*a)/rhs.func(*b)
        return lhs/rhs

    def as_powers_dict(self):
        """Return self as a dictionary of factors with each factor being
        treated as a power.

        See Also
        ========

        diofant.core.expr.Expr.as_powers_dict

        """
        d = defaultdict(int)
        for term in self.args:
            b, e = term.as_base_exp()
            d[b] += e
        return d

    def _eval_as_numer_denom(self):
        """Expression -> a/b -> a, b.

        See Also
        ========

        diofant.core.expr.Expr.as_numer_denom

        """
        # don't use _from_args to rebuild the numerators and denominators
        # as the order is not guaranteed to be the same once they have
        # been separated from each other
        numers, denoms = list(zip(*[f.as_numer_denom() for f in self.args]))
        return self.func(*numers), self.func(*denoms)

    def as_base_exp(self):
        """Return base and exp of self.

        See Also
        ========

        diofant.core.expr.Expr.as_base_exp

        """
        e1 = None
        bases = []
        nc = 0
        for m in self.args:
            b, e = m.as_base_exp()
            if not b.is_commutative:
                nc += 1
            if e1 is None:
                e1 = e
            elif e != e1 or nc > 1:
                return self, S.One
            bases.append(b)
        return self.func(*bases), e1

    def _eval_is_polynomial(self, syms):
        return all(term._eval_is_polynomial(syms) for term in self.args)

    def _eval_is_rational_function(self, syms):
        return all(term._eval_is_rational_function(syms) for term in self.args)

    def _eval_is_algebraic_expr(self, syms):
        return all(term._eval_is_algebraic_expr(syms) for term in self.args)

    def _eval_is_commutative(self):
        return _fuzzy_group(a.is_commutative for a in self.args)

    def _eval_is_finite(self):
        return _fuzzy_group(a.is_finite for a in self.args)

    def _eval_is_complex(self):
        return _fuzzy_group((a.is_complex for a in self.args), quick_exit=True)

    def _eval_is_infinite(self):
        if any(a.is_infinite for a in self.args):
            if not any(not a.is_nonzero for a in self.args):
                return True

    def _eval_is_rational(self):
        r = _fuzzy_group((a.is_rational for a in self.args), quick_exit=True)
        if r:
            return r
        elif r is False:
            return self.is_zero

    def _eval_is_algebraic(self):
        r = _fuzzy_group((a.is_algebraic for a in self.args), quick_exit=True)
        if r:
            return r
        elif r is False:
            return self.is_zero

    def _eval_is_zero(self):
        if any(a.is_zero for a in self.args):
            if all(a.is_finite for a in self.args):
                return True
        elif all(a.is_nonzero for a in self.args):
            return False

    def _eval_is_integer(self):
        is_rational = self.is_rational

        if is_rational:
            n, d = self.as_numer_denom()
            if d is S.One:
                return True
            elif d == 2:
                return n.is_even
        else:
            return is_rational

    def _eval_is_polar(self):
        if all(arg.is_polar or arg.is_positive for arg in self.args):
            return True

    def _eval_is_extended_real(self):
        real = True
        finite = True
        zero = one_neither = False

        for t in self.args:
            if t.is_finite and not t.is_complex:
                return t.is_complex
            elif t.is_imaginary or t.is_extended_real:
                if t.is_imaginary:
                    real = not real
                elif not t.is_real:
                    if zero is not False:
                        return
                    finite = False

                z = t.is_zero
                if not z and zero is False:
                    zero = z
                elif z:
                    if all(a.is_finite for a in self.args):
                        return True
                    return

                if not finite and t.is_real and z is not False:
                    return
            elif t.is_complex and t.is_real is False:
                if one_neither:
                    return  # complex terms might cancel
                one_neither = True
            else:
                return

        if one_neither:  # self is a+I*b or I*b
            if real:
                return zero  # real*self is like self: neither is real
        elif zero is False:
            return real  # can't be trumped by 0
        elif real:
            return real  # doesn't matter what zero is

    def _eval_is_imaginary(self):
        obj = I*self
        if obj.is_Mul:
            return fuzzy_and([obj._eval_is_extended_real(),
                              obj._eval_is_finite()])
        else:
            return obj.is_real

    def _eval_is_irrational(self):
        for t in self.args:
            a = t.is_irrational
            if a:
                if all(x.is_rational and x.is_nonzero
                       for x in self.args if x != t):
                    return True
                return
            elif a is None:
                return
        return False

    def _eval_is_positive(self):
        """Return True if self is positive, False if not, and None if it
        cannot be determined.

        This algorithm is non-recursive and works by keeping track of the
        sign which changes when a negative or nonpositive is encountered.
        Whether a nonpositive or nonnegative is seen is also tracked since
        the presence of these makes it impossible to return True, but
        possible to return False if the end result is nonpositive. e.g.

            pos * neg * nonpositive -> pos or zero -> None is returned
            pos * neg * nonnegative -> neg or zero -> False is returned

        """
        sign = 1
        saw_NON = False
        for t in self.args:
            if t.is_positive:
                continue
            elif t.is_negative:
                sign = -sign
            elif t.is_zero:
                if self.is_finite:
                    return False
                else:
                    return
            elif t.is_nonpositive:
                sign = -sign
                saw_NON = True
            elif t.is_nonnegative:
                saw_NON = True
            else:
                return
        if sign == 1 and saw_NON is False:
            return True
        if sign < 0:
            return False

    def _eval_is_negative(self):
        obj = -self
        if obj.is_Mul:
            return obj._eval_is_positive()
        else:
            return obj.is_positive

    def _eval_is_odd(self):
        is_integer = self.is_integer

        if is_integer:
            r, acc = True, 1
            for t in self.args:
                if not t.is_integer:
                    return
                elif t.is_even or (acc + t).is_odd:
                    r = False
                elif r is False:
                    pass
                elif r and t.is_odd is None:
                    r = None
                acc = t
            return r
        else:
            return is_integer

    def _eval_subs(self, old, new):
        from . import Integer
        from ..functions.elementary.complexes import sign
        from ..ntheory.factor_ import multiplicity
        from ..simplify.powsimp import powdenest
        from ..simplify.radsimp import fraction

        if not old.is_Mul:
            return

        # try keep replacement literal so -2*x doesn't replace 4*x
        if old.args[0].is_Number and old.args[0] < 0:
            if self.args[0].is_Number:
                if self.args[0] < 0:
                    return self._subs(-old, -new)
                return

        def base_exp(a):
            # if I and -1 are in a Mul, they get both end up with
            # a -1 base (see issue sympy/sympy#6421); all we want here are the
            # true Pow separated into base and exponent
            if a.is_Pow:
                return a.as_base_exp()
            return a, S.One

        def breakup(eq):
            """Break up powers of eq when treated as a Mul::

                b**(Rational*e) -> b**e, Rational

            commutatives come back as a dictionary {b**e: Rational}
            noncommutatives come back as a list [(b**e, Rational)]

            """
            c, nc = defaultdict(int), []
            for a in Mul.make_args(eq):
                a = powdenest(a)
                b, e = base_exp(a)
                if e is not S.One:
                    co, _ = e.as_coeff_mul()
                    b = Pow(b, e/co)
                    e = co
                if a.is_commutative:
                    c[b] += e
                else:
                    nc.append([b, e])
            return c, nc

        def rejoin(b, co):
            """
            Put rational back with exponent; in general this is not ok, but
            since we took it from the exponent for analysis, it's ok to put
            it back.

            """
            b, e = base_exp(b)
            return Pow(b, e*co)

        def ndiv(a, b):
            """If b divides a in an extractive way (like 1/4 divides 1/2
            but not vice versa, and 2/5 does not divide 1/3) then return
            the integer number of times it divides, else return 0.

            """
            if not b.denominator % a.denominator or not a.denominator % b.denominator:
                return int(a/b)
            return 0

        # give Muls in the denominator a chance to be changed (see issue sympy/sympy#5651)
        # rv will be the default return value
        rv = None
        n, d = fraction(self)
        self2 = self
        if d is not S.One:
            self2 = n._subs(old, new)/d._subs(old, new)
            if not self2.is_Mul:
                return self2._subs(old, new)
            if self2 != self:
                rv = self2

        # Now continue with regular substitution.

        # handle the leading coefficient and use it to decide if anything
        # should even be started; we always know where to find the Rational
        # so it's a quick test

        co_self = self2.args[0]
        co_old = old.args[0]
        co_xmul = None
        if co_old.is_Rational and co_self.is_Rational:
            # if coeffs are the same there will be no updating to do
            # below after breakup() step; so skip (and keep co_xmul=None)
            if co_old != co_self:
                co_xmul = co_self.extract_multiplicatively(co_old)
        elif co_old.is_Rational:
            return rv

        # break self and old into factors

        c, nc = breakup(self2)
        old_c, old_nc = breakup(old)

        # update the coefficients if we had an extraction
        # e.g. if co_self were 2*(3/35*x)**2 and co_old = 3/5
        # then co_self in c is replaced by (3/5)**2 and co_residual
        # is 2*(1/7)**2

        if co_xmul and co_xmul.is_Rational and abs(co_old) != 1:
            mult = Integer(multiplicity(abs(co_old), co_self))
            c.pop(co_self)
            if co_old in c:
                c[co_old] += mult
            else:
                c[co_old] = mult
            co_residual = co_self/co_old**mult
        else:
            co_residual = 1

        # do quick tests to see if we can't succeed

        ok = True
        if len(old_nc) > len(nc):
            # more non-commutative terms
            ok = False
        elif len(old_c) > len(c):
            # more commutative terms
            ok = False
        elif {i[0] for i in old_nc} - {i[0] for i in nc}:
            # unmatched non-commutative bases
            ok = False
        elif set(old_c) - set(c):
            # unmatched commutative terms
            ok = False
        elif any(sign(c[b]) != sign(old_c[b]) for b in old_c):
            # differences in sign
            ok = False
        if not ok:
            return rv

        if not old_c:
            cdid = None
        else:
            rat = []
            for (b, old_e) in old_c.items():
                c_e = c[b]
                rat.append(ndiv(c_e, old_e))
                if not rat[-1]:
                    return rv
            cdid = min(rat)

        if not old_nc:
            ncdid = None
            for i in range(len(nc)):
                nc[i] = rejoin(*nc[i])
        else:
            ncdid = 0  # number of nc replacements we did
            take = len(old_nc)  # how much to look at each time
            limit = cdid or oo  # max number that we can take
            failed = []  # failed terms will need subs if other terms pass
            i = 0
            while limit and i + take <= len(nc):
                hit = False

                # the bases must be equivalent in succession, and
                # the powers must be extractively compatible on the
                # first and last factor but equal inbetween.

                rat = []
                for j in range(take):
                    if nc[i + j][0] != old_nc[j][0]:
                        break
                    elif j == 0:
                        rat.append(ndiv(nc[i + j][1], old_nc[j][1]))
                    elif j == take - 1:
                        rat.append(ndiv(nc[i + j][1], old_nc[j][1]))
                    elif nc[i + j][1] != old_nc[j][1]:
                        break
                    else:
                        rat.append(1)
                    j += 1
                else:
                    ndo = min(rat)
                    if ndo:
                        if take == 1:
                            assert cdid
                            ndo = min(cdid, ndo)
                            nc[i] = Pow(new, ndo)*rejoin(nc[i][0],
                                                         nc[i][1] - ndo*old_nc[0][1])
                        else:
                            ndo = 1

                            # the left residual

                            l = rejoin(nc[i][0], nc[i][1] - ndo*old_nc[0][1])

                            # eliminate all middle terms

                            mid = new

                            # the right residual (which may be the same as the middle if take == 2)

                            ir = i + take - 1
                            r = (nc[ir][0], nc[ir][1] - ndo*old_nc[-1][1])
                            if r[1]:
                                if i + take < len(nc):
                                    nc[i:i + take] = [l*mid, r]
                                else:
                                    r = rejoin(*r)
                                    nc[i:i + take] = [l*mid*r]
                            else:

                                # there was nothing left on the right

                                nc[i:i + take] = [l*mid]

                        limit -= ndo
                        ncdid += ndo
                        hit = True
                if not hit:

                    # do the subs on this failing factor

                    failed.append(i)
                i += 1
            else:

                if not ncdid:
                    return rv

                # although we didn't fail, certain nc terms may have
                # failed so we rebuild them after attempting a partial
                # subs on them

                failed.extend(range(i, len(nc)))
                for i in failed:
                    nc[i] = rejoin(*nc[i]).subs({old: new})

        # rebuild the expression

        if cdid is None:
            do = ncdid
        elif ncdid is None:
            do = cdid
        else:
            do = min(ncdid, cdid)

        margs = []
        for b in c:
            if b in old_c:

                # calculate the new exponent

                e = c[b] - old_c[b]*do
                margs.append(rejoin(b, e))
            else:
                margs.append(rejoin(b.subs({old: new}), c[b]))
        if cdid and not ncdid:

            # in case we are replacing commutative with non-commutative,
            # we want the new term to come at the front just like the
            # rest of this routine

            margs = [Pow(new, cdid)] + margs
        return co_residual*self2.func(*margs)*self2.func(*nc)

    def _eval_nseries(self, x, n, logx):
        from ..simplify import powsimp
        terms = [t.nseries(x, n=n, logx=logx) for t in self.args]
        return powsimp(self.func(*terms).expand(), combine='exp', deep=True)

    def _eval_as_leading_term(self, x):
        return self.func(*[t.as_leading_term(x) for t in self.args])

    def _eval_conjugate(self):
        return self.func(*[t.conjugate() for t in self.args])

    def _eval_transpose(self):
        return self.func(*[t.transpose() for t in self.args[::-1]])

    def _eval_adjoint(self):
        return self.func(*[t.adjoint() for t in self.args[::-1]])

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self.

        Examples
        ========

        >>> (-3*sqrt(2)*(2 - 2*sqrt(2))).as_content_primitive()
        (6, -sqrt(2)*(-sqrt(2) + 1))

        See Also
        ========

        diofant.core.expr.Expr.as_content_primitive

        """
        coef = S.One
        args = []
        for i, a in enumerate(self.args):
            c, p = a.as_content_primitive(radical=radical)
            coef *= c
            if p is not S.One:
                args.append(p)
        # don't use self._from_args here to reconstruct args
        # since there may be identical args now that should be combined
        # e.g. (2+2*x)*(3+3*x) should be (6, (1 + x)**2) not (6, (1+x)*(1+x))
        return coef, self.func(*args)

    def as_ordered_factors(self, order=None):
        """Transform an expression into an ordered list of factors.

        Examples
        ========

        >>> (2*x*y*sin(x)*cos(x)).as_ordered_factors()
        [2, x, y, sin(x), cos(x)]

        """
        cpart, ncpart = self.args_cnc()
        cpart.sort(key=lambda expr: expr.sort_key(order=order))
        return cpart + ncpart

    @property
    def _sorted_args(self):
        return tuple(self.as_ordered_factors())


def _keep_coeff(coeff, factors, clear=True, sign=False):
    """Return ``coeff*factors`` unevaluated if necessary.

    If ``clear`` is False, do not keep the coefficient as a factor
    if it can be distributed on a single factor such that one or
    more terms will still have integer coefficients.

    If ``sign`` is True, allow a coefficient of -1 to remain factored out.

    Examples
    ========

    >>> _keep_coeff(S.Half, x + 2)
    (x + 2)/2
    >>> _keep_coeff(S.Half, x + 2, clear=False)
    x/2 + 1
    >>> _keep_coeff(S.Half, (x + 2)*y, clear=False)
    y*(x + 2)/2
    >>> _keep_coeff(Integer(-1), x + y)
    -x - y
    >>> _keep_coeff(Integer(-1), x + y, sign=True)
    -(x + y)

    """
    from . import Integer

    if not coeff.is_Number:
        if factors.is_Number:
            factors, coeff = coeff, factors
        else:
            return coeff*factors
    if coeff is S.One:
        return factors
    elif coeff is S.NegativeOne and not sign:
        return -factors
    elif factors.is_Add:
        if not clear and coeff.is_Rational and coeff.denominator != 1:
            q = Integer(coeff.denominator)
            for i in factors.args:
                c, t = i.as_coeff_Mul()
                r = c/q
                if r == int(r):
                    return coeff*factors
        return Mul._from_args((coeff, factors))
    elif factors.is_Mul:
        margs = list(factors.args)
        if margs[0].is_Number:
            margs[0] *= coeff
            if margs[0] == 1:
                margs.pop(0)
        else:
            margs.insert(0, coeff)
        return Mul._from_args(margs)
    else:
        return coeff*factors


def expand_2arg(e):
    from ..simplify.simplify import bottom_up

    def do(e):
        if e.is_Mul:
            c, r = e.as_coeff_Mul()
            if c.is_Number and r.is_Add:
                return Add(*[c*ri for ri in r.args], evaluate=False)
        return e

    return bottom_up(e, do)


from .numbers import I, Rational, nan, oo, zoo
from .power import Pow
from .add import Add
