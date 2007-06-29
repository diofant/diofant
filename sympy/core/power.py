from basic import Basic
from symbol import Symbol, O
from numbers import Rational, Real, ImaginaryUnit
from functions import log, exp
from sympy.core.stringPict import prettyForm, stringPict


def integer_nthroot(y, n):
    """
    Usage
    =====
        Return a tuple containing x = floor(y**(1/n))
        and a boolean indicating whether the result is exact (that is,
        whether x**n == y).

    Examples
    ========

        >>> integer_nthroot(16,2)
        (4, True)
        >>> integer_nthroot(26,2)
        (5, False)
        >>> integer_nthroot(1234567**7, 7)
        (1234567L, True)
        >>> integer_nthroot(1234567**7+1, 7)
        (1234567L, False)
    """

    y = int(y); n = int(n)
    if y < 0: raise ValueError, "y must not be negative"
    if n < 1: raise ValueError, "n must be positive"
    if y in (0, 1): return y, True
    if n == 1: return y, True

    # Search with Newton's method, starting from floating-point
    # approximation. Care must be taken to avoid overflow.
    from math import log as _log
    guess = 2**int(_log(y, 2)/n)
    xprev, x = -1, guess
    while abs(x - xprev) > 1:
        t = x**(n-1)
        xprev, x = x, x - (t*x-y)//(n*t)
    # Compensate
    while x**n > y:
        x -= 1
    return x, x**n == y

class pole_error(ZeroDivisionError):
    pass

class Pow(Basic):
    """
    Usage
    =====
        This class represent's the power of two elements. so whenever you call '**', an
        instance of this class is created.

    Notes
    =====
        When an instance of this class is created, the method .eval() is called and will
        preform some inexpensive symplifications.


    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> type(x**2)
        <class 'sympy.core.power.Pow'>
        >>> (x**2)[:]
        [x, 2]

    See also
    ========
        L{Add.eval}
    """

    mathml_tag = "power"

    def __init__(self,a,b):
        Basic.__init__(self)
        self._args = [Basic.sympify(a), Basic.sympify(b)]


    def __str__(self):
        from addmul import Pair
        if self.exp == -1:
            if isinstance(self.base, Pair):
                return "(1/(%s))" % str(self.base)
            else:
                return "(1/%s)" % str(self.base)

        f = ""
        if isinstance(self.base,Pair) or isinstance(self.base,Pow):
            f += "(%s)"
        else:
            f += "%s"
        f += "**"
        if isinstance(self.exp,Pair) or isinstance(self.exp,Pow) \
            or (isinstance(self.exp,Rational) and \
            (not self.exp.is_integer or (self.exp.is_integer and \
            int(self.exp) < 0)) ):
            f += "(%s)"
        else:
            f += "%s"
        return f % (str(self.base), str(self.exp))


    def __pretty__(self):
        if self.exp == Rational(1,2): # if it's a square root
            bpretty = self.base.__pretty__()
            bl = int((bpretty.height() / 2.0) + 0.5)

            s2 = stringPict("\\/")
            for x in xrange(1, bpretty.height()):
                s3 = stringPict(" " * (2*x+1) + "/")
                s2 = stringPict(*s2.top(s3))
            s2.baseline = -1

            s = stringPict("__" + "_" * bpretty.width())
            s = stringPict(*s.below("%s" % str(bpretty)))
            s = stringPict(*s.left(s2))
            return prettyForm(str(s), baseline=bl)
        elif self.exp == -1:
            # things like 1/x
            return prettyForm("1") / self.base.__pretty__()
        a, b = self._args
        return a.__pretty__()**b.__pretty__()


    def __latex__(self):
        from addmul import Pair
        f = ""
        if isinstance(self.base,Pair) or isinstance(self.base,Pow):
            f += "{(%s)}"
        else:
            f += "{%s}"
        f += "^"
        if isinstance(self.exp,Pair) or isinstance(self.exp,Pow) \
            or (isinstance(self.exp,Rational) and \
            (not self.exp.is_integer or (self.exp.is_integer and \
            int(self.exp) < 0)) ):
            f += "{(%s)}"
        else:
            f += "{%s}"
        return f % (self.base.__latex__(),self.exp.__latex__())

    def __mathml__(self):
        import xml.dom.minidom
        if self._mathml:
            return self._mathml
        dom = xml.dom.minidom.Document()
        x = dom.createElement("apply")
        x_1 = dom.createElement(self.mathml_tag)
        x.appendChild(x_1)
        for arg in self._args:
            x.appendChild( arg.__mathml__() )
        self._mathml = x
        return self._mathml

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    def eval(self):
        from addmul import Mul
        from numbers import oo
        if isinstance(self.exp, Rational) and self.exp.iszero():
            return Rational(1)

        if isinstance(self.exp, Rational) and self.exp.isone():
            return self.base

        if isinstance(self.base, Rational) and self.base.iszero():
            if isinstance(self.exp,Rational):# and self.exp.is_integer:
                if self.exp.iszero():
                    raise pole_error("pow::eval(): 0^0.")
                elif self.exp < 0:
                    raise pole_error("%s: Division by 0." % str(self))
            return Rational(0)

        if isinstance(self.base, Rational) and self.base.isone():
            return Rational(1)

        if isinstance(self.base, Real) and isinstance(self.exp,Real):
            return self

        if isinstance(self.base, Rational) and isinstance(self.exp, Rational):

            if self.exp.is_integer:
                if self.exp > 0:
                    return Rational(self.base.p ** self.exp.p , self.base.q ** self.exp.p)
                else:
                    return Rational(self.base.q ** (-self.exp.p) , self.base.p ** (-self.exp.p) )

            if self.base == -1:
                # calculate the roots of -1
                from sympy.modules.trigonometric import sin, cos
                from sympy.core.numbers import pi
                r = cos(pi/self.exp.q) + ImaginaryUnit()*sin(pi/self.exp.q)
                return r**self.exp.p

            # Rational ** Rational, general case. Calculate explicitly whenever possible.
            a = self.base
            b = self.exp
            if a >= 0:
                x, xexact = integer_nthroot(a.p, b.q)
                y, yexact = integer_nthroot(a.q, b.q)
                if xexact and yexact:
                    res = Rational(x ** abs(b.p), y ** abs(b.p))
                    if b >= 0:
                        return res
                    else:
                        return 1/res
            else:
                return Pow(-1, b) * Pow(-a, b)

            # left out:
            # case base negative && not a perfect number, like sqrt(-2)
            # TODO: implement for exponent of 1/4, 1/6, 1/8, etc.
            # return ((-1)**self.exp)*Pow(-self.base, self.exp, evaluate=False)

        if isinstance(self.base, Pow):
            return Pow(self.base.base,self.base.exp*self.exp)

        if isinstance(self.base, exp):
            if self.base.is_number:
                return exp(self.exp*self.base._args)

        if isinstance(self.base, Mul):
            a,b = self.base.getab()
            if self.exp==-1 or (isinstance(a,Rational) and a.evalf()>0):
                return (Pow(a,self.exp) * Pow(b,self.exp))

        if isinstance(self.base,ImaginaryUnit):
            if isinstance(self.exp,Rational) and self.exp.is_integer:
                if int(self.exp) % 2 == 0:
                    return Rational(-1) ** ((int(self.exp) % 4)/2)

        if isinstance(self.exp,Rational) and self.exp.is_integer:
            if isinstance(self.base,Mul):
                if int(self.exp) % 2 == 0:
                    n = self.base[0]
                    if n.is_number and n < 0:
                        return (-self.base)**self.exp

        if isinstance(self[0],Real) and self[1].is_number:
            return Real(self[0]**self[1].evalf())

        if not self.base.is_commutative:
            if isinstance(self.exp, Rational) and self.exp.is_integer:
                    n = int(self.exp)
                    #only try to simplify it for low exponents (for speed
                    #reasons).
                    if n > 1 and n < 10:
                        r = self.base
                        for i in range(n-1):
                            r = r * self.base
                        return r

        if self.exp == -1 and self.base == oo:
            return Rational(0)
        return self


    def evalf(self):
        if self.base.is_number and self.exp.is_number:
            return Real(float(self.base)**float(self.exp))
            #FIXME: we need a way of raising a decimal to the power of a decimal (it doesen't work if self.exp is not an integer
        else:
            raise ValueError

    def diff(self,sym):
        f = self.base
        g = self.exp
        return (self*(g*log(f)).diff(sym))

    def series(self,sym,n):
        from addmul import Add
        if not self.exp.has(sym):
            if isinstance(self.base,Symbol): return self
            try:
                return Basic.series(self,sym,n)
            except pole_error:
                if isinstance(self.exp,Rational) and self.exp.isminusone():
                    g = self.base.series(sym,n)
                    if isinstance(g, Add):
                        g = g.removeO()
                    elif isinstance(g, O):
                        g = Rational(0)
                    #write g as g=c0*w^e0*(1+Phi)
                    #1/g is then 1/g=c0*w^(-e0)/(1+Phi)
                    #plus we expand 1/(1+Phi)=1-Phi+Phi**2-Phi**3...
                    c0,e0 = g.leadterm(sym)
                    Phi = (g/(c0*sym**e0)-1).expand()
                    e = 0
                    #n-=1
                    for i in range(n):
                        e += (-1)**i * Phi**i
                    e+=O(Phi**n)
                    e *= sym ** (-e0) / c0
                    return e.expand()
                if not isinstance(self.exp,Rational):
                    e = exp(self.exp * log(self.base))
                    return e.series(sym,n)
                #self.base is kind of:  1/x^2 + 1/x + 1 + x + ...
                e = self.base.series(sym,n)
                if isinstance(e, Add):
                    e = e.removeO()
                ldeg = e.ldegree(sym)
                s= ((e*sym**(-ldeg)).expand()**self.exp).series(sym,n+
                        int(ldeg.evalf()))
                return (s * sym**(ldeg * self.exp)).expand()
        try:
            if self.has(O(sym)):
                e = exp(((self.exp*log(self.base).series(sym,n)).expand())).series(sym,n)
                return e
            return Basic.series(self,sym,n)
        except pole_error:
            e = exp((self.exp*log(self.base)))
            return e.series(sym,n)
            try:
                a=self.base.series(sym,n)
                b=self.exp.series(sym,n)
                return Basic.series((a**b),sym,n)
            except pole_error:
                e = exp((self.exp*log(self.base)))
                return e.series(sym,n)

    def expand(self):
        from addmul import Add, Mul

        def _expand_bin(a, b, n):
            """calculates the expansion of (a+b)**n using newton's binomial
            formula (also called triangle of Pascal)

            See L{http://en.wikipedia.org/wiki/Binomial_theorem}
            """
            s = a**n
            cur_coeff = Rational(1)
            cur_exp = 1
            for i in range(1, n+1):
                # we could speed this using that the coefficients are
                # symetrical
                cur_coeff = cur_coeff * (n-i+1) / i
                s += cur_coeff * (a**(n-i)) * (b**(i))
            return s

        def _expand_multi(**args):
            """calculate the expansion of (a1 + a2 + ...)**n
            TODO
            """
            pass

        if isinstance(self.exp, (Real, Rational)):
            if self.exp.is_integer:
                n = int(self.exp)
                if n > 1:
                    base = self.base.expand()
                    if isinstance(base, Add) and self.base.is_commutative:
                        # try to expand using the binomial formula
                        if len(base[:]) == 2:
                            a, b = base.getab()
                            return _expand_bin(a, b, n)
                        else:
                            #implement the multinomial formula
                            pass
                    a = self.base
                    while n > 1:
                        a = Mul(a,self.base,evaluate=False)
                        #a *= self.base
                        n -= 1
                    return a.expand()
                if n < 0:
                    p = Pow(self.base, -self.exp).expand()
                    if isinstance(p, Mul):
                        return Mul(*(a**(-1) for a in p[:]))

        return Pow(self[0].expand(),self[1].expand())
        return self

    def combine(self):
        from functions import exp
        a = self[0].combine()
        b = self[1].combine()
        if isinstance(a, exp):
            return exp(a[0]*b)
        return self

    def evalc(self):
        if self.exp.is_number:
            c = self.base.evalc()
            if self.exp.is_integer:
                r = c ** self.exp
                re = r.expand()
                if re == r:
                    return re
                else:
                    return re.evalc()
            else:
                from sympy.modules.trigonometric import atan, cos, sin
                from sympy.core.numbers import I
                re,im = c.get_re_im()
                r = (re**2 + im**2)**Rational(1,2)
                t = atan(im / re)

                #tp = ((t + 2*k*pi) / self.exp)*I  # this is right, for k=...,-2,-1,0,1,2,...
                rp = r**self.exp
                tp = t*self.exp

                return rp*cos(tp) + rp*sin(tp)*I
        else:
            return self

    def subs(self,old,new):
        if self == old:
            return new
        elif exp(self.exp * log(self.base)) == old:
            return new
        else:
            return (self.base.subs(old,new) ** self.exp.subs(old,new))

    def match(self, pattern, syms=None):
        from symbol import Symbol
        if syms == None:
            syms = pattern.atoms(type=Symbol)
        def addmatches(r1,r2):
            l1 = list(r1)
            l2 = list(r2)
            if l1 == l2:
                if l1 != []:
                    #fix it in a general case
                    assert len(l1)==1
                    p = l1[0]
                    if r1[p] != r2[p]:
                        return None
            r1.update(r2)
            return r1
        from addmul import Mul
        if isinstance(pattern, Mul):
            return Mul(Rational(1),self,evaluate = False).match(pattern,syms)
        if not isinstance(pattern, Pow):
            return None
        r1 = self[0].match(pattern[0],syms)
        if r1!=None:
            r2 = self[1].match(pattern[1],syms)
            #print r1,r2,"<--",self[1],pattern[1],syms
            if r2!=None:
                return addmatches(r1,r2)
        return None

    @property
    def is_commutative(self):
        return self.base.is_commutative and self.exp.is_commutative

    @property
    def is_integer(self):
        if self.base.is_integer and self.exp.is_integer:
            return self.exp.is_nonnegative
        else:
            return None

    @property
    def is_even(self):
        return self.base.is_even and self.exp.is_nonnegative_integer

    @property
    def is_odd(self):
        return self.base.is_odd and self.exp.is_nonnegative_integer

#    @property
#    def is_negative(self):
#        pass

#    @property
#    def is_positive(self):
#        pass

#    @property
#    def is_nonpositive(self):
#        pass

#    @property
#    def is_nonnegative(self):
#        pass

