import hashing
from basic import Basic
from symbol import Symbol, NCSymbol
from numbers import Rational, Real, Number, ImaginaryUnit
from functions import log, exp

class pole_error(Exception):
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
        
        In some cases, the eval() method will return an object that is not an instance of the
        class Add, so for example if x is a Symbol, (x+x) will create a class Add with arguments
        (x,x) , that will be evaluated via the .eval() method, and this method will return a 
        class Mul with arguments (2,x), that is how x+x --> 2*x is done
        
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> type(1+x)
        <class 'sympy.core.addmul.Add'>
        >>> (1+x)[:]
        (1, x)
    
    See also
    ========
        L{Add.eval}
    """
        
    mathml_tag = "power"

    def __init__(self,a,b):
        Basic.__init__(self)
        self._args = [Basic.sympify(a), Basic.sympify(b)]
        
    def hash(self):
        if self._mhash:
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.add(self.base.hash())
        self._mhash.add(self.exp.hash())
        return self._mhash.value

        
    def __str__(self):
        from addmul import Pair
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
    
    @property
    def mathml(self):
        s = "<apply>" + "<" + self.mathml_tag + "/>"
        for a in self._args:
                s += a.mathml
        s += "</apply>"
        return s
        
    @property
    def base(self):
        return self._args[0]
    
    @property
    def exp(self):
        return self._args[1]
        
    def eval(self):
        from addmul import Mul
        if isinstance(self.exp,Rational) and self.exp.iszero():
            return Rational(1)
        if isinstance(self.exp,Rational) and self.exp.isone():
            return self.base
        if isinstance(self.base,Rational) and self.base.iszero():
            if isinstance(self.exp,Rational):# and self.exp.is_integer:
                if self.exp.iszero():
                    raise pole_error("pow::eval(): 0^0.")
                elif self.exp < 0:
                    raise pole_error("pow::eval(): Division by 0.")
            return Rational(0)
        
        if isinstance(self.base,Rational) and self.base.isone():
            return Rational(1)
        
        if isinstance(self.base,Real) and isinstance(self.exp,Real):
            return self
        
        if isinstance(self.base, Rational) and isinstance(self.exp, Rational):
            if self.exp.is_integer:
                if self.exp > 0: 
                    return Rational(self.base.p ** self.exp.p , self.base.q ** self.exp.p)
                else:
                    return Rational(self.base.q ** (-self.exp.p) , self.base.p ** (-self.exp.p) )
                
            if self.base.is_integer:
                a = int(self.base)
                bq = self.exp.q
                if a>0:
                    x = int(a**(1./bq)+0.5)
                    if x**bq == a:
                        assert isinstance(x,int)
                        return Rational(x)**self.exp.p
        if isinstance(self.base,Pow): 
            return Pow(self.base.base,self.base.exp*self.exp)
        if isinstance(self.base,exp): 
            if self.base.is_number:
                return exp(self.exp*self.base._args)
        if isinstance(self.base,Mul): 
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
        if isinstance(self.base, NCSymbol):
            if isinstance(self.exp, Rational) and self.exp.is_integer:
                    n = int(self.exp)
                    #only try to simplify it for low exponents (for speed
                    #reasons).
                    if n > 1 and n < 10:
                        r = self.base
                        for i in range(n-1):
                            r = r * self.base
                        return r
        return self
        

    def evalf(self):
        if self.base.is_number and self.exp.is_number:
            return Real(float(self.base)**float(self.exp))
            #FIXME: we need a way of raising a decimal to the power of a decimal (it doesen't work if self.exp is not an integer
        else:
            raise ValueError

    @property
    def is_commutative(self):
        return self.base.is_commutative and self.exp.is_commutative
        
    def diff(self,sym):
        f = self.base
        g = self.exp
        return (self*(g*log(f)).diff(sym))
        
    def series(self,sym,n):
        if not self.exp.has(sym):
            if isinstance(self.base,Symbol): return self
            try:
                return Basic.series(self,sym,n)
            except pole_error:
                if isinstance(self.exp,Rational) and self.exp.isminusone():
                    g = self.base.series(sym,n)
                    #write g as g=c0*w^e0*(1+Phi)
                    #1/g is then 1/g=c0*w^(-e0)/(1+Phi)
                    #plus we expand 1/(1+Phi)=1-Phi+Phi**2-Phi**3...
                    c0,e0 = g.leadterm(sym)
                    Phi = (g/(c0*sym**e0)-1).expand()
                    e = 0
                    for i in range(n):
                        e += (-1)**i * Phi**i
                    e *= sym ** (-e0) / c0
                    #print n,Phi,c0,e0,g,self.base
                    return e.expand()
                if not isinstance(self.exp,Rational):
                    e = exp(self.exp * log(self.base))
                    return e.series(sym,n)
                #self.base is kind of:  1/x^2 + 1/x + 1 + x + ...
                e = self.base.series(sym,n)
                ldeg = e.ldegree(sym)
                #print "power:",e,self.exp,ldeg,e.eval()
                s= ((e*sym**(-ldeg)).expand()**self.exp).series(sym,n+
                        int(ldeg.evalf()))
                return (s * sym**(ldeg * self.exp)).expand()
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            try:
                a=self.base.series(sym,n)
                b=self.exp.series(sym,n)
                return Basic.series((a**b),sym,n)
            except pole_error:
                e = exp((self.exp*log(self.base)))
                return e.series(sym,n)
            
    def expand(self):
        from addmul import Mul
        if isinstance(self.exp,Number):
            if self.exp.is_integer:
                n = int(self.exp)
                if n > 1:
                    a = self.base
                    while n > 1:
                        a = Mul(a,self.base,evaluate=False)
                        #a *= self.base
                        n -= 1
                    return a.expand()
        return self

    def evalc(self):
        e=self.expand()
        if e!=self:
            return e.evalc()
        if isinstance(e.base, Symbol):
            #this is wrong for nonreal exponent
            return self
        print self
        raise NotImplementedError
        
    def subs(self,old,new):
        if self == old:
            return new
        elif exp(self.exp * log(self.base)) == old:
            return new
        else:
            return (self.base.subs(old,new) ** self.exp.subs(old,new))

    def match(self, pattern, syms):
        def addmatches(r1,r2):
            #print r1,r2
            l1 = list(r1)
            l2 = list(r2)
            if l1 == l2:
                p = l1[0]
                if r1[p] != r2[p]:
                    return None
            r1.update(r2)
            return r1
        assert isinstance(pattern, Pow)
        r1 = self[0].match(pattern[0],syms)
        if r1!=None:
            r2 = self[1].match(pattern[1],syms)
            #print r1,r2,self[1],pattern[1],syms
            if r2!=None:
                return addmatches(r1,r2)
        return None
