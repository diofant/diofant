import hashing
from basic import Basic,AutomaticEvaluationType
from symbol import Symbol
from numbers import Rational,Real,Number,ImaginaryUnit
from functions import log,exp

class pole_error(Exception):
    pass

class Pow(Basic):

    #__metaclass__ = AutomaticEvaluationType
    
    def __init__(self,a,b,evaluated=False):
        Basic.__init__(self,evaluated)
        self.base = self.sympify(a)
        self.exp = self.sympify(b)
        
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.add(self.base.hash())
        self.mhash.add(self.exp.hash())
        return self.mhash.value
        
    def print_normal(self):
        from addmul import Pair
        f = ""
        if isinstance(self.base,Pair) or isinstance(self.base,Pow):
            f += "(%s)"
        else:
            f += "%s"
        f += "^"
        if isinstance(self.exp,Pair) or isinstance(self.exp,Pow):
            f += "(%s)"
        else:
            f += "%s"
        return f % (self.base,self.exp)
        
    def __str__(self):
        return self.print_normal()
        
    def get_baseandexp(self):
        return (self.base,self.exp)
        
    def eval(self):
        from addmul import Mul
        if self.evaluated: return self
        self.base = self.base.eval()
        self.exp = self.exp.eval()
        if isinstance(self.exp,Rational) and self.exp.iszero():
            return Rational(1)
        if isinstance(self.exp,Rational) and self.exp.isone():
            return self.base
        if isinstance(self.base,Rational) and self.base.iszero():
            if isinstance(self.exp,Rational):# and self.exp.isinteger():
                if self.exp.iszero():
                    raise pole_error("pow::eval(): 0^0.")
                elif self.exp < 0:
                    raise pole_error("pow::eval(): Division by 0.")
            return Rational(0)
        if isinstance(self.base,Rational) and self.base.isone():
            return Rational(1)
        if isinstance(self.base,Real) and isinstance(self.exp,Real):
            return self.base ** self.exp
        if isinstance(self.base,Rational) and isinstance(self.exp,Rational):
            if self.exp.isinteger(): 
                return self.base ** self.exp
            if self.base.isinteger():
                a = self.base.getinteger()
                bq = self.exp.q
                x = int(a**(1./bq)+0.5)
                if x**bq == a:
                    assert isinstance(x,int)
                    return (Rational(x)**self.exp.p).eval()
        if isinstance(self.base,Pow): 
            return Pow(self.base.base,self.base.exp*self.exp).eval()
        if isinstance(self.base,Mul): 
            a,b = self.base.getab()
            return (Pow(a,self.exp) * Pow(b,self.exp)).eval()
        if isinstance(self.base,ImaginaryUnit):
            if isinstance(self.exp,Rational) and self.exp.isinteger():
                if self.exp.getinteger()==2:
                    return (-Rational(1)).eval()
        return Pow(self.base,self.exp,evaluated=True)
        
    def evalf(self):
        if hasattr(self.base, 'evalf') and hasattr(self.exp, 'evalf'):
            return self.base.evalf()**self.exp.evalf()
        else: 
            raise ValueError('Can not evaluate a symbolic value')
        
    def diff(self,sym):
        f = self.base
        g = self.exp
        return (self*(g*log(f)).diff(sym)).eval()
        
    def series(self,sym,n):
        from addmul import Add
        #if isinstance(self.exp,Rational):
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
                    return e.eval()
                if not isinstance(self.exp,Rational):
                    e = exp(self.exp * log(self.base)).eval()
                    return e.series(sym,n)
                #self.base is kind of:  1/x^2 + 1/x + 1 + x + ...
                e = self.base.series(sym,n).eval()
                ldeg = e.ldegree(sym)
                #print "power:",e,self.exp,ldeg,e.eval()
                s= ((e*sym**(-ldeg)).expand()**self.exp).series(sym,n+
                        int(ldeg.evalf()))
                return (s * sym**(ldeg * self.exp)).expand()
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            e = exp(self.exp*log(self.base)).eval()
            return e.series(sym,n)
            
    def expand(self):
        if isinstance(self.exp,Number):
            if self.exp.isinteger():
                n=self.exp.getinteger()
                if n > 1:
                    a = self.base
                    while n > 1:
                        a *= self.base
                        n -= 1
                    return a.expand()
        return self
        
    def subs(self,old,new):
        if self == old:
            return new
        elif exp(self.exp * log(self.base)) == old:
            return new
        else:
            return (self.base.subs(old,new) ** self.exp.subs(old,new)).eval()
