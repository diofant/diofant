from utils import sign

class basic(object):
    def __init__(self):
        self.evaluated=False;
        self.mhash=0
    def __repr__(self):
        return str(self)
    def __neg__(self):
        from numbers import rational
        return domul(rational(-1),self)
    def __pos__(self):
        return self
    def __add__(self,a):
        return doadd(self,a)
    def __radd__(self,a):
        return doadd(a,self)
    def __sub__(self,a):
        return doadd(self,-a)
    def __rsub__(self,a):
        return doadd(a,-self)
    def __mul__(self,a):
        return domul(self,a)
    def __rmul__(self,a):
        return domul(a,self)
    def __div__(self,a):
        from numbers import rational
        return domul(self,dopow(a,rational(-1)))
    def __rdiv__(self,a):
        from numbers import rational
        return domul(a,dopow(self,rational(-1)))
    def __pow__(self,a):
        return dopow(self,a)
    def __rpow__(self,a):
        return dopow(a,self)
    def __eq__(self,a):
        return self.eval().isequal(c(a).eval())
    def __ne__(self,a):
        return not self.__eq__(a)
    def __lt__(self,a):
        raise "'<' not supported."

    def eval(self):
        """Returns canonical form of myself. 
        
        If we are evaluated (i.e. in the canonical form), the hold 
        method should be called.

        the eval() method should alway return a new object (following the
        general rule of not changing)
        
        """
        return self
    def hold(self):
        """Sets "evaluated" flag. This means, we are in the canonical form,
        and eval don't have to do anything."""
        self.evaluated=True
        return self
    def isequal(self,a):
        return self.hash()==a.hash()
    def cmphash(a,b):
        return sign(a.hash()-b.hash())
    def diffn(self,sym,n):
        while n:
            self=self.diff(sym)
            n-=1
        return self
    def series(self,sym,n):
        from numbers import rational
        f=self
        e=f.subs(sym,rational(0))
        fact=rational(1)
        for i in range(1,n+1):
            fact*=rational(i)
            f=f.diff(sym)
            e+=f.subs(sym,rational(0))*sym**i/fact
        return e.eval()
    def subs(self,old,new):
        if self.isequal(old):
            return new
        else:
            return self
    def has(self,sub):
        from symbol import symbol
        n=symbol("dummy")
        return self.subs(sub,n)!=self
    def leadterm(self,x):
        """Returns the leading term c0*x^e0 of the power series 'self' in x
        with the lowest power of x in a form (c0,e0)
        """
        if not self.evaluated:
            return self.eval().leadterm(x)
        from numbers import rational
        from power import pow
        from add import add,mul
        from symbol import symbol
        def domul(x):
            if len(x)>1:
                return mul(x)
            return x[0]
        def extract(t,x):
            if not t.has(x):
                return t,rational(0)
            if isinstance(t,pow):
                return  rational(1),  t.b
            elif isinstance(t,symbol):
                return  rational(1),  rational(1)
            assert isinstance(t,mul)
            for i,a in enumerate(t.args):
                if a.has(x):
                    if isinstance(a,pow):
                        return  domul(t.args[:i]+t.args[i+1:]),  a.b
                    elif isinstance(a,symbol):
                        return  domul(t.args[:i]+t.args[i+1:]),  rational(1)
                    assert False
            return t,s.rational(0)
        if not isinstance(self,add):
            return extract(self,x)
        lowest=[0,(rational(10)**10).eval()]
        for t in self.args:
            t2=extract(t,x)
            #if t2[1]<lowest[1]:
            if (lowest[1]-t2[1]).evalf()>0:
                lowest=t2
            elif t2[1]==lowest[1]:
                lowest=((lowest[0]+t2[0]).eval(),lowest[1])
        return lowest
    def ldegree(self,sym):
        n =self.leadterm(sym)[1]
        return n.getinteger()
    def expand(self):
        return self

def c(a):
    """for "a" int, returns rational(a), for "a" float returns real, 
    otherwise "a"."""
    from numbers import rational,real
    if isinstance(a,int):
        return rational(a)
    elif isinstance(a,float):
        return real(a)
    else:
        assert isinstance(a,basic)
        return a
def doadd(a,b):
    from add import add
    return add(c(a),c(b))
def domul(a,b):
    from add import mul
    return mul(c(a),c(b))
def dopow(a,b):
    from power import pow
    return pow(c(a),c(b))
