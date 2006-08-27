import hashing
from basic import basic
from symbol import symbol
from numbers import rational,real,number
from functions import ln

class pole_error(Exception):
    pass

class pow(basic):
    def __init__(self,a,b):
        basic.__init__(self)
        self.a=a
        self.b=b
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.add(self.a.hash())
        self.mhash.add(self.b.hash())
        return self.mhash.value
    def printnormal(self):
        from add import pair
        f=""
        if isinstance(self.a,pair) or isinstance(self.a,pow):
            f+="(%s)"
        else:
            f+="%s"
        f+="^"
        if isinstance(self.b,pair) or isinstance(self.b,pow):
            f+="(%s)"
        else:
            f+="%s"
        return f%(self.a,self.b)
    def __str__(self):
        return self.printnormal()
    def getbaseandexp(self):
        return (self.a,self.b)
    def eval(self):
        from add import mul
        if self.evaluated: return self
        self.a=self.a.eval()
        self.b=self.b.eval()
        if isinstance(self.b,rational) and self.b.iszero():
            return rational(1)
        if isinstance(self.b,rational) and self.b.isone():
            return self.a
        if isinstance(self.a,rational) and self.a.iszero():
            if isinstance(self.b,rational) and self.b.isinteger():
                if self.b.iszero():
                    raise pole_error("pow::eval(): 0^0.")
                elif self.b.getinteger()<0:
                    raise pole_error("pow::eval(): Division by 0.")
            return rational(0)
        if isinstance(self.a,rational) and self.a.isone():
            return rational(1)
        if isinstance(self.a,real) and isinstance(self.b,real):
            return self.a.pownumber(self.b)
        if isinstance(self.a,rational) and isinstance(self.b,rational) and \
                self.b.isinteger():
            return self.a.pownumber(self.b)
        if isinstance(self.a,pow): 
            return pow(self.a.a,self.a.b*self.b).eval()
        if isinstance(self.a,mul): 
            a,b=self.a.getab()
            return (pow(a,self.b)*pow(b,self.b)).eval()
        return pow(self.a,self.b).hold()
    def diff(self,sym):
        f=self.a
        g=self.b
        return (self*(g*ln(f)).diff(sym)).eval()
    def series(self,sym,n):
        from add import add
        if isinstance(self.b,rational):# and self.b.isinteger():
            if isinstance(self.a,symbol): return self
            try:
                return basic.series(self,sym,n)
            except pole_error:
                #self.a is kind of:  1/x^2 + 1/x + 1 + x + ...
                e=self.a.series(sym,n).eval()
                ldeg=e.ldegree(sym)
                #print "power:",e,self.b,ldeg,e.eval()
                s= ((e*sym**(-ldeg)).expand()**self.b).series(sym,n+ldeg)
                return (s*sym**(ldeg*self.b)).expand()
        return basic.series(self,sym,n)
    def expand(self):
        if isinstance(self.b,number):
            if self.b.isinteger():
                n=self.b.getinteger()
                if n>1:
                    a=self.a
                    while n>1:
                        a*=self.a
                        n-=1
                    return a.expand()
        return self
    def subs(self,old,new):
        return (self.a.subs(old,new)**self.b.subs(old,new)).eval()
