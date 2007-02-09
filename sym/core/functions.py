import hashing
from basic import basic,c
from numbers import rational

class function(basic):
    def __init__(self,arg):
        basic.__init__(self)
        self.arg=c(arg)
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.arg.hash())
        return self.mhash.value
    def diff(self,sym):
        return (self.derivative()*self.arg.diff(sym)).eval()
    def subs(self,old,new):
        e=basic.subs(self,old,new)
        #if e==self:
        if e.isequal(self):
            return (type(self)(self.arg.subs(old,new))).eval()
        else:
            return e
    def __str__(self):
        f="%s(%s)"
        return f%(self.getname(),str(self.arg))
    def series(self,sym,n):
        from numbers import rational
        from power import pole_error
        from symbol import symbol
        try:
            return basic.series(self,sym,n)
        except pole_error:
            pass
        arg=self.arg.series(sym,n)
        w=symbol("w",True)
        e=type(self)(w)
        if arg.has(sym):
            e=e.series(w,n)
        e=e.subs(w,arg)
        return e.eval().expand()

class exp(function):
    def getname(self):
        return "exp"
    def derivative(self):
        return exp(self.arg)
    def expand(self):
        return exp(self.arg.expand()).eval()
    def eval(self):
        if self.evaluated: return self
        arg=self.arg.eval()
        if isinstance(arg,rational) and arg.iszero():
            return rational(1)
        if isinstance(arg,ln):
            return arg.arg
        return exp(arg).hold()

class ln(function):
    def getname(self):
        return "ln"
    def derivative(self):
        return rational(1)/self.arg
    def eval(self):
        from add import mul
        from power import pow
        if self.evaluated: return self
        arg=self.arg.eval()
        if isinstance(arg,rational) and arg.isone():
            return rational(0)
        elif isinstance(arg,exp):
            return arg.arg.hold()
        elif isinstance(arg,mul):
            a,b=arg.getab()
            return (ln(a)+ln(b)).eval()
        elif isinstance(arg,pow):
            return (arg.b*ln(arg.a)).eval()
        return ln(arg).hold()
    def evalf(self):
        import math
        #print type(self.arg)
        return math.log(self.arg.evalf())
    def series(self,sym,n):
        from numbers import rational
        from power import pole_error
        try:
            return basic.series(self,sym,n)
        except pole_error:
            pass
        arg=self.arg.series(sym,n)
        #write arg as=c0*w^e0*(1+Phi)
        #ln(arg)=ln(c0)+e0*ln(w)+ln(1+Phi)
        #plus we expand ln(1+Phi)=Phi-Phi**2/2+Phi**3/3...
        w=sym
        c0,e0=arg.leadterm(w)
        Phi=(arg/(c0*w**e0)-1).expand()
        #print "  LN:",self,c0,w,e0,Phi
        e=ln(c0)+e0*ln(w)
        import limits
        e=e.subs(ln(w),limits.whattosubs)
        #print "    LN2:",e
        for i in range(1,n+1):
            e+=(-1)**(i+1) * Phi**i /i
        #print "    LN3:",e.eval()
        return e.eval()
