import hashing
from basic import Basic,AutomaticEvaluationType
from numbers import Rational, Real

class Function(Basic):
    """Abstract class representing a mathematical function. 
    It is the base class for common fuctions such as exp, log, sin, tan, etc.
    """
    
#    __metaclass__ = AutomaticEvaluationType

    def __init__(self,arg,evaluated=False):
        Basic.__init__(self,evaluated)
        self.arg = self.sympify(arg)
        
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
        e = Basic.subs(self,old,new)
        #if e==self:
        if e.isequal(self):
            return (type(self)(self.arg.subs(old,new))).eval()
        else:
            return e
        
    def __str__(self):
        f = "%s(%s)"
        return f%(self.getname(),str(self.arg))
    
    def series(self,sym,n):
        from numbers import Rational
        from power import pole_error
        from symbol import Symbol
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            pass
        arg = self.arg.series(sym,n)
        w = Symbol("w",True)
        e = type(self)(w)
        if arg.has(sym):
            e = e.series(w,n)
        e = e.subs(w,arg)
        return e.eval().expand()

class exp(Function):
    """Return e raised to the power of x
    """ 
    
    def getname(self):
        return "exp"
        
    def derivative(self):
        return exp(self.arg)
        
    def expand(self):
        return exp(self.arg.expand()).eval()
        
    def eval(self):
        if self.evaluated: return self
        arg = self.arg.eval()
        if isinstance(arg,Rational) and arg.iszero():
            return Rational(1)
        if isinstance(arg,log):
            return arg.arg
        return exp(arg,evaluated=True)

class log(Function):
    """Return the natural logarithm (base e) of x
    """
    
    def getname(self):
        return "log"
        
    def derivative(self):
        return Rational(1)/self.arg
        
    def eval(self):
        from addmul import Mul
        from power import Pow
        if self.evaluated: return self
        arg=self.arg.eval()
        if isinstance(arg,Rational) and arg.isone():
            return Rational(0)
        elif isinstance(arg,exp):
            return arg.arg.hold()
        elif isinstance(arg,Mul):
            a,b = arg.getab()
            return (log(a)+log(b)).eval()
        elif isinstance(arg,Pow):
            return (arg.exp * log(arg.base)).eval()
        return log(arg,evaluated=True)
        
    def evalf(self):
        import math
        #print type(self.arg)
        return math.log(self.arg.evalf())
        
    def series(self,sym,n):
        from numbers import Rational
        from power import pole_error
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            pass
        arg=self.arg.series(sym,n)
        #write arg as=c0*w^e0*(1+Phi)
        #log(arg)=log(c0)+e0*log(w)+log(1+Phi)
        #plus we expand log(1+Phi)=Phi-Phi**2/2+Phi**3/3...
        w = sym
        c0,e0 = arg.leadterm(w)
        Phi=(arg/(c0*w**e0)-1).expand()
        #print "  LN:",self,c0,w,e0,Phi
        e=log(c0)+e0*log(w)
        #FIXME a huge hack. needs fixing.....
        from sym import limits
        e=e.subs(log(w),limits.whattosubs)
        #print "    LN2:",e
        for i in range(1,n+1):
            e+=(-1)**(i+1) * Phi**i /i
        #print "    LN3:",e.eval()
        return e.eval()

ln = log
