"""This module provides an abstract class Function, as well as some mathematical
functions that use Function as its base class. 
"""

import hashing
from basic import Basic
from numbers import Rational, Real
import decimal

class Function(Basic):
    """Abstract class representing a mathematical function. 
    It is the base class for common fuctions such as exp, log, sin, tan, etc.
    """
    
    def __init__(self, arg):
        Basic.__init__(self, is_commutative=True)
        self._args = self.sympify(arg)

    def getname(self):
        return self.__class__.__name__
        
    def __getitem__(self, iter):
        return (self._args,)[iter]
        # do this to force extra nesting and so [:] is coherent across sympy
    
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addint(self._args.hash())
        return self._mhash.value
    
    def diff(self, sym):
        return (self.derivative()*self._args.diff(sym))
    
    def derivative(self):
        return Derivative(self,self._args)
    
    def subs(self, old, new):
        e = Basic.subs(self,old,new)
        #if e==self:
        if e.isequal(self):
            return (type(self)(self._args.subs(old,new)))
        else:
            return e
        
    def __str__(self):
        f = "%s(%s)"
        return f % (self.getname(),str(self._args))

    @property
    def mathml(self):
        return "<apply><%s/> %s </apply>" % (self.mathml_tag, self._args.mathml)
        
    def series(self, sym, n):
        from power import pole_error
        from symbol import Symbol
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            pass
        #this only works, if arg(0) -> 0, otherwise we are in trouble
        arg = self._args.series(sym,n)
        l = Symbol("l",dummy=True)
        #the arg(0) goes to z0
        z0 = arg.subs(log(sym),l).subs(sym,0)
        w = Symbol("w",True)
        e = type(self)(w)
        if arg.has(sym):
            e = e.series(w,n)
        e = e.subs(w,arg-z0)

        #this only works for exp 
        #generally, the problem is with expanding around other point
        #than arg == 0.
        assert isinstance(self,exp)
        e= (exp(z0)*e).expand().subs(l,log(sym))
        return e.expand()
    
    def evalf(self, precision=28):
        """
        Evaluate the current function to a real number.
        
        @param precision: the precision used in the calculations, 
        @type precision: C{int}
        @return: Real number
        
        """
        raise NotImplementedError

class exp(Function):
    """Return e raised to the power of x
    """ 
    
    def derivative(self):
        return exp(self._args)
        
    def expand(self):
        return exp(self._args.expand())
        
    def eval(self):
        arg = self._args
        if isinstance(arg,Rational) and arg.iszero():
            return Rational(1)
        if isinstance(arg,log):
            return arg._args
        return self

    def evalc(self):
        from numbers import I
        from addmul import Mul
        #we will need to move sin,cos to core
        from sympy.modules import cos,sin
        x,y = self._args.get_re_im()
        return exp(x)*cos(y)+I*exp(x)*sin(y)
    
    def evalf(self, precision=28):
        if not self._args.is_number:
            raise ValueError 
        x = Real(self._args) # argument to decimal (full precision)
        decimal.getcontext().prec = precision + 2
        i, lasts, s, fact, num = 0, 0, 1, 1, 1
        while s != lasts:
            lasts = s    
            i += 1
            fact *= i
            num *= x     
            s += num / fact   
        decimal.getcontext().prec = precision - 2        
        return +s

class log(Function):
    """Return the natural logarithm (base e) of x
    """
    
    def derivative(self):
        return Rational(1)/self._args
        
    def eval(self):
        from addmul import Mul
        from power import Pow
        arg = self._args
        if isinstance(arg,Rational) and arg.isone():
            return Rational(0)
        elif isinstance(arg,exp):
            return arg._args
        elif isinstance(arg,Mul):
            a,b = arg.getab()
            return log(a)+log(b)
        elif isinstance(arg,Pow):
            return arg.exp * log(arg.base)
        return self
        
    def evalf(self):
        #TODO: add precision
        import math
        return Real(math.log(self._args.evalf()) )
        
    def series(self,sym,n):
        from numbers import Rational
        from power import pole_error
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            pass
        arg=self._args.series(sym,n)
        #write arg as=c0*w^e0*(1+Phi)
        #log(arg)=log(c0)+e0*log(w)+log(1+Phi)
        #plus we expand log(1+Phi)=Phi-Phi**2/2+Phi**3/3...
        w = sym
        c0,e0 = arg.leadterm(w)
        Phi=(arg/(c0*w**e0)-1).expand()
        if c0.is_number:
            assert c0.evalf()>0
        e=log(c0)+e0*log(w)
        for i in range(1,n+1):
            e+=(-1)**(i+1) * Phi**i /i
        return e

ln = log
    
class abs_(Function):
    """Return the absolute value of x"""
    
    mathml_tag = "abs"
   
    def getname(self):
        return "abs"
 
    def eval(self):
        from addmul import Mul,Add
        from symbol import Symbol
        from numbers import I
        
        arg = self._args
        if arg.is_number or (isinstance(arg, Symbol) and arg.is_real):
            return (arg*arg.conjugate()).expand()**Rational(1,2)
        elif isinstance(arg, Mul):
            _t = arg.getab()[0]
            if _t.is_number and _t < 0:
                return abs(-self._args)
        elif isinstance(arg, Add):
            b,a = arg.getab()
            if isinstance(a, Symbol) and a.is_real:
                if isinstance(b, Mul):
                    a,b=b.getab()
                    if a == I:
                        if isinstance(b, Symbol) and b.is_real:
                            return (arg*arg.conjugate()).expand()**Rational(1,2)
        return self
        
    def evalf(self):
        if self._args.is_number:
            return self.eval()
        else:
            raise ValueError
        
    def derivative(self):
        return sign(self._args)
    
    def series(self):
        pass
    
    def x__eq__(self, a):
        #FIXME: currently this does not work
        # here we are checking for function equality, like in
        # abs(x) == abs(-x)
        if isinstance(a, abs_): 
            if a._args**2 == self._args**2:
                return true
            else:
                return False
        raise ArgumentError("Wrong function arguments")
    
class sign(Function):
    
    def eval(self):
        if self._args.is_number:
            if self._args < 0:
                return Rational(-1)
            elif self._args == 0:
                return Rational(0)
            else:
                return Rational(1)
        return self
            
    def evalf(self, precision=28):
        if self._args.is_number:
            return self.eval()
        else:
            raise ArgumentError
        
    def derivative(self):
        return Rational(0)

class Derivative(Basic):

    def __init__(self,f,x):
        Basic.__init__(self)
        self.f=self.sympify(f)
        self.x=self.sympify(x)

    def eval(self):
        if isinstance(self.f, Derivative):
            if self.f.x != self.x and not self.f.has(self.x):
                return Rational(0)
        return self

    def doit(self):
        return self.f.diff(self.x)

    def diff(self,x):
        return Derivative(self,x)

    def __str__(self):
        if isinstance(self.f,Function):
            return "%s'(%r)"%(self.f.getname(),self.f._args)
        else:
            return "(%r)'"%self.f

    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        self._mhash.addint(self.f.hash())
        self._mhash.addint(self.x.hash())
        return self._mhash.value

