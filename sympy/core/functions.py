"""This module provides an abstract class Function, as well as some mathematical
functions that use Function as its base class. 
"""

from sympy.core.basic import Basic
from sympy.core.symbol import Symbol
from sympy.core.numbers import Rational, Real
import decimal
import math
from sympy.core.stringPict import stringPict, prettyForm

class Function(Basic):
    """Abstract class representing a mathematical function. 
    It is the base class for common fuctions such as exp, log, sin, tan, etc.
    """
    
    def __init__(self, arg):
        Basic.__init__(self, is_commutative=True)
        self._args = self.sympify(arg)

    def __float__(self):
        if self._args.is_number:
            try:
                return eval("math.%s( self._args )" % self.__class__.__name__ )
            except NameError:
                return float(self.evalf())
        else:
            raise ValueError("Cannot evaluate at a symbolic value")

    def getname(self):
        return self.__class__.__name__
        
    def __getitem__(self, iter):
        return (self._args,)[iter]
        # do this to force extra nesting and so [:] is coherent across sympy
    
    def diff(self, sym):
        return (self.derivative()*self._args.diff(sym))
    
    def derivative(self):
        return Derivative(self,self._args)
    
    def subs(self, old, new):
        e = Basic.subs(self,old,new)
        if e == self:
            return (type(self)(self._args.subs(old,new)))
        else:
            return e
        
    def __str__(self):
        f = "%s(%s)"
        return f % (self.getname(),str(self._args))

    def __mathml__(self):
        import xml.dom.minidom
        if self._mathml:
            return self._mathml
        dom = xml.dom.minidom.Document()
        x = dom.createElement("apply")
        x.appendChild(dom.createElement(self.mathml_tag))
        x.appendChild( self._args.__mathml__() )
        self._mathml = x
        return self._mathml

    def __latex__(self):
        s = Symbol(self.getname()).__latex__()
        s += r"\left(%s\right)" % self._args.__latex__()
        return s
        
    def series(self, sym, n=6):
        from power import pole_error
        from symbol import Symbol, O
        from addmul import Add
        if not self.has(O(sym)):
            try:
                return Basic.series(self,sym,n)
            except pole_error:
                pass
        #this only works, if arg(0) -> 0, otherwise we are in trouble
        if not self.has(O(sym)):
            arg = self._args.series(sym,n)
        else:
            arg = self._args
        argorig = arg
        if isinstance(arg,Add):
            arg = arg.removeO()
        l = Symbol("l", dummy=True)
        #the arg(0) goes to z0
        z0 = arg.subs(log(sym),l).subs(sym,0)
        w = Symbol("w",dummy=True)
        e = type(self)(w)
        from addmul import Add
        if arg.has(sym):
            e = e.series(w,n)
            e = e.removeO()
        e = e.subs(w,argorig-z0)

        #this only works for exp 
        #generally, the problem is with expanding around other point
        #than arg == 0.
        assert isinstance(self,exp)
        e= (exp(z0)*e).expand().subs(l,log(sym))
        if isinstance(e,Add) and e.has(O(sym)):
            return e
        return e.expand().series(sym,n)
    
    def evalf(self, precision=28):
        """
        Evaluate the current function to a real number.
        
        @param precision: the precision used in the calculations, 
        @type precision: C{int}
        @return: Real number
        
        """
        raise NotImplementedError
    
    def __pretty__(self):
        """
        Function application.
        Some functions are optimized to omit parentheses.
        They must have a single argument.
        """
        return prettyForm.apply(self.getname(), self._args)

class exp(Function):
    """Return e raised to the power of x
    """ 
    
    def __pretty__(self):
        return prettyForm('e', binding=prettyForm.ATOM)**self._args.__pretty__()
    
    def __latex__(self):
        return "{e}^{%s}" % self[0].__latex__()

    def derivative(self):
        return exp(self._args)

    def combine(self):
        return exp(self[0].combine())
        
    def expand(self):
        from addmul import Add
        a = self[0].expand()
        if isinstance(a, Add):
            r = 1
            for x in a:
                r*=exp(x)
            return r
        return exp(a)
        
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
        from sympy import cos,sin
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
        #elif isinstance(arg,Mul):
        #    a,b = arg.getab()
        #    return log(a)+log(b)
        #elif isinstance(arg,Pow):
        #    return arg.exp * log(arg.base)
        return self

    def expand(self):
        from addmul import Mul
        from power import Pow
        arg = self[0]
        if isinstance(arg,Mul):
            a,b = arg.getab()
            return log(a).expand()+log(b).expand()
        elif isinstance(arg,Pow):
            return arg.exp * log(arg.base).expand()
        else:
            return self
        
    def evalf(self):
        #TODO: add precision
        import math
        return Real(math.log(self._args.evalf()) )
        
    def series(self,sym,n):
        from numbers import Rational
        from power import pole_error
        from symbol import O
        if not self.has(O(sym)):
            try:
                return Basic.series(self,sym,n)
            except pole_error:
                pass
        if not self.has(O(sym)):
            arg=self._args.series(sym,n)
        else:
            arg=self._args
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
            if isinstance(arg, Rational):
                return Rational(abs(arg.p), arg.q)
            else:
                return (arg*arg.conjugate()).expand()**Rational(1,2)
        elif isinstance(arg, Mul):
            _t = arg.getab()[0]
            if _t.is_number and _t < 0:
                return abs(-self._args)
        a = Symbol('a')
        b = Symbol('b')
        match = arg.match(a+I*b, [a,b])
        if  (match is not None) and match[a].is_real and match[b].is_real:
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
    
def sqrt(x):
    return x**(Rational(1,2))
    
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
    
    mathml_tag = 'diff'

    def __init__(self,f,x):
        Basic.__init__(self)
        self.f = self.sympify(f)
        self.x = self.sympify(x)
        self._args = (self.f, self.x)
        #i.e. self[:] = (f, x), which means self = f'(x)
        
    def __pretty__(self):
         f, x = [a.__pretty__() for a in (self.f, self.x)]
         a = prettyForm('d')
         a = prettyForm(*a.below(stringPict.LINE, 'd%s' % str(x)))
         a.baseline = a.baseline + 1
         a = prettyForm(binding=prettyForm.FUNC, *stringPict.next(a, f))
         return a
     
    def __mathml__(self):
        if self._mathml:
            return self._mathml
        import xml.dom.minidom
        dom = xml.dom.minidom.Document()
        x = dom.createElement("apply")
        x.appendChild(dom.createElement(self.mathml_tag))
        
        x_1 = dom.createElement('bvar')
        
        x.appendChild(x_1)
        x.appendChild( self.f.__mathml__() )
        x.childNodes[1].appendChild( self.x.__mathml__() )
        self._mathml = x
        return self._mathml

    def __latex__(self):
        from sympy.core.addmul import Add
        s = r"\frac{\partial}{\partial %s} " % self.x.__latex__()
        if isinstance(self.f, Add):
            s += r"\left(" + self.f.__latex__() + r"\right)"
        else:
            s += self.f.__latex__()
        return s

    def eval(self):
        from addmul import Mul
        if isinstance(self.f, Derivative):
            if self.f.x != self.x and not self.f.has(self.x):
                return Rational(0)
        if isinstance(self.f, Mul):
            #(2*x)' -> 2*x'
            atoms = self.f[:]
            with_x = []
            without_x = []
            for x in atoms:
                if x.has(self.x):
                    with_x.append(x)
                else:
                    without_x.append(x)
            if len(without_x) == 0:
                return self
            elif len(without_x) == 1:
                a = without_x[0]
            else:
                a = Mul(*without_x)
            if len(with_x) == 0:
                b = 1
            elif len(with_x) == 1:
                b = with_x[0]
            else:
                b = Mul(*with_x)
            return a*Derivative(b, self.x)
        return self

    def doit(self):
        return self.f.doit().diff(self.x)

    def diff(self,x):
        return Derivative(self,x)

    def __str__(self):
        if isinstance(self.f,Function):
            return "%s'(%r)" % (self.f.getname(),self.f._args)
        else:
            return "(%r)'" % self.f

    def subs(self, old, new):
        e = Basic.subs(self,old,new)
        #if e==self:
        if e == self:
            return Derivative(self[0].subs(old,new), self[1])
        else:
            return e

def diff(f, x, times = 1, evaluate=True):
    """Derivate f with respect to x
    
    It's just a wrapper to unify .diff() and the Derivative class, 
    it's interface is similar to that of integrate()
    
    see http://documents.wolfram.com/v5/Built-inFunctions/AlgebraicComputation/Calculus/D.html
    """
    f = Basic.sympify(f)
    if evaluate == True:
        for i in range(0,times):
            f = f.diff(x)
        return f
    else:
        return Derivative(f, x)
