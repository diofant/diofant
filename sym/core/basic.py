from utils import sign

class AutomaticEvaluationType(type):
    def __call__(self,*args,**kwargs):
        obj=type.__call__(self,*args,**kwargs)
        return obj.eval()

class Basic(object):
    
#    __metaclass__ = AutomaticEvaluationType
    
    def __init__(self,evaluated=False):
        self.evaluated = evaluated;
        self.mhash = 0
        
    def __repr__(self):
        return str(self)
    
    def __neg__(self):
        from numbers import Rational
        return self._domul(Rational(-1),self)
        
    def __pos__(self):
        return self
        
    def __add__(self,a):
        return self._doadd(self, a)
    
    def __abs__(self):
        """Returns the absolute value of self. 
        
        Example usage: 
        >>> from sym import *
        >>> abs(1+2*I)
        5^1/2
        >>> x = Symbol('x')
        >>> abs(-x)
        x
        """
        return (self*self.conjugate()).expand().sqrt()
        
    def __radd__(self,a):
        return self._doadd(a, self)
        
    def __sub__(self,a):
        return self._doadd(self, -a)
        
    def __rsub__(self,a):
        return self._doadd(a, -self)
        
    def __mul__(self,a):
        return self._domul(self, a)
        
    def __rmul__(self,a):
        return self._domul(a, self)
        
    def __div__(self,a):
        from numbers import Rational
        return self._domul(self,self._dopow(a,Rational(-1)))
        
    def __rdiv__(self,a):
        from numbers import Rational
        return self._domul(self.sympify(self) ** Rational(-1), self.sympify(a))
        
    def __pow__(self,a):
        return self._dopow(self, a)
        
    def __rpow__(self,a):
        return self._dopow(a, self)
        
    def __eq__(self,a):
        return self.eval().isequal(self.sympify(a).eval())
        
    def __ne__(self,a):
        return not self.__eq__(a)
        
    def __lt__(self,a):
        raise NotImplementedError("'<' not supported.")

    @staticmethod
    def _doadd(a,b):
        from addmul import Add
        return Add(Basic.sympify(a), Basic.sympify(b))

    @staticmethod
    def _domul(a, b):
        from addmul import Mul
        return Mul(Basic.sympify(a), Basic.sympify(b))

    @staticmethod
    def _dopow(a, b):
        from addmul import Pow
        return Pow(Basic.sympify(a), Basic.sympify(b))

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
        self.evaluated = True
        return self
    
    @staticmethod
    def sympify(a):
        """for "a" int, returns Rational(a), for "a" float returns real, 
        otherwise "a" (=it's a Basic subclass)."""
        from numbers import Rational, Real
        if isinstance(a,int):
            return Rational(a)
        elif isinstance(a,float):
            return Real(a)
        else:
            assert isinstance(a,Basic)
            return a
        
    def isequal(self,a):
        return self.hash() == a.hash()
        
    def cmphash(a,b):
        return sign(a.hash()-b.hash())
        
    def diffn(self,sym,n):
        while n:
            self = self.diff(sym)
            n -= 1
        return self
        
    def series(self,sym,n):
        from numbers import Rational
        f = self
        e = f.subs(sym,Rational(0))
        fact = Rational(1)
        for i in range(1,n+1):
            fact *= Rational(i)
            f = f.diff(sym)
            e += f.subs(sym,Rational(0))*(sym**i)/fact
        return e.eval()
        
    def subs(self,old,new):
        if self.isequal(old):
            return self.sympify(new)
        else:
            return self
            
    def has(self,sub):
        from symbol import Symbol
        n = Symbol("dummy")
        return self.subs(sub,n)!=self
        
    def leadterm(self,x):
        """Returns the leading term c0*x^e0 of the power series 'self' in x
        with the lowest power of x in a form (c0,e0)
        """
        
        if not self.evaluated:
            return self.eval().leadterm(x)
        from numbers import Rational
        from power import Pow
        from addmul import Add,Mul
        from symbol import Symbol
        def domul(x):
            if len(x) > 1:
                return Mul(x)
            return x[0]
        def extract(t,x):
            if not t.has(x):
                return t,Rational(0)
            if isinstance(t,Pow):
                return  Rational(1),  t.exp
            elif isinstance(t,Symbol):
                return  Rational(1),  Rational(1)
            assert isinstance(t,Mul)
            for i,a in enumerate(t.args):
                if a.has(x):
                    if isinstance(a,Pow):
                        return  domul(t.args[:i] + t.args[i+1:]),  a.exp
                    elif isinstance(a,Symbol):
                        return  domul(t.args[:i] + t.args[i+1:]),  Rational(1)
                    assert False
            return t,s.Rational(0)
        if not isinstance(self,Add):
            return extract(self,x)
        lowest = [0,(Rational(10)**10).eval()]
        for t in self.args:
            t2 = extract(t,x)
            #if t2[1]<lowest[1]:
            if (lowest[1] - t2[1]).evalf()>0:
                lowest=t2
            elif t2[1] == lowest[1]:
                lowest=((lowest[0] + t2[0]).eval(),lowest[1])
        return lowest
        
    def ldegree(self,sym):
        """Returns the lowest power of the sym
        """
        return self.leadterm(sym)[1]
        
    def expand(self):
        return self

    def conjugate(self):
        """Returns a  complex conjugate of self. 
        
        Note: this implementeation assumes that all Symbols are real,
        so we just need to change the sign at "i".
        """
        from numbers import I
        return self.subs(I,-I)

    def sqrt(self):
        """Returns square root of self."""
        from numbers import Rational
        return (self**(Rational(1)/2)).eval()

        
    def print_tree(self):
        """The canonical tree representation
        """
        return str(self)
