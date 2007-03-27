"""Base class for all objects in sympy"""

from sympy.core import hashing

class AutomaticEvaluationType(type):
    """Metaclass for all objects in sympy
     It evaluates the object just after creation, so that for example
     x+2*x becomes 3*x
     """
    def __call__(self,*args,**kwargs):
        if kwargs.has_key("evaluate"):
            evaluate = kwargs["evaluate"]
            del kwargs["evaluate"]
        else:
            evaluate = True
        obj = type.__call__(self,*args,**kwargs)
        if evaluate: 
            return obj.eval()
        else: 
            return obj
        

class Basic(object):
    """
    Base class for all objects in sympy
    
    possible assumptions are: 
        
        - is_real
        
        - is_commutative
        
        - is_bounded
        
    Assumptions can have 3 possible values: 
    
        - True, when we are sure about a property. For example, when we are
        working only with real numbers:
        >>> from sympy import *
        >>> Symbol('x', is_real = True)
        x
        
        - False
        
        - None (if you don't know if the property is True or false)
    """

    __metaclass__ = AutomaticEvaluationType
    
    @property
    def mathml_tag(self):
        """Return the mathml tag of the current object. 
        
        For example, if symbol x has a mathml representation as::
        
           <ci>x</ci>
           
        then x.mathml should return "ci"

        Basic.mathml_tag() returns the class name as the mathml_tag, this is
        the case sometimes (sin, cos, exp, etc.). Otherwise just override this
        method in your class.
        """
        
        return self.__class__.__name__.lower()
    
    def __init__(self, *args, **kwargs):
        self._assumptions = {
                 'is_real' : None, 
                 'is_integer' : None,
                 'is_commutative' : None, 
                 'is_bounded' : None, 
                 }
        self._mhash = 0
        self._args = []
        for k in kwargs.keys():
            if self._assumptions.has_key(k):
                self._assumptions[k] = kwargs[k]
            else:
                raise NotImplementedError ( "Assumption not implemented" )
        
    def __add__(self,a):
        from addmul import Add
        return Add(self, self.sympify(a))
    
    def __radd__(self,a):
        from addmul import Add
        return Add(self.sympify(a), self)
        
    def __getattr__(self, name):
        if self._assumptions.has_key(name):
            return self._assumptions[name]
        else:
            raise AttributeError("'%s' object has no attribute '%s'"%
                (self.__class__.__name__, name))
    
    def __getitem__(self, iter):
        return self._args[iter]
    
    def __len__(self):
        return len(self._args)
    
    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(type(self))
    
    def __neg__(self):
        from numbers import Rational
        return self._domul(Rational(-1),self)
        
    def __pos__(self):
        return self
    
    def __float__(self):
        return float(self.evalf())
 
    def __abs__(self):
        """Returns the absolute value of self. 
        
        Example usage: 
          >>> from sympy import *
          >>> abs(1+2*I)
          5**(1/2)
          >>> x = Symbol('x')
          >>> abs(-x)
          abs(x)
        """
        from functions import abs_
        return abs_(self)
        
    def __radd__(self,a):
        return self._doadd(a, self)
        
    def __sub__(self,a):
        return self._doadd(self, -a)
        
    def __rsub__(self,a):
        return self._doadd(a, -self)
        
    def __mul__(self,a):
        try:
            a=self.sympify(a)
        except:
            return a.__rmul__(self)
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
        if a == None: return False
        return self.isequal(self.sympify(a))
        
    def __ne__(self,a):
        return not self.__eq__(a)
        
    def __lt__(self,a):
        from sympy.core.numbers import Real
        if self._isnumber(self) and self._isnumber(a): 
            return self.evalf() < Real(a).evalf()
        else:
            raise NotImplementedError("'<' not supported.")
        
    def __gt__(self,a):
        from numbers import Real
        if self._isnumber(self) and self._isnumber(a): 
            return self.evalf() > Real(a).evalf()
        else:
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
        
        The eval() method should alway return a new object (following the
        general rule of not changing)
        
        """
        return self
   
    def evalf(self):
        raise ValueError

    def evalc(self):
        """Rewrites self in the form x+i*y.

        It should raise an exceptin, if this is not possible.
        
        """
        raise NotImplementedError

    def get_re_im(self):
        """Returns (x,y) where self=x+i*y""" 
        from numbers import I
        e=self.evalc()
        x = e.subs(I,0)
        y = (e+(-x).expand()).subs(I,1)
        return x,y
 
    @staticmethod
    def sympify(a):
        """for "a" int, returns Rational(a), for "a" float returns real, 
        otherwise "a" (=it's a Basic subclass)."""
        import decimal
        from numbers import Rational, Real
        
        if isinstance(a,int) or isinstance(a, long):
            return Rational(a)
        elif isinstance(a,(float, decimal.Decimal)):
            return Real(a)
        else:
            assert isinstance(a,Basic)
            return a
        
    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        return self._mhash.value
        
    def isequal(self,a):
        return self.hash() == (self.sympify(a)).hash()
        
    @staticmethod
    def cmphash(a,b):
        return Basic._sign(a.hash()-b.hash())
        
    def diffn(self,sym,n):
        while n:
            self = self.diff(sym)
            n -= 1
        return self
        
    def series(self,sym,n):
        from numbers import Rational
        from symbol import Symbol
        from functions import log
        w=Symbol("l",dummy=True)
        f = self.subs(log(sym),-w)
        e = f.subs(sym,Rational(0))
        fact = Rational(1)
        for i in range(1,n+1):
            fact *= Rational(i)
            f = f.diff(sym)
            e += f.subs(sym,Rational(0))*(sym**i)/fact
        e=e.subs(w,-log(sym))
        return e
        
    def subs(self,old,new):
        if self == old:
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
        #TODO: move out of Basic, maybe to polynomials.py?
        
        from numbers import Rational
        from power import Pow
        from addmul import Add,Mul
        from symbol import Symbol
        
        def domul(x):
            if len(x) > 1:
                return Mul(*x)
            return x[0]
        
        def extract(t,x):
            """Parses "t(x)", which is expected to be in the form c0*x^e0,
            and returns (c0,e0). It raises an exception, if "t(x)" is not
            in this form.
            """
            if not t.has(x):
                return t,Rational(0)
            if isinstance(t,Pow):
                return  Rational(1),  t.exp
            elif isinstance(t,Symbol):
                return  Rational(1),  Rational(1)
            assert isinstance(t,Mul)
            for i,a in enumerate(t._args):
                if a.has(x):
                    if isinstance(a,Pow):
                        return  domul(t[:i] + t[i+1:]),  a.exp
                    if isinstance(a,Symbol):
                        return  domul(t[:i] + t[i+1:]),  Rational(1)
                    assert False
            return t,s.Rational(0)
        
        if not isinstance(self,Add):
            return extract(self,x)
        lowest = [0,(Rational(10)**10)]
        l = Symbol("l",dummy=True)
        from functions import log
        for t in self[:]:
            t2 = extract(t.subs(log(x),-l),x)
            if (lowest[1] - t2[1]).evalf()>0:
                lowest=t2
            elif t2[1] == lowest[1]:
                lowest = ((lowest[0] + t2[0]),lowest[1])
        return lowest[0].subs(l,-log(x)), lowest[1].subs(l,-log(x))
        
    def ldegree(self,sym):
        """Returns the lowest power of the sym
        """
        #TODO: move out of Basic
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
        #TODO: move to functions
        from numbers import Rational
        return (self**(Rational(1)/2))

    @staticmethod
    def muleval(x,y):
        """See
        http://groups.google.com/group/sympy/browse_thread/thread/aadbef6e2a4ae335

        Try to simplify x*y. You can either return a simplified expression
        or None.
        
        
        """
        return None

    @staticmethod
    def addeval(x,y):
        """
        Try to simplify x+y in this order. You can either return a simplified
        expression or None.
        
        """
        return None

    def isnumber(self):
        """Return True if self is a number. False otherwise. 
        """
        
        return not 'ci' in self.mathml
        # ci is the mathml notation for symbol, so we assume that 
        # if it's mathml has not the ci tag, then it has no symbols

    @property
    def mathml(self):
        """Returns a MathML expression representing the current object"""
        return "<%s> %s </%s>" % (self.mathml_tag, str(self), self.mathml_tag)
    
    def print_tree(self):
        """The canonical tree representation"""
        return str(self)
    
    def atoms(self, s = [], type=None):
        """Returns the atoms (objects of length 1) that form current
        object. 
        
        Example: 
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> (x+y**2+ 2*x*y).atoms()
        [y, 2, x]
        
        You can also filter the results by a given type of object
        >>> (x+y+2+y**2*sin(x)).atoms(type=sin)
        [sin(x)]
        
        >>> (x+y+2+y**2*sin(x)).atoms(type=Symbol)
        [y, x]
        
        >>> (x+y+2+y**2*sin(x)).atoms(type=Number)
        [2]
        """
        s_temp = s[:] # make a copy to avoid collision with global s
        for arg in self:
            if len(arg) == 1:
                if not arg in s_temp:
                    s_temp.append(arg)
            else:
                # recursive
                s_temp = arg.atoms(s_temp)
        if type is not None:
            # sort only the atoms of a given type
            return filter(lambda x : isinstance(x, type), s_temp)
        return s_temp

    @staticmethod
    def _isnumber(x):
        # TODO: remove
        #don't use this function. Use x.isnumber() instead
        from numbers import Number
        from basic import Basic
        from decimal import Decimal
        if isinstance(x, (Number, int, float, long, Decimal)):
            return True
        assert isinstance(x, Basic)
        return x.isnumber()
    
    @staticmethod
    def _sign(x):
        """Return the sign of x, that is, 
        1 if x is positive, 0 if x == 0 and -1 if x is negative
        """
        if x < 0: return -1
        elif x==0: return 0
        else: return 1

    def match(self, pattern, syms):
        if len(syms) == 1:
            if pattern == syms[0]:
                return self
        if type(self) != type(pattern):
            return None
        for a,b in zip(self,pattern):
            r = a.match(b, syms)
            if r==None:
                return None
        return r
