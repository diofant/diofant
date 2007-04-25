"""Base class for all objects in sympy"""

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
        
        - is_number
        
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
                 'is_dummy' : None, 
                 }
        self._mhash = 0
        self._args = []
        for k in kwargs.keys():
            #TODO: maybe use dict.update() for this
            if self._assumptions.has_key(k):
                self._assumptions[k] = kwargs[k]
            else:
                raise NotImplementedError ( "Assumption %s not implemented" % str(k))
        
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
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return self.__class__.__name__ + "(" + str(self[:])[1:-1] + ")"
    
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
        """Test for equality of two expressions
        Currently this is done by comparing the hashes of the two expressions
        """
        if a is None: 
            return False
        return hash(self) == hash(self.sympify(a))
        
    def __ne__(self,a):
        return not self.__eq__(a)
        
    def __lt__(self,a):
        from sympy.core.numbers import Real
        a = Basic.sympify(a)
        if self.is_number and a.is_number: 
            return self.evalf() < a.evalf()
        else:
            raise NotImplementedError("'<' not supported.")
        
    def __gt__(self,a):
        a = Basic.sympify(a)
        if self.is_number and a.is_number: 
            return self.evalf() > a.evalf()
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
        """
        Usage
        =====
            Converts an arbitrary expression to a type that can be used
            inside sympy. For example, it will convert python int's into
            instance of sympy.Rational, floats into intances of sympy.Real, etc.
            
        Notes
        =====
            It currently accepts as arguments: 
                - any object defined in sympy (except maybe matrices [TODO])
                - standard numeric python types: int, long, float, Decimal
                - strings (like "0.09" or "2e-19")
            
            If the argument is already a type that sympy understands, it will do
            nothing but return that value - this can be used at the begining of a
            method to ensure you are workint with the corerct type. 
            
        Examples
        ========
            >>> def is_real(a):
            ...     a = Basic.sympify(a)
            ...     return a.is_real
            >>> is_real(2.0)
            True
            >>> is_real(2)
            True
            >>> is_real("2.0")
            True
            >>> is_real("2e-45")
            True
        """
        
        if isinstance(a, Basic):
            #most common case
            return a
        
        import decimal
        from numbers import Rational, Real
        
        if isinstance(a,int) or isinstance(a, long):
            return Rational(a)
        elif isinstance(a,(float, decimal.Decimal, str)):
            return Real(a)
        else:
            assert isinstance(a,Basic)
            return a
        
    def __hash__(self):
        return self.hash()
    
    def hash(self):
        if self._mhash: 
            return self._mhash
        self._mhash = hash(str(self))
        return self._mhash
    
        
    def isequal(self,a):
        return self.hash() == (self.sympify(a)).hash()
        
    @staticmethod
    def cmphash(a,b):
        return Basic._sign(a.hash()-b.hash())
        
    def series(self,sym, n=6):
        """
        Usage
        =====
            Return the Taylor series of self with respect to sym until the n-th term. 
        
        Notes
        =====
            If you don't specify n, the default is 6
        
        Examples
        ========
        
            >>> from sympy import *
            >>> x = Symbol('x')
            >>> sin(x).series(x, 5)
            x-1/6*x**3+1/120*x**5
        """
        from numbers import Rational
        from symbol import Symbol
        from functions import log
        w=Symbol("l", is_dummy=True)
        f = self.subs(log(sym),-w)
        e = f.subs(sym,Rational(0))
        fact = Rational(1)
        for i in range(1,n+1):
            fact *= Rational(i)
            f = f.diff(sym)
            e += f.subs(sym,Rational(0))*(sym**i)/fact
        e=e.subs(w,-log(sym))
        return e

    def subs_dict(self, di):
        """Substitutes all old -> new defined in the dictionary "di"."""
        x = self
        for d in di:
            x = x.subs(d,di[d])
        return x
        
    def subs(self,old,new):
        """Substitutes an expression old -> new."""
        if self == old:
            return self.sympify(new)
        else:
            return self
            
    def has(self,sub):
        from symbol import Symbol
        n = Symbol("dummy", is_dummy = True)
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
        l = Symbol("l", is_dummy=True)
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

    @staticmethod
    def muleval(x,y):
        """
        Usage
        =====
            Try to simplify x*y. You can either return a simplified expression
            or None.
        Notes
        =====
            This method is called from Add.eval() and Mul.eval(), so that means it
            is called each time you create an instance of those classes
    
            See
            http://groups.google.com/group/sympy/browse_thread/thread/aadbef6e2a4ae335
            
        See also
        ========
            L{sympy.addmul.Add.eval}
            L{sympy.addmul.Mul.eval}
            L{sympy.power.Pow.eval}
        """
        return None

    @staticmethod
    def addeval(x,y):
        """
        Try to simplify x+y in this order. You can either return a simplified
        expression or None.
        
        see docs for muleval
        
        """
        return None

    @property
    def is_number(self):
        """Return True if self is a number. False otherwise. 
        """
        
        return not 'ci' in self.mathml
        # ci is the mathml notation for symbol, so we assume that 
        # if it's mathml has not the ci tag, then it has no symbols

    @property
    def mathml(self):
        """Returns a MathML expression representing the current object"""
        assumptions = ""
        for a in self._assumptions: 
            if self._assumptions[a] is not None:
                assumptions += str(a) + ":" + str(self._assumptions[a]) + ";"
        return "<%s sympy:assumptions='%s'> %s </%s>" % (self.mathml_tag, assumptions, str(self), self.mathml_tag)
    
    def print_tree(self):
        """The canonical tree representation"""
        return str(self)
    
    def atoms(self, s = [], type=None):
        """Returns the atoms that form current
        object. 
        
        Example: 
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> (x+y**2+ 2*x*y).atoms()
        [2, x, y]
        
        You can also filter the results by a given type of object
        >>> (x+y+2+y**2*sin(x)).atoms(type=Symbol)
        [x, y]
        
        >>> (x+y+2+y**2*sin(x)).atoms(type=Number)
        [2]
        """
        from sympy.core.numbers import Number
        from sympy.core.symbol import Symbol

        atoms_class = (Number, Symbol)

        s_temp = s[:] # make a copy to avoid collision with global s
        for arg in self:
            if isinstance(arg, atoms_class):
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
        #don't use this function. Use x.is_number instead
        from numbers import Number
        from basic import Basic
        from decimal import Decimal
        if isinstance(x, (Number, int, float, long, Decimal)):
            return True
        assert isinstance(x, Basic)
        return x.is_number
    
    @staticmethod
    def _sign(x):
        """Return the sign of x, that is, 
        1 if x is positive, 0 if x == 0 and -1 if x is negative
        """
        if x < 0: return -1
        elif x==0: return 0
        else: return 1

    def doit(self):
        """Calls recursively doit() on every element in the expression tree. """
        x = [a.doit() for a in self]
        e = (type(self))(*x)
        return e


    def match(self, pattern, syms=None):
        #print "B",self,pattern,syms,type(self),type(pattern)
        from symbol import Symbol
        if syms == None:
            syms = pattern.atoms(type=Symbol)
        if len(syms) == 1:
            if pattern == syms[0]:
                return {syms[0]: self}
        if type(self) != type(pattern):
            from addmul import Mul
            from numbers import Rational
            if isinstance(pattern, Mul):
                return Mul(Rational(1),self,
                        evaluate = False).match(pattern,syms)
            return None
        r2 = None
        #print "aaaa",self,pattern
        for a,b in zip(self,pattern):
            r = a.match(b, syms)
            #print "A",a,b,syms,"-->",r
            if r==None:
                return None
            if r2 == None:
                r2 = r
            else:
                r2.update(r)
        return r2

    def __pretty__(self):
        """Make representation as prettyForm: to be overridden
        for everything that looks better in 2D.
        """
        from sympy.core.stringPict import prettyForm
        return prettyForm('[%s]'%self)