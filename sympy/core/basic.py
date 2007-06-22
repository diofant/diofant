"""Base class for all objects in sympy"""

from sympy.core.hashing import mhash


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
    """Base class for all objects in SymPy. It is possible to specify
       objects' behaviour using assumptions mechanism. This is most
       useful in pure symbolic algorithms which eg. would be able to
       perform more specific rewriting to intermediate expressions,
       giving simplified and probably more readable results.

       Possible assumptions are:

        - is_number, is_real, is_integer

        - is_negative, is_positive

        - is_nonnegative, is_nonpositive

        - is_odd, is_even, is_prime

        - is_nonzero

        - is_commutative

        - is_bounded

        - is_dummy

       Assumptions are defined in triple-valued logic (in Python sense)
       using True, when we are certain that property applies to given
       object, None, when we aren't sure and False otherwise:

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> sin(x).is_bounded
       True

       >>> print Symbol('x').is_negative
       None

       >>> Rational(1, 2).is_integer
       False

       To assume something about SymPy's object, it is necessary to put
       appropriate flags in selectd object's constructor:

       >>> Symbol('k', integer = True).is_integer
       True
       >>> Symbol('k', is_integer = True).is_integer
       True

       Once object is created and assumptions set it is not allowed to
       change their value without creating new object.

       Sometimes the list of required assumtions can be really long. In
       this case you can use abbreviations to reduce the amount of typing:

       >>> k = Symbol('k', integer = True, nonnegative = True)
       >>> k.is_integer and k.is_nonnegative
       True
       >>> k.is_nonnegative_integer
       True

       >>> k = Symbol('k', nni = True)
       >>> k.is_integer and k.is_nonnegative
       True
       >>> k.is_nonnegative_integer
       True

    """

    __metaclass__ = AutomaticEvaluationType

    def __init__(self, *args, **kwargs):
        self._mhash = 0
        self._mathml = None
        self._args = []

        self._assumptions = {
            'is_real'        : None,
            'is_integer'     : None,
            'is_negative'    : None,
            'is_positive'    : None,
            'is_nonnegative' : None,
            'is_nonpositive' : None,
            'is_nonzero'     : None,
            'is_commutative' : None,
            'is_bounded'     : None,
            'is_dummy'       : None,
            'is_prime'       : None,
            'is_odd'         : None,
            'is_even'        : None,
            }

        dependencies = {
            'is_integer'     : lambda x: { 'is_real'        : x },

            'is_negative'    : lambda x: { 'is_positive'    : (not x) and None,
                                           'is_nonnegative' : not x },

            'is_positive'    : lambda x: { 'is_negative'    : (not x) and None,
                                           'is_nonpositive' : not x },

            'is_nonnegative' : lambda x: { 'is_negative'    : not x },
            'is_nonpositive' : lambda x: { 'is_positive'    : not x },

            'is_odd'         : lambda x: { 'is_integer'     : x or None,
                                           'is_real'        : x or None,
                                           'is_even'        : (not x) and None },

            'is_even'        : lambda x: { 'is_integer'     : x or None,
                                           'is_real'        : x or None,
                                           'is_odd'         : (not x) and None },

            'is_prime'       : lambda x: { 'is_integer'     : x or None,
                                           'is_real'        : x or None,
                                           'is_positive'    : x or None,
                                           'is_negative'    : (not x) and None,
                                           'is_nonpositive' : (not x) and None },

            }

        abbreviations = {
            'is_nni'         : [ 'is_integer', 'is_nonnegative' ],
            'is_npi'         : [ 'is_integer', 'is_nonpositive' ],
            }

        def update_assumptions(key, value):
            if value is not None:
                if dependencies.has_key(key):
                    self._assumptions.update(dependencies[key](value))

                self._assumptions[key] = value

        for key in kwargs.keys():
            value = kwargs[key]

            # this is for compatibility reason
            if key[0:3] != 'is_':
                key = 'is_' + key

            if self._assumptions.has_key(key):
                update_assumptions(key, value)
            elif abbreviations.has_key(key):
                for key in abbreviations[key]:
                    update_assumptions(key, value)
            else:
                raise NotImplementedError("Assumption %s not implemented" % str(key))

    def __getattr__(self, name):
        if self._assumptions.has_key(name):
            return self._assumptions[name]
        else:
            raise AttributeError("'%s' object has no attribute '%s'"%
                (self.__class__.__name__, name))

    def __setattr__(self, name, value):
        if name[0:3] != 'is_':
            object.__setattr__(self, name, value)
        else:
             raise AttributeError("Modification of assumptions is not allowed")

    @property
    def is_zero(self):
        return None

    @property
    def is_unit(self):
        return None

    @property
    def is_negative_integer(self):
        return self.is_integer and self.is_negative

    @property
    def is_nonnegative_integer(self):
        return self.is_integer and self.is_nonnegative

    @property
    def is_positive_integer(self):
        return self.is_integer and self.is_positive

    @property
    def is_nonpositive_integer(self):
        return self.is_integer and self.is_nonpositive

    @property
    def is_number(self):
        """Return True if self is a number or False otherwise."""

        from sympy.core.symbol import Symbol
        return self.atoms(type=Symbol) == []

    def __add__(self,a):
        from addmul import Add
        return Add(self, self.sympify(a))

    def __radd__(self,a):
        from addmul import Add
        return Add(self.sympify(a), self)

    def __div__(self,a):
        from numbers import Rational
        return self._domul(self,self._dopow(a,Rational(-1)))

    def __rdiv__(self,a):
        from numbers import Rational
        return self._domul(self.sympify(self) ** Rational(-1), self.sympify(a))

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

    def __pretty__(self):
        """Make representation as prettyForm: to be overridden
        for everything that looks better in 2D.
        """
        from sympy.core.stringPict import prettyForm
        return prettyForm('[%s]'%self)

    def __pow__(self,a):
        return self._dopow(self, a)

    def __rpow__(self,a):
        return self._dopow(a, self)

    def __eq__(self, a):
        """Return new equation with lhs set to 'self' and rhs set to 'a'.
           If You wanted just an boolean equality test then Equation class
           will handle it in __nonzero__ and __eq__ routines.
        """

        if a is None:
            return False
        else:
            from sympy.modules.solvers import Equation
            return Equation(self, self.sympify(a))

    def __ne__(self, a):
        if a is None:
            return True
        else:
            return self.hash() != self.sympify(a).hash()

    def __lt__(self, a):
        """Make boolean comparison when both 'self' and 'a' are numbers or
           else return new inequlaity with lhs and rhs set respectively. The
           same applies to __le__, __gt__ and __ge__.
        """

        if a is None:
            return False
        else:
            other = Basic.sympify(a)

            if self.is_number and other.is_number:
                return self.evalf() < other.evalf()
            else:
                from sympy.modules.solvers import Inequality
                return Inequality(self, other)

    def __le__(self, a):
        if a is None:
            return False
        else:
            other = Basic.sympify(a)

            if self.is_number and other.is_number:
                return self.evalf() <= other.evalf()
            else:
                from sympy.modules.solvers import StrictInequality
                return StrictInequality(self, other)

    def __gt__(self, a):
        if a is None:
            return False
        else:
            other = Basic.sympify(a)

            if self.is_number and other.is_number:
                return self.evalf() > other.evalf()
            else:
                from sympy.modules.solvers import Inequality
                return Inequality(other, self)

    def __ge__(self, a):
        if a is None:
            return False
        else:
            other = Basic.sympify(a)

            if self.is_number and other.is_number:
                return self.evalf() >= other.evalf()
            else:
                from sympy.modules.solvers import StrictInequality
                return StrictInequality(other, self)

    @property
    def mathml_tag(self):
        """Return the mathml tag of the current object.

        For example, if symbol x has a mathml representation as::

           <ci>x</ci>

        then x.mathml_tag should return "ci"

        Basic.mathml_tag() returns the class name as the mathml_tag, this is
        the case sometimes (sin, cos, exp, etc.). Otherwise just override this
        method in your class.
        """

        return self.__class__.__name__.lower()

    def __mathml__(self):
        """Returns a MathML expression representing the current object"""
        import xml.dom.minidom
        if self._mathml:
            return self._mathml
        dom = xml.dom.minidom.Document()
        x = dom.createElement(self.mathml_tag)
        for arg in self._args:
            x.appendChild( arg.__mathml__() )
        self._mathml = x

        return self._mathml


    def __latex__(self):
        return str(self) # override this for a custom latex representation

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

        It will raise an exception, if this is not possible.

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
        elif isinstance(a, complex):
            from numbers import I
            real,imag = Basic.sympify(a.real), Basic.sympify(a.imag)
            ireal, iimag = int(real), int(imag)
            t = ireal+iimag*1j
            if t == a:
                return ireal + iimag*I
            else:
                return real + imag*I
        else:
            if not isinstance(a,Basic):
                raise ValueError("%s must be a subclass of basic" % str(a))
            return a

    def __hash__(self):
        return hash(str(self.hash()))
        # needed by sets and other python functions
        # we do not return .hash() since this is a long

    def hash(self):
        if self._mhash:
            return self._mhash.value
        self._mhash = mhash()
        self._mhash.addstr(self.__class__.__name__)
        for item in self[:]:
            if isinstance(item, Basic):
                self._mhash.add(item.hash())
            else:
                self._mhash.addstr(str(item))
        if self.is_dummy:
            self._mhash.value += 1
        return self._mhash.value

    @staticmethod
    def cmphash(a,b):
        return Basic._sign(a.hash()-b.hash())

    def series(self,sym, n=6):
        """
        Usage
        =====
            Return the Taylor series around 0+ (i.e. 0 from the right) of self
            with respect to sym until the n-th term. Use substitution if you
            want to get a series around a different point or from the left.

        Notes
        =====
            If you don't specify n, the default is 6

        Examples
        ========

            >>> from sympy import *
            >>> x = Symbol('x')
            >>> sin(x).series(x, 5)
            O(x**5)+x-1/6*x**3
        """
        from numbers import Rational
        from symbol import Symbol, O
        from functions import log
        w=Symbol("l", dummy=True)
        f = self.subs(log(sym),-w)
        e = f.subs(sym,Rational(0))
        fact = Rational(1)
        for i in range(1,n):
            fact *= Rational(i)
            #print "x1",f
            f = f.diff(sym)
            #print "x2",f
            e += f.subs(sym,Rational(0))*(sym**i)/fact
        e=e.subs(w,-log(sym))
        #print self,e
        if e == self:
            return e
        else:
            return e+O(sym**n)

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
        n = Symbol("dummy", dummy = True)
        return self.subs(sub,n)!=self

    def has_any(self, subs):
        return len ([sub for sub in subs if self.has(sub) == True]) > 0

    def leadterm(self,x):
        """Returns the leading term c0*x^e0 of the power series 'self' in x
        with the lowest power of x in a form (c0,e0)
        """
        #TODO: move out of Basic, maybe to polynomials.py?

        from numbers import Rational
        from power import Pow
        from addmul import Add,Mul
        from symbol import Symbol, O
        if isinstance(self,Add):
            self = self.removeO()
        if isinstance(self,O):
            self = Rational(0)

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
        l = Symbol("l", dummy=True)
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

    def combine(self):
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

        if isinstance(self, atoms_class):
            if type is None:
                return [self]
            return filter(lambda x : isinstance(x, type), [self])

        if isinstance(self, atoms_class):
            if type:
                if not isinstance(self, type):
                    return []
            return [self]

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
        if isinstance(pattern, Symbol):
            # case pattern is just a symbol: there's nothing to match
            return {pattern: self}
        if type(self) != type(pattern):
            from addmul import Mul
            from numbers import Rational
            if isinstance(pattern, Mul):
                return Mul(Rational(1),self,
                        evaluate=False).match(pattern,syms)
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

