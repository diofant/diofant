
from sympy.core.hashing import mhash
from sympy.core.basic import Basic
from sympy.core.numbers import Rational
from sympy.core.stringPict import prettyForm

dummycount = 0

class Symbol(Basic):
    """
    Assumptions::
       is_real = True
       is_commutative = True

    You can override the default assumptions in the constructor::
       >>> A = Symbol('A', is_commutative = False)
       >>> B = Symbol('B', is_commutative = False)
       >>> A*B != B*A
       True
       >>> (A*B*2 == 2*A*B) == True # multiplication by scalars is commutative
       True
    """
    
    mathml_tag = "ci"
    
    dummy_num = 0

    def __init__(self, name, commutative=True, dummy=False, real=False, 
                 *args, **kwargs):
        """if dummy == True, then this Symbol is totally unique, i.e.::
        
        >>> (Symbol("x") == Symbol("x")) == True
        True
        
        but with the dummy variable ::
        
        >>> (Symbol("x", dummy = True) == Symbol("x", dummy = True)) == True
        False

        """
        
        self._assumptions = {
                         'is_commutative' : commutative,
                         'is_dummy': dummy, 
                         'is_real': real, 
                         }
        
        for k in kwargs.keys():
            self._assumptions[k] = kwargs[k]
        
        Basic.__init__(self, **self._assumptions)
        self.name = name
        if self.is_dummy:
            global dummycount
            self.dummy_num = dummycount
            dummycount += 1
        #self._args = [name]

    def __str__(self):
        if not self.is_dummy:
            return str(self.name)
        else:
            # if x is dummy
            return str(self.name + '__' + str(self.dummy_num))
        
    def __mathml__(self):
        import xml.dom.minidom
        if self._mathml:
            return self._mathml
        dom = xml.dom.minidom.Document()
        x = dom.createElement(self.mathml_tag)
        x.appendChild(dom.createTextNode(self.name))
        self._mathml = x
        
        return self._mathml

    def hash(self):
        if self._mhash: 
            return self._mhash.value
        self._mhash = mhash()
        self._mhash.addstr(self.__class__.__name__)
        self._mhash.addstr(self.name)
        if self.is_dummy:
            global dummycount
            self._mhash.value += dummycount
            dummycount += 1
        return self._mhash.value
    

    def diff(self,sym):
        if not self.is_commutative:
            raise NotImplementedError("Differentiation of non-commutative objects. " + \
                                      + "Doesn't have a meaning.")
        if self == sym:
            return Rational(1)
        else:
            return Rational(0)

    def evalc(self):
        if self.is_real:
            return self
        raise NotImplementedError

    def doit(self):
        return self

    def match(self, pattern, syms):
        if self == pattern:
            return {}
        if len(syms) == 1:
            if pattern == syms[0]:
                return {syms[0]: self}
            if self == pattern:
                return {}
        if isinstance(pattern, Symbol):
            try:
                return {syms[syms.index(pattern)]: self}
            except ValueError:
                pass
        from addmul import Mul
        if isinstance(pattern, Mul):
            return Mul(Rational(1),self,evaluate = False).match(pattern,syms)
        return None
    
    def __pretty__(self): 
        return prettyForm(self.name, binding=prettyForm.ATOM)


class Order(Basic):
    """
    Represents O(f(x)) at the point x = 0.

    Definition
    ==========

    g(x) = O(f(x)) as x->a  if and only if
    |g(x)|<=M|f(x)| near x=a                     (1)

    In our case Order is for a=0. An equivalent way of saying (1) is:

    lim_{x->a}  |g(x)/f(x)|  < oo
    
    Let's illustrate it on the following example:

    sin x = x - x**3/3! + O(x**5)

    where in this case O(x**5) = x**5/5! - x**7/7! + .... and the definition
    of O means:

    |x**5/5! - x**7/7! + ....| <= M|x**5|      near x=0

    or equivalently:

    lim_{x->0} |x**5/5! - x**7/7! + .... / x**5| < oo

    which surely is true, because 
    
    lim_{x->0} |x**5/5! - x**7/7! + .... / x**5| = 1/5!


    So intuitively O(x**3) means (in our case): all terms x**3, x**4 and
    higher. But not x**2, x or 1.

    Examples:
    =========
    >>> from sympy import *
    >>> x = Symbol("x")
    >>> Order(x)
    O(x)
    >>> Order(x)*x
    O(x**2)
    >>> Order(x)-Order(x)
    O(x)

       External links
       --------------

         U{Big O notation<http://en.wikipedia.org/wiki/Big_O_notation>}
    """

    def __init__(self, f, sym=None):
        """O(f) at the point x = 0"""
        Basic.__init__(self)
        self._args = [self.sympify(f)]
        if sym:
            self.sym = sym
        else:
            self.sym = self._args[0].atoms(type = Symbol)
            if len(self.sym) == 1:
                self.sym = self.sym[0]
            else:
                #well, let's try to guess
                if self.sym == []:
                    self.sym = Rational(1)
                else:
                    raise "Don't know the variable in Order"

    def eval(self):
        from addmul import Mul, Add
        from numbers import Real, Rational
        f = self[0]
        if isinstance(f, Mul):
            if isinstance(f[0], (Real, Rational)):
                assert len(f[:]) == 2
                return Order(f[1])
            if not f[0].has(self.sym):
                assert len(f[:]) == 2
                return Order(f[1])
            if not f[1].has(self.sym):
                assert len(f[:]) == 2
                return Order(f[0])
            e = f.expand()
            if isinstance(e, Add):
                r=0
                for x in e:
                    r+=Order(x)
                return r
        if isinstance(f, Add):
            if isinstance(f[0], (Real, Rational)):
                if len(f[:]) == 2:
                    return Order(f[1])
            r=0
            for x in f:
                r+=Order(x)
            return r
        if isinstance(f, (Real, Rational)) and f!=0 and f!=1:
            return Order(Rational(1))
        return self

    def __str__(self):
        return "O(%s)"%str(self[0])

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Order) and isinstance(y, Order):
            return Order(x[0]*y[0])
        if isinstance(y, Order):
            return Order(x*y[0],sym = y.sym)
        return None

    @staticmethod
    def addeval(x, y):
        from power import pole_error
        if isinstance(x, Order) and isinstance(y, Order):
            if isinstance(x.sym, Symbol):
                sym = x.sym
            else:
                sym = y.sym
            #calculate inf = True if this limit is oo:
            #inf = lim_{x->a}  |g(x)/f(x)| == oo
            inf = False
            try:
                #we don't want to depend on the limit module, thus
                #we use the pole_error way, which works in most cases
                (y[0]/x[0]).subs(sym,0)
            except pole_error:
                inf = True
            #print x,y,inf
            if inf:
                return y
            else:
                return x
        if isinstance(x, Order):
            #calculate inf = True if this limit is oo:
            #inf = lim_{x->a}  |g(x)/f(x)| == oo
            inf = False
            try:
                #we don't want to depend on the limit module, thus
                #we use the pole_error way, which works in most cases
                (y/x[0]).subs(x.sym,0)
            except pole_error:
                inf = True

            if inf:
                return None
            else:
                return Order(x[0])

        if isinstance(y, Order):
            return Order.addeval(y,x)
        return None

    def subs(self,old,new):
        """Substitutes an expression old -> new."""
        e = Basic.subs(self,old,new)
        if e == self:
            if old == self.sym:
                if new == 0:
                    return Rational(0)
		elif isinstance(new, Symbol):
		    return Order(new)
                else:
                    raise ValueError("Cannot substitute (%s, %s) in Order" % (new, old) )
        return e

    def diff(self, var):
        e = self[0].diff(var)
        if e == 0:
            return Order(1)
        else:
            return Order(e)
