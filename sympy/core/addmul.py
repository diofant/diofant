"""
This module has the necessary code to represent addition and multiplication of elements, via
the classes Add, Mul and it's base class, Pair. 

This is a central part of the core
"""

import hashing
from sympy.core.basic import Basic
from sympy.core.numbers import Number, Rational, Real
from sympy.core.power import Pow, pole_error

class Pair(Basic):
    """Abstract class containing common code to add and mul classes.
    Should not be used directly
    """
    
    def __init__(self, *args):
        Basic.__init__(self)
        for arg in args:
            assert isinstance(arg, Basic)
        self._args = args

    def __lt__(self, a):
        return self.evalf() < a
    
    @property
    def mathml(self):
        s = "<apply>" + "<" + self.mathml_tag + "/>"
        for a in self._args:
            s += a.mathml
        s += "</apply>"
        return s
    
    def hash(self):
        if self._mhash:
            return self._mhash.value
        self._mhash = hashing.mhash()
        self._mhash.addstr(str(type(self)))
        for i in self[:]:
            self._mhash.add(i.hash())
        return self._mhash.value

    
    def tryexpand(self, a):
        if isinstance(a,Mul) or isinstance(a,Pow):
            return a.expand()
        else:
            return a
            
    def flatten(self, a ):
        """::
            flatten([add(x,4),Mul(a,5),add(x,b),x]) ->
                  [x,4,Mul(a,5),x,b,x] if self is add
                  [add(x,4),a,5,add(x,b),x] if self is Mul

        returns a copy of "a", where the the classes of the same type as this
        (e.g. if we are Mul, then all the Muls, if we are add, then all the
        adds) are substituted for their arguments. all other classes are
        left intact.
        """
        b=[]
        for x in a:
            if isinstance(x,type(self)):
                b.extend(x[:])
            else:
                b.append(x)
        return b
        
    @staticmethod
    def coerce(a, action):
        """coerce([x,y,z],action) -> action(action(action([],x),y),z)"""
        #equivalent code:
        #exp=[]
        #for x in a:
        #    exp=action(exp,x)
        #return exp
        return reduce(action, a, [])
        
    @staticmethod
    def coerce_numbers(a, action, default):
        """coercenumbers([x,4,a,10],action,Rational(1)) -> 
        (action(action(Rational(1),4),10),[x,a])

        picks out the numbers of the list "a" and applies the action on them
        (add or Mul).
        """
        n = default
        b = []
        for x in a:
            if isinstance(x,Number):
                n = action(n,x)
            else:
                b.append(x)
        return (n,b)
    
    def print_tree(self):
        def indent(s,type=1):
            x = s.split("\n")
            r = "+--%s\n"%x[0]
            for a in x[1:]:
                if a=="": continue
                if type==1:
                    r += "|  %s\n"%a
                else:
                    r += "   %s\n"%a
            return r
        if isinstance(self,Mul):
            f="Mul\n"
        else:
            assert(isinstance(self,Add))
            f="Add\n"
        for a in self._args[:-1]:
            f += indent(a.print_tree())
        f += indent(self._args[-1].print_tree(),2)
        return f
    
    @property        
    def is_commutative(self):
        for x in self[:]:
            if not x.is_commutative:
                return False
        return True


class Mul(Pair):
    
    mathml_tag = "times"
     
    def __str__(self):
        f = ""
        a = self._args
        if isinstance(a[0],Rational):
            if a[0].isminusone():
                f = "-"
                a = self._args[1:]
            elif a[0].isone():
                f = ""
                a = self._args[1:]
        for x in a:
            if isinstance(x,Pair):
                f += "(%s)*"
            else:
                f += "%s*"
        f = f[:-1]
        return f % tuple([str(x) for x in a])
    

    @staticmethod
    def get_baseandexp(a):
        # TODO: remove
        if isinstance(a,Pow):
            return a[:]
        else:
            return (a,Rational(1))

    @staticmethod
    def try_to_coerce(x, y):
        """Tries to multiply x * y in this order and see if it simplifies. 
        
        If it succeeds, returns (x*y, True)
        otherwise (x, False)
        where x is the original x 
        """
        #TODO: See also: L{Add.eval}
        z1 = y.muleval(x,y)
        z2 = x.muleval(x,y)

        if z1 or z2:
            if (z1 and z2):
                #sanity check
                assert z1==z2
            if z1:
                return z1, True
            if z2:
                return z2, True

        if isinstance(x,Number) and isinstance(y, Number):
            return x*y, True
        xbase,xexp = Mul.get_baseandexp(x)
        ybase,yexp = Mul.get_baseandexp(y)
        if xbase.isequal(ybase):
            #this whole "if" is to correctly cooperate with Pow.eval()
            #so we don't get infinite recursion. It's not elegant, but it
            #works.
            if not xbase.is_commutative:
                e = Add(xexp,yexp)
                if e != 0:
                    return Pow(xbase,e,evaluate=False), True

            return Pow(xbase,Add(xexp,yexp)), True
        else:
            return x, False
            
    def eval(self):
        "Flatten, put all Rationals in the front, sort arguments"

        def _mul_c(exp,x):
            e = []
            for i,y in enumerate(exp):
                z,ok = self.try_to_coerce(y,x)
                if isinstance(z, Number) and i != 0:
                    #c and 1/c could have been coerced to 1 or i^2 to -1
                    #or 2^(1/2)^2 to 2, etc.
                    #z == 0 is probably a bug
                    assert z != 0
                    e[0] *= z
                else:
                    e.append(z)
                if ok: 
                    e.extend(exp[i+1:])
                    return e
            e.append(x)
            return e

        def _mul_nc(exp,x):
            if exp == []: return [x]
            #try to join only last and the one before last object
            z,ok = self.try_to_coerce(exp[-1], x)
            if ok:
                return exp[:-1]+[z]
            else:
                return exp[:-1]+[z]+[x]

        def _nc_separate(a):
            """Separates commutative and non-commutative parts of "a" """
            c_part = []
            nc_part = []
            for x in a:
                if x.is_commutative:
                    c_part.append(x)
                else:
                    nc_part.append(x)
            return c_part, nc_part

        #(((a*4)*b)*a)*5  -> a*4*b*a*5:
        a = self.flatten(self._args)
        #separate C and NC parts
        c_part, nc_part = _nc_separate(a)
        nc_part_tmp = self.coerce(nc_part,_mul_nc)
        #the coerce method could generate some C things, or nested Muls, so
        #flatten and separate C and NC parts again
        nc_part_tmp = self.flatten(nc_part_tmp)
        c_part2, nc_part = _nc_separate(nc_part_tmp)
        #we put 1 in front of everything
        a = self.coerce([Rational(1)]+c_part+c_part2,_mul_c)
        n,c_part = a[0], a[1:]
        #so that now "n" is a Number and "c_part" doesn't contain any number
        if n == 0: return Rational(0)
        c_part.sort(Basic.cmphash)
        #this if is for multiplying Symbol*Matrix and Number*Matrix
        #I think it's not needed anymore, since Matrix is implemented
        #differently now. But I am leaving it here for now, because
        #it works and the noncommutative objects are not well tested yet.
        if len(nc_part) == 1:
            if len(c_part) == 1:
                z, ok = self.try_to_coerce(c_part[0], nc_part[0])
                if ok: return z
            if len(c_part) == 0:
                z, ok = self.try_to_coerce(n, nc_part[0])
                if ok: return z
        a=c_part+nc_part
        #put the number in front of all the other args
        if n != 1: a=[n]+a
        if len(a) > 1:
            #construct self again, but evaluated this time
            return Mul(evaluate=False, *a)
        elif len(a) == 1:
            return a[0]
        else:
            return Rational(1)
            
    def evalf(self):
        a, b = self.getab()
        if a.isnumber() and b.isnumber():
            return Real(a)*Real(b)
        else: 
            raise ValueError("Cannot evaluate a symbolic value")

    def evalc(self):
        a, b = self.getab()
        return (a.evalc() * b.evalc()).expand()
            
    def getab(self):
        """Pretend that self=a*b and return a,b
        
        in general, self=a*b*c*d*..., but in many algorithms, we 
        want to have just 2 arguments to Mul. Use this function to 
        simulate this interface. (the returned b = b*c*d.... )
        """
        a = self._args[0]
        if len(self._args) == 1:
            b = Rational(1)
        elif len(self._args) == 2:
            b = self._args[1]
        else:
            b = Mul(*self._args[1:])
        return (a,b)
        
    def diff(self,sym):
        r = Rational(0)
        for i in range(len(self._args)):
            d = self._args[i].diff(sym)
            for j in range(len(self._args)):
                if i != j:
                    d *= self._args[j]
            r+=d
        return r
        
    def series(self,sym,n):
        """expansion for Mul
        need to handle this correctly:
        (e^x-1)/x  -> 1+x/2+x^2/6
        the first term is (e^x-1)/x evaluated at x=0. normally, we would use
        a limit. but in cas, the limit is computed using series. so we must
        calculate series differently - using the bottom up approach:
        first expand x, then e^x, then e^x-1, and finally (e^x-1)/x
        """
        a,b=self.getab()
        x=a.series(sym,n)
        try:
            y=b.series(sym,n)
        except pole_error:
            #we are not able to expand b, 
            #but if a goes to 0 and b is bounded, 
            #the result is just a*const, so we just return a
            a0 = x.subs(sym,0)
            if a0==0 and b.bounded():
                return x
            #we cannot expand x*y
            raise
        return (x*y).expand()
        
    def expand(self):
        a,b = self.getab()
        a = self.tryexpand(a)
        b = self.tryexpand(b)
        if isinstance(a,Add):
            d = Rational(0)
            for t in a[:]:
                d += (t*b).expand()
            return d
        elif isinstance(b,Add):
            d = Rational(0)
            for t in b[:]:
                d += (a*t).expand()
            return d
        else:
            return a*b
        
    def subs(self,old,new):
        a,b = self.getab()
        e = a.subs(old,new)*b.subs(old,new)
        if isinstance(e, Basic):
            return e
        else:
            return e

    def match(self, pattern, syms):
        if len(self) == 2:
            r = self[0].match(pattern[0], syms)
            if r!=None: return r
            r = self[1].match(pattern[1], syms)
            if r!=None: return r
            r = self[0].match(pattern[1], syms)
            if r!=None: return r
            r = self[1].match(pattern[0], syms)
            if r!=None: return r
            return None
        return Basic.match(self, pattern, syms)


class Add(Pair):
    """
    Usage
    =====
        This class represent's the addition of two elements. so whenever you call '+', an 
        instance of this class is created. 
        
    Notes
    =====
        When an instance of this class is created, the method .eval() is called and will
        preform some inexpensive symplifications. 
        
        In some cases, the eval() method will return an object that is not an instance of the
        class Add, so for example if x is a Symbol, (x+x) will create a class Add with arguments
        (x,x) , that will be evaluated via the .eval() method, and this method will return a 
        class Mul with arguments (2,x), that is how x+x --> 2*x is done
        
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> type(1+x)
        <class 'sympy.core.addmul.Add'>
        >>> (1+x)[:]
        (1, x)
    
    See also
    ========
        L{Add.eval}
    """
    
    mathml_tag = "plus"
    
    def __str__(self):
        """Returns a string representation of the expression in self."""
        
        f = "%s" % str(self._args[0])
        for i in range(1,len(self._args)):
            num_part = _extract_numeric(self._args[i])[0]
            if num_part < 0:
                f += "%s" % str(self._args[i])
            else:
                f += "+%s" % str(self._args[i])
        return f    

                
    def getab(self):
        """Pretend that self = a+b and return a,b
        
        in general, self=a+b+c+d+..., but in many algorithms, we 
        want to ha+ve just 2 arguments to add. Use this function to 
        simulate this interface. (the returned b = b+c+d.... )
        
        If you want to obtain all the arguments of a given expression, use
        the slices syntax, like in (1+x)[:]
        """
        a=self._args[0]
        if len(self._args) == 1:
            b = Rational(0)
        elif len(self._args) == 2:
            b = self._args[1]
        else:
            assert len(self._args) > 2
            b = Add(*self._args[1:])
        return (a,b)
    

        
    def eval(self):
        """
        Usage
        =====
            This method is called automatically when an instance of this class
            is created, and will always return a new object. 
            
        Notes
        =====
            Currently, what this method does is: 
            
                - flatten all arguments, i.e, substitute all instances of Add by their arguments
                  for example, self.flatten(1, 1+x ) --> [1,1,x]
                
                - adds it's arguments beying aware of some identities, like that x+x -> 2*x and 
                  that numbers can be added without restrictions using their own __add__ method
                
                - sort
        See also
        ========
            L{Mul.eval}, L{Pow.eval}
            
        TODO
        ====
            - Perform a complexity analysis
            - probably optimizations can be done (algorithmic optimizations)
        """

        def _add(exp,x):
            an, a = _extract_numeric(x)
            e = []
            ok = False
            for y in exp:
                bn, b = _extract_numeric(y)
                if (not ok) and a.isequal(b):
                    e.append(Mul(an + bn,a))
                    ok = True
                else:
                    z1 = x.addeval(y,x)
                    z2 = y.addeval(y,x)

                    if z1 or z2:
                        if (z1 and z2):
                            #sanity check, only when z1 and z2 are not Order
                            from symbol import Order
                            if not isinstance(z1,Order):
                                assert z1 == z2
                        if z1:
                            e.append(z1)
                            ok = True
                        elif z2:
                            e.append(z2)
                            ok = True
                    else:
                        e.append(y)
            if not ok: e.append(x)
            return e

        def _add_Number(a,b):
            """Adds two Real or Rational Numbers"""
            if isinstance(a,Number):
                return a+b
        
        a = self.flatten(self._args)
        a = self.coerce(a,_add)
        #n,a = self.coerce_numbers(a, Rational.__add__, Rational(0))
        n,a = self.coerce_numbers(a, _add_Number, Rational(0))
        a.sort(Basic.cmphash)
        if n != 0:
            a = [n] + a
        if len(a)>1:
            return Add(evaluate=False, *a)
        elif len(a) == 1:
            return a[0]
        else:
            return Rational(0)
        
    def evalf(self):
        a,b = self.getab()
        if hasattr(a, 'evalf') and hasattr(b, 'evalf'):
            return a.evalf() + b.evalf()
        else:
            raise ValueError('Can not evaluate a symbolic value')

    def evalc(self):
        a, b = self.getab()
        return (a.evalc() + b.evalc()).expand()
    
    def diff(self,sym):
        d = Rational(0)
        for x in self._args:
            d += x.diff(sym)
        return d
    
    def expand(self):
        """Tries to expands all the terms in the sum."""
        d = Rational(0)
        for x in self._args:
            d += self.tryexpand(x)
        return d
    
    def subs(self,old,new):
        d = Rational(0)
        for x in self._args:
            d += x.subs(old,new)
        return d
    
    def series(self,sym,n):
        """expansion for add
        need to handle this correctly:
        x+1/x
        tries to use Basic.series (which substitutes x->0), if it fails,
        expands term by term
        """
        try:
            return Basic.series(self,sym,n)
        except pole_error:
            a,b = self.getab()
            #there is a cancelation problem here:
            #implement the class Order
            return (a.series(sym,n)+b.series(sym,n))

    def match(self, pattern, syms):
        if len(syms) == 1:
            p = syms[0]
            if not pattern.has(p):
                if self == pattern:
                    return {}
                else:
                    return None
            rest = pattern - p
            from symbol import Symbol
            if isinstance(rest, Symbol):
                if rest in self[:]:
                    return {p: self - rest}
                else:
                    return None
            for x in rest[:]:
                if not (x in self[:]):
                    return None
            return {p: (self - rest).expand()}
        return NotImplementedError("More than one pattern symbol")

def _extract_numeric(x):
    """Returns the numeric and symbolic part of x.
    For example, 1*x -> (1,x)
    Works only with simple expressions. 
    """
    if isinstance(x, Mul) and isinstance(x._args[0], Number):
        return x.getab()
    else:
        return (Rational(1), x)
