"""
This module has the necessary code to represent addition and multiplication of elements, via
the classes Add, Mul and it's base class, Pair. 

This is a central part of the core
"""

from sympy.core.basic import Basic
from sympy.core.numbers import Rational, Real, Infinity
from sympy.core.power import Pow, pole_error
from sympy.core.functions import Function

from sympy.core.stringPict import prettyForm

#from decorator import decorator

#@decorator
#def debug(f,*args,**kw):
    #print "-"*40
    #print "starting %s%s" % (f.func_name,args)
#    r = f(*args,**kw)
    #print "done %s%s = %s" % (f.func_name,args,r)
    #print "-"*40
#    return r

class Pair(Basic):
    """Abstract class containing common code to Add and Mul classes.
    Should not be used directly.
    """
    
    def __init__(self, *args):
        Basic.__init__(self)
        for arg in args:
            assert isinstance(arg, Basic)
        self._args = args
        
    def __mathml__(self):
        """Returns a MathML expression representing the current object"""
        import xml.dom.minidom
        if self._mathml:
            return self._mathml
        dom = xml.dom.minidom.Document()
        x = dom.createElement('apply')
        x_1 = dom.createElement(self.mathml_tag)
        x.appendChild(x_1)
        for arg in self._args:
            x.appendChild( arg.__mathml__() )
        self._mathml = x
        return self._mathml
    
    def tryexpand(self, a):
        if isinstance(a, (Mul, Pow, Function)):
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
            if isinstance(x,(Real, Rational)):
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
            #checks wether all arguments are commutative
            if not x.is_commutative:
                return False
        return True

    def match(self, pattern, syms=None, exclude = "None"):
        """
        Imagine that we are matching (3*x**2).match(a*x,[a])

        The optional parameter "exclude" is a list of variables ("x") that we
        don't want to occur in the "a". This way, the user can specify, what
        kind of matches he is interested in. The default behavior (if the user
        doesn't specify this argument, i.e. exclude = "None") is:

            1)if there is just one more variable besides the matched ones, for
            example "x", then exclude = x.

            2)otherwise exclude=None, which means, that 
            Rational(2).match(a*x,[a],exlude = None) returns {a: 2/x}

        """
        #print self, pattern, syms
        from sympy.core.symbol import Symbol
        from sympy.core.numbers import Constant
        if syms == None:
            syms = pattern.atoms(type=Symbol)
            #print syms
        if exclude == "None":
            #extract the rest of the symbols in pattern (excluding syms)
            exclude = list(set(pattern.atoms(type=Symbol))-set(syms))
            if exclude == []:
                exclude = None
        if len(syms) == 1:
            if pattern.has(syms[0]):
                if pattern == syms[0]:
                    return {syms[0]: self}
            else:
                if self == pattern:
                    return {}
                else:
                    return None
        if isinstance(pattern, (Symbol, Constant)):
            try:
                return {syms[syms.index(pattern)]: self}
            except ValueError:
                pass
        if not isinstance(pattern, type(self)):
            return None
        ops = list(self[:])[:]
        pat = list(pattern[:])[:]
        global_wildcard = None
        for p in pat:
            if p in syms:
                global_wildcard = p
                break
        if global_wildcard:
            pat.remove(global_wildcard)
        #print ops,pat
        r2 = dict()
        for p in pat:
            for o in ops:
                if isinstance(o,Pair):
                    r = o.match(p,syms, exclude = exclude)
                else:
                    r = o.match(p,syms)
                #print o,p,syms,"->",r
                if r!= None:
                    ops.remove(o)
                    break
            if r == None:
                #print self, pattern, ops, pat
                if type(self) == type(pattern):
                    if isinstance(self,Mul):
                        r = (self/p).match(pattern/p, syms)
                        #print "  ",r
                        if exclude:
                            if r:
                                for x in r:
                                    if r[x].has(exclude[0]):
                                        return None
                                    if len(exclude) == 2:
                                        if r[x].has(exclude[1]):
                                            return None
                        return r
                return None
            r2.update(r)
        if global_wildcard:
            if len(ops) == 0:
                #return None
                if isinstance(self, Add):
                    rst = Rational(0)
                else:
                    rst = Rational(1)
            else:
                rst = type(self)(*ops)
            r2.update({global_wildcard: rst})
        #print "XX",r2,self,pat,global_wildcard
        return r2


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

        num = [] # items in numerator
        den = [] # items in denominator
        for x in a:
            if isinstance(x, Pow) and x.exp == -1:
                if isinstance(x.base, Pair):
                    den.append( "(" + str(x.base) + ")" )
                else:
                    den.append( str(x.base) )
            else:
                if isinstance(x, Pair):
                    num.append( "(" + str(x) + ")" )
                else:
                    num.append( str(x) )

        if len(num) == 0:
            num.append("1")

        snum = str.join('*', num)
        if len(den) == 0:
            return "%s%s" % (f, snum)
        else:
            sden = str.join('*', den)
            if len(den) > 1: sden = "(" + sden + ")"
            return "%s%s/%s" % (f, snum, sden)
    
    def __float__(self):
        a = 1
        for arg in self._args[:]:
            a *= float(arg)
        return a

    def __latex__(self):
        from symbol import Symbol
        from power import Pow
        f = []
        a = self._args
        if isinstance(a[0],Rational):
            if a[0].isminusone():
                f = ["-"]
                a = self._args[1:]
            elif a[0].isone():
                f = []
                a = self._args[1:]
        multsymb = "\cdot"
        for x in a:
            xs = x.__latex__()
            if isinstance(x,Pair):
                f += ["(%s)" % xs]
            # Insert explicit multiplication sign in some cases
            elif isinstance(x, Symbol) and "mathrm" in xs or \
                 isinstance(x, Pow) and "mathrm" in x.base.__latex__():
                f += [multsymb] + [xs] + [multsymb]
            else:
                f += [xs]
        # Remove extra multiplication signs
        for i in range(len(f)-1):
            if f[i] == f[i+1] == multsymb:
                f[i] = ""
        f = [x for x in f if x]
        if f[0] == multsymb: f = f[1:]
        if f[-1] == multsymb: f = f[:-1]
        return " ".join(f)

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

        if z1 is not None or z2 is not None:
            if z1 is not None and z2 is not None:
                #sanity check
                assert z1==z2
            if z1 is not None:
                return z1, True
            if z2 is not None:
                return z2, True
        
        # if it contains infinity
        if (x.atoms(type=Infinity) != []) or y.atoms(type=Infinity) != []:
            return x, False

        if isinstance(x,(Real, Rational)) and isinstance(y, (Real, Rational)):
            return x*y, True
        xbase,xexp = Mul.get_baseandexp(x)
        ybase,yexp = Mul.get_baseandexp(y)
        if xbase == ybase:
            #this whole "if" is to correctly cooperate with Pow.eval()
            #so we don't get infinite recursion. It's not elegant, but it
            #works.
            if not xbase.is_commutative:
                e = Add(xexp,yexp)
                if e != 0:
                    return Pow(xbase,e,evaluate=False), True
            return Pow(xbase,Add(xexp,yexp)), True
        elif xexp == 1 or yexp == 1 or xexp == -1 or yexp == -1:
            return x, False
        # hack to fix polynomial expansion
        elif xbase.is_number and y.base.is_number or \
             isinstance(x.base, Function) and isinstance(y.base, Function):
            if xexp == yexp:
                expr = Mul(xbase, ybase)
                # Don't evaluate Pow to avoid infinite recursion!
                # return Pow(expr, xexp, evaluate=True), True
                return Pow(expr, xexp, evaluate=False), True
            elif xexp == -yexp:
                den = Pow(ybase, -1)
                expr = Mul(xbase, den)
                return Pow(expr, xexp, evaluate=True), True
            return x, False
        else:
            return x, False
            
    def eval(self):
        "Flatten, put all Rationals in the front, sort arguments"

        def _mul_c(exp,x):
            e = []
            for i,y in enumerate(exp):
                z,ok = self.try_to_coerce(y,x)
                if isinstance(z, (Rational, Real)) and i != 0:
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
        if n == 0: 
            if self.atoms(type=Infinity) == []:
                return Rational(0)
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
        if n != 1: 
            a = [n]+a
        if len(a) > 1:
            #construct self again, but non-evaluated this time
            return Mul(evaluate=False, *a)
        elif len(a) == 1:
            return a[0]
        else:
            return Rational(1)
            
    def evalf(self):
        a, b = self.getab()
        if a.is_number and b.is_number:
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
        """Series expansion for Mul.

        It needs to handle this correctly:

        (e^x-1)/x  -> 1+x/2+x^2/6

        the first term is (e^x-1)/x evaluated at x=0, normally, we would use
        a limit. But in a CAS, the limit is computed using series, so we must
        calculate series differently - using the bottom up approach:
        first expand x, then e^x, then e^x-1, and finally (e^x-1)/x.
        """

        a,b=self.getab()
        try:
            x=a.series(sym,n)
        except pole_error:
            y=b.series(sym,n)
            a0 = y.subs(sym,0)
            if a0==0 and a.is_bounded:
                return y
            #we cannot expand x*y
            raise
        try:
            y=b.series(sym,n)
        except pole_error:
            #we are not able to expand b, 
            #but if a goes to 0 and b is bounded, 
            #the result is just a*const, so we just return a
            a0 = x.subs(sym,0)
            if a0==0 and b.is_bounded:
                return x
            #we cannot expand x*y
            raise
        #print x,y,x*y,(x*y).expand()
        return (x*y).expand()
        
    def expand(self):
        a,b = self.getab()
        a = self.tryexpand(a)
        b = self.tryexpand(b)
        if isinstance(a, Add):
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

    def combine(self):
        from functions import exp, log
        a,b = self.getab()
        a = a.combine()
        b = b.combine()
        if isinstance(a, exp) and isinstance(b, exp):
            return exp(a[0]+b[0])
        if isinstance(a, log):
            return log(a[0]**b)
        if isinstance(b, log):
            return log(b[0]**a)
        if isinstance(a, exp):
            for x in b:
                if isinstance(x, exp):
                    return exp(a[0]+x[0])*b/x
        return a*b
        
    def subs(self,old,new):
        a,b = self.getab()
        e = a.subs(old,new)*b.subs(old,new)
        if isinstance(e, Basic):
            return e
        else:
            return e

    def __pretty__(self):
        a = [] # items in the numerator
        b = [] # items that are in the denominator (if any)

        # Gather terms for numerator/denominator
        for item in self._args:
            if isinstance(item, Pow) and item.exp == -1:
                b.append(item.base)
            elif isinstance(item, Rational):
                if item.p != 1:
                    a.append( Rational(item.p) )
                if item.q != 1:
                    b.append( Rational(item.q) )
            else:
                a.append(item)

        # Convert to pretty forms
        for i in xrange(0, len(a)):
            if isinstance(a[i], Add) and len(a) > 1:
                a[i] = prettyForm(*a[i].__pretty__().parens())
            else:
                a[i] = a[i].__pretty__()

        for i in xrange(0, len(b)):
            if isinstance(b[i], Add) and len(b) > 1:
                b[i] = prettyForm(*b[i].__pretty__().parens())
            else:
                b[i] = b[i].__pretty__()

        # Construct a pretty form
        if len(b) == 0:
            return prettyForm.__mul__(*a)
        else:
            if len(a) == 0:
                a.append( Rational(1).__pretty__() )
            return prettyForm.__mul__(*a) / prettyForm.__mul__(*b)

class Add(Pair):
    """
    Usage
    =====
        This class represent's the addition of two elements. so whenever you call '+', an 
        instance of this class is created. 
        
    Notes
    =====
        When an instance of this class is created, the method .eval()
        is called and will perform some inexpensive simplifications. 
        
        In some cases, the eval() method will return an object that is
        not an instance of the class Add. So, for example, if x is a
        Symbol, (x+x) will create a class Add with arguments (x,x),
        that will be evaluated via the .eval() method, and this method
        will return a class Mul with arguments (2,x). That is how
        x+x --> 2*x is done.
        
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
    
    
    def __float__(self):
        a = 0
        for arg in self._args[:]:
            a += float(arg)
        return a  
    
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
    
    def __latex__(self):
        f = "%s" % self[0].__latex__()
        for i in range(1,len(self[:])):
            num_part = _extract_numeric(self[i])[0]
            if num_part < 0:
              f += "%s" % self[i].__latex__()
            else:
              f += "+%s" % self[i].__latex__()
        return f    

    def getab(self):
        """Pretend that self = a+b and return a,b
        
        in general, self=a+b+c+d+..., but in many algorithms, we 
        want to have just 2 arguments to add. Use this function to 
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
    

        
#    @debug
    def eval(self):
        """
        Usage
        =====
            This method is called automatically when an instance of
            this class is created, and will always return a new object. 
            
        Notes
        =====
            Currently, what this method does is: 
            
                - Flatten all arguments, i.e, substitute all instances
                  of Add by their arguments for example:
                  self.flatten(1, 1+x ) --> [1,1,x]
                
                - Add its arguments being aware of some identities,
                  like that x+x -> 2*x and that numbers can be added
                  without restrictions using their own __add__ method.
                
                - Sort
        See also
        ========
            L{Mul.eval}, L{Pow.eval}
            
        TODO
        ====
            - Perform a complexity analysis.
            - Probably optimizations can be done (algorithmic optimizations).
            - The problem is with addeval() and muleval(), that some classes
              implement, like Pauli and non-commutativity of some objects.
        """

#        @debug
        def _add(exp,x):
            an, a = _extract_numeric(x)
            e = []
            ok = False
            for y in exp:
            # try to put all numeric parts together
                bn, b = _extract_numeric(y)
                if (not ok) and a == b:
                    if isinstance(a, Infinity) or isinstance(b, Infinity):
                        if x != y:
                            # case oo - oo
                            raise ArithmeticError("Cannot compute oo-oo")
                    _m = Mul(an+bn, a)
                    if _m != 0:
                        e.append(_m)
                    ok = True
                else:
                    #print "I",x,y
                    z1 = x.addeval(y,x)
                    z2 = y.addeval(y,x)
                    #print "R",z1,z2,e,y,ok

                    if z1 is not None or z2 is not None:
                        if z1 is not None and z2 is not None:
                            #sanity check, only when z1 and z2 are not O
                            from symbol import O
                            if not isinstance(z1,O):
                                assert z1 == z2
                        if not ok:
                            if z1 is not None:
                                e.append(z1)
                                ok = True
                            elif z2 is not None:
                                e.append(z2)
                                ok = True
                    else:
                        e.append(y)
            if not ok: e.append(x)
            return e

        def _add_Number(a,b):
            """Adds two Real or Rational Numbers"""
            if isinstance(a,(Real, Rational)):
                return a + b
            else:
                raise ArgumentError
        
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

    def combine(self):
        r = Rational(0)
        for x in self:
            r+=x.combine()
        return r

    def removeO(self):
        return self.removeO()

    def removeO(self):
        """Removes the O(...) from the expression.
        
        Example:

        assert (2+O(x)) != 2
        assert (2+O(x)).removeO() == 2
        assert (2+x+O(x**2)).removeO() == x+2
        """
        from symbol import O
        for x in self:
            if isinstance(x, O):
                a = list(self[:])
                a.remove(x)
                return Add(*a)
        return self
        raise "the O not found, cannot remove"

    
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
            #implement the class O
            #print a.series(sym,n)
            #print b.series(sym,n)
            return (a.series(sym,n)+b.series(sym,n))
        
    def __pretty__(self):
        #return prettyForm.__add__(*[arg.__pretty__() for arg in self._args])
        from numbers import Number
        from stringPict import stringPict
        pforms = []
        for ind in xrange(0, len(self._args)):
            x = self._args[ind]
            if isinstance(x, Mul) and isinstance(x[0], Number) \
              and x[0] < 0:
                pform1 = (-x).__pretty__()
                if ind == 0:
                    if '\r' in str(pform1):
                        pform2 = '-'
                    else:
                        pform2 = '- '
                else:
                    pform2 = ' - '
                pform = stringPict.next(pform2, pform1)
                pforms.append(prettyForm(pform[0], pform[1], binding=prettyForm.NEG))
            else:
                pforms.append(x.__pretty__())
        return prettyForm.__add__(*pforms)

def _extract_numeric(x):
    """Returns the numeric and symbolic part of x.
    For example, 1*x -> (1,x)
    Works only with simple expressions. 
    """
    if isinstance(x, Mul) and isinstance(x._args[0], (Rational, Real)):
        return x.getab()
    else:
        return (Rational(1), x)
