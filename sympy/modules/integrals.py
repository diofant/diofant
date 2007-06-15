from sympy import Basic, Symbol, Number, Mul, Pow, log, Add
from sympy import cos, sin
from sympy.core.stringPict import stringPict, prettyForm

class IntegralError(Exception):
    pass

class Integral(Basic):
    """
    Definite integral.

    Integral(f, (x,a,b) represents \int_a^b f(x) dx

    Usage
    =====

    print Integral(1/t, (t,1,x))  will print::
    
         int_{1}^{x} (t^(-1)) dt

    print Integral(1/t, (t,1,x)).doit() will print::
    
        log(x)

    print Integral(1/t, (t,1,x)).diff(x) will print::
    
        1/x

    Currently can only integrate very simple functions, like polynoms.
    You can however implement as many formulas as you want at the end
    of the primitive_function() function.

    """
    
    mathml_tag = 'int'

    def __init__(self, f, args):
        "int_a^b f(x)  dx"
        Basic.__init__(self)
        self.f=self.sympify(f)
        
        if isinstance(args, tuple):
            #case definite integral
            if len(args) != 3:
                print "Wrong number of arguments"
            self.x = self.sympify(args[0])
            self.a = self.sympify(args[1])
            self.b = self.sympify(args[2])
        else:
            assert isinstance(args, Basic)
            self.x = args
            self.a , self.b = None, None
            # todo self[:] ?
            
    def __latex__(self):
        if self.a is None:
            # if this is an indefinite integral
            return "\int %s\,d%s" % (self.f.__latex__(), self.x.__latex__())
        else:
            return "\int^%s_%s %s\,d%s" % (self.a.__latex__(), self.b.__latex__(), 
                                           self.f.__latex__(), self.x.__latex__() )
    
        
    def __pretty__(self):
        arg = prettyForm( *stringPict.next(self.f.__pretty__(), " d" , self.x.__pretty__()) )
        arg.baseline = 0

        # Create top and bottom parts based on definite/indefinite
        if self.a is not None:
            t = "/ %s" % self.b
            b = "/ %s" % self.a
        else:
            t = "/"
            b = "/"

        # Pad things so the integral sign is left-aligned properly
        width = max(len(t), len(b))
        t = prettyForm( t + (' ' * (width - len(t))) )
        b = prettyForm( b + (' ' * (width - len(b))) )

        # Create bar based on the height of the argument
        c = '|' + (' ' * (width -1))
        for x in xrange(1, arg.height()):
            c += '\r|' + (' ' * (width -1))

        # Construct the pretty form with the sign and the argument
        a = prettyForm(c)
        a = prettyForm(*a.below(b))
        a = prettyForm(*a.top(t))
        return prettyForm(*stringPict.right(a, " ", arg))
            
    def __mathml__(self):
        if self._mathml:
            return self._mathml
        import xml.dom.minidom
        dom = xml.dom.minidom.Document()
        x = dom.createElement("apply")
        x.appendChild(dom.createElement(self.mathml_tag))
        
        x_1 = dom.createElement('bvar')
        x_2 = dom.createElement('lowlimit')
        x_3 = dom.createElement('uplimit')
        
        x.appendChild(x_1)
        x.appendChild(x_2)
        x.appendChild(x_3)
        x.appendChild( self.f.__mathml__() )
        #TODO: add lowlimit, uplimit
        self._mathml = x
        
        return self._mathml
            
    def diff(self,sym):
        if sym==self.x:
            raise IntegralError("Cannot differentiate the integration variable")
        return (self.b.diff(sym)*self.f.subs(self.x,self.b)-\
            self.a.diff(sym)*self.f.subs(self.x,self.a))

    def __str__(self):
        if not isinstance(self.a, type(None)):
            # case definite integral
            return "int_{%r}^{%r} (%r) d%r"%(self.a,self.b,self.f,self.x)
        else:
            #case undefinite integral
            return "int(%r) d%r" % (self.f, self.x)

    def doit(self):
        """Try to do the integral."""
        F = self.primitive_function(self.f,self.x)
        if isinstance(self.a, type(None)):
            return F
        else:
            return (F.subs(self.x,self.b)-F.subs(self.x,self.a))

    @staticmethod
    def primitive_function(f,x):
        """Try to calculate a primitive function to "f(x)".
        
        Use heuristics and an integral table.
        """
        if isinstance(f,Mul):
            a,b = f.getab()
            if not a.has(x): return a*Integral.primitive_function(b,x)
            if not b.has(x): return b*Integral.primitive_function(a,x)
        if isinstance(f,Add):
            a,b = f.getab()
            return Integral.primitive_function(a,x)+Integral.primitive_function(b,x)
        if not f.has(x): return f*x
        if f==x: return x**2/2
        if isinstance(f,Pow):
            if isinstance(f.exp,Number):
                if x == f.base:
                    if f.exp==-1: return log(abs(x))
                    else: return x**(f.exp+1)/(f.exp+1)
                elif x in f.base and isinstance(f.base, Mul):
                    other = 1
                    for b in f.base:
                        if b != x: other *= b
                    other = other ** f.exp
                    
                    if f.exp==-1: return log(abs(x)) * other
                    else: return x**(f.exp+1)/(f.exp+1) * other

        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]
        integral_table = {
                a/(b*x+c): a/b * log(abs(b*x+c)),
                a*sin(b*x): -a/b * cos(b*x),
                a*cos(b*x): a/b * sin(b*x),
                log(x): x*log(x)-x
                }
        for k in integral_table:
            r = f.match(k, [a,b,c])
            if r != None:
                return integral_table[k].subs_dict(r)

        raise IntegralError("Don't know how to do this integral. :(")
        
    
def integrate(f, *args, **kargs):
    """
    Usage
    =====
    
      Indefinite integrals
      --------------------
    
      integrate(f, x) -> Returns the indefinite integral S{int} f(x) dx
      
      integrate(f, x, y) -> Return the indefinite double integral 
      S{int} S{int} f(x, y) dy dx 
      
      integrate(f, x, y, z, ...) -> Return the indefinite multiple integral
      (arbitrary number of variables) S{int} S{int} ... S{int} f(x, y, z, ...) dx ... dy dz
    
      
      Definite Integrals
      ------------------
    
      integrate(f, (x, a, b)) -> Returns the definite integral with integration 
      limits a, b
      
      integrate(f, (x, a, b), (y, c, d)) -> Returns the definite double integral
      
      
    Notes
    =====
    
      Currently only very simple integrals are computed.The general algorithm 
      for calculating integrals is described U{here<http://sympy.googlecode.com/svn/trunk/doc/issac98.pdf>}
      Someone just needs to implement it. :)
      
      Has an optional parameter evaluate, which can have value True or False. 
      If set to False, the integral will not be evaluated. Default is set to True.
      
   
    Further examples
    ================
      >>> from sympy import Symbol
      >>> x, y = Symbol('x'), Symbol('y')
      >>> integrate(2*x*y, (x,0,1), (y,-1,2))
      3/2
      >>> integrate(y, y)
      1/2*y**2
      >>> integrate(y*x, y)
      1/2*x*y**2
      >>> integrate(y*x, x, y)
      1/4*x**2*y**2
      >>> integrate(x*y**2 , (x,1,2), y)
      1/2*y**3
      >>> integrate(x , (x,1,2), evaluate=False)
      int_{1}^{2} (x) dx
      
    See also
    ========
    
      - L{limit<sympy.modules.limits.limit>}
       
      - External links
        - U{Riemman integral<http://planetmath.org/encyclopedia/RiemannIntegral.html>}
      
    """

    def _integrate_one_var(f, args, kargs):
        """Integrate over one variable"""
        if kargs.has_key('evaluate') and kargs['evaluate'] == False:
            return Integral(f, args)
        else:
            return Integral(f, args).doit()
    
    last_int = f
    for a in args:
        last_int= _integrate_one_var(last_int, a, kargs)
    return last_int
