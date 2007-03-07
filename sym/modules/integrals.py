from sym import Basic, Symbol, Number, Mul, Pow, log, Add

class IntegralError(Exception):
    pass

class Integral(Basic):
    """
    Definite integral.

    Integral(a,b,f,x) represents \int_a^b f(x) dx

    Usage:

    print Integral(1,x, 1/t ,t)
    
        will print: int_{1}^{x} (t^(-1)) dt

    print Integral(1,x, 1/t ,t).doit()
    
        will print: log(x)

    print Integral(1,x, 1/t ,t).diff(x)
    
        will print: 1/x

    
    Currently can only integrate very simple functions, like polynoms.
    You can however implement as many formulas as you want at the end
    of the primitive_function() function.

    The general algorithm for calculating integrals is described here:

    http://sympy.googlecode.com/svn/trunk/doc/issac98.pdf

    Someone just needs to implement it. :)
    """

    def __init__(self, a, b, arg, x):
        "int_a^b arg  dx"
        Basic.__init__(self)
        self.arg=self.sympify(arg)
        self.a=self.sympify(a)
        self.b=self.sympify(b)
        assert isinstance(x, Symbol)
        self.x=x

    def eval(self):
        return Integral(self.a,self.b,self.arg,self.x,evaluated=True)

    def diff(self,sym):
        if sym==self.x:
            raise IntegralError("Cannot differentiate the integration variable")
        return (self.b.diff(sym)*self.arg.subs(self.x,self.b)-\
            self.a.diff(sym)*self.arg.subs(self.x,self.a))

    def __str__(self):
        return "int_{%r}^{%r} (%r) d%r"%(self.a,self.b,self.arg,self.x)

    def doit(self):
        """Try to do the integral."""
        F=self.primitive_function(self.arg,self.x)
        return (F.subs(self.x,self.b)-F.subs(self.x,self.a))

    def primitive_function(self,f,x):
        """Try to calculate a primitive function to "f(x)".
        
        Use heuristics.
        """
        if isinstance(f,Mul):
            a,b=f.getab()
            if not a.has(x): return a*self.primitive_function(b,x)
            if not b.has(x): return b*self.primitive_function(a,x)
        if isinstance(f,Add):
            a,b=f.getab()
            return self.primitive_function(a,x)+self.primitive_function(b,x)
        if not f.has(x): return f*x
        if f==x: return x**2/2
        if isinstance(f,Pow):
            if f.base==x and isinstance(f.exp,Number):
                if f.exp==-1: return log(x)
                else: return x**(f.exp+1)/(f.exp+1)

        #Implement any other formula here

        raise IntegralError("Don't know how to do this integral. :(")
