from sym import Basic, Symbol

class Integral(Basic):

    def __init__(self, a, b, arg, x, evaluated=False):
        "int_a^b arg  dx"
        Basic.__init__(self,evaluated)
        self.arg=arg
        self.a=self.sympify(a)
        self.b=self.sympify(b)
        assert isinstance(x, Symbol)
        self.x=x

    def eval(self):
        Integral(self.a,self.b,self.arg.eval(),self.x,evaluated=True)

    def diff(self,sym):
        return (self.b.diff(sym)*self.arg.subs(self.x,self.b)-\
            self.a.diff(sym)*self.arg.subs(self.x,self.a)).eval()

        if sym==self.x:
            return self.arg
        raise "Sorry, unimplemented yet."

    def __str__(self):
        return "int_{%r}^{%r} (%r) d%r"%(self.a,self.b,self.arg,self.x)
