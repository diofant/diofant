import hashing
from basic import basic
from numbers import rational

class function(basic):
    def __init__(self,arg):
        basic.__init__(self)
        self.arg=arg
    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash=hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.arg.hash())
        return self.mhash.value
    def diff(self,sym):
        return (self.derivative()*self.arg.diff(sym)).eval()
    def subs(self,old,new):
        e=basic.subs(self,old,new)
        #if e==self:
        if e.isequal(self):
            return (type(self)(self.arg.subs(old,new))).eval()
        else:
            return e
    def __str__(self):
        f="%s(%s)"
        return f%(self.getname(),str(self.arg))

class exp(function):
    def getname(self):
        return "exp"
    def derivative(self):
        return exp(self.arg)
    def eval(self):
        if self.evaluated: return self
        arg=self.arg.eval()
        if isinstance(arg,rational) and arg.iszero():
            return rational(1)
        if isinstance(self.arg,ln):
            return self.arg.arg
        return exp(arg).hold()

class ln(function):
    def getname(self):
        return "ln"
    def derivative(self):
        return rational(1)/self.arg
    def eval(self):
        if self.evaluated: return self
        self.arg=self.arg.eval()
        if isinstance(self.arg,rational) and self.arg.isone():
            return rational(0)
        if isinstance(self.arg,exp):
            return self.arg.arg
        return self.hold()

class sin(function):
    def getname(self):
        return "sin"
    def derivative(self):
        return cos(self.arg)
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,rational) and self.arg.iszero():
            return rational(0)
        return self.hold()

class cos(function):
    def getname(self):
        return "cos"
    def derivative(self):
        return -sin(self.arg)
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,rational) and self.arg.iszero():
            return rational(1)
        return self.hold()

class tan(function):
    def getname(self):
        return "tan"
    def derivative(self):
        return rational(1)/cos(self.arg)**rational(2)
    def eval(self):
        return (sin(self.arg)/cos(self.arg)).eval()
