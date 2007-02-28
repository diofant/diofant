from sym import Rational,Function

class sin(Function):
    
    def getname(self):
        return "sin"
        
    def derivative(self):
        return cos(self.arg)
        
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,Rational) and self.arg.iszero():
            return Rational(0)
        return self.hold()

class cos(Function):
    
    def getname(self):
        return "cos"
        
    def derivative(self):
        return -sin(self.arg)
        
    def eval(self):
        if self.evaluated: return self
        if isinstance(self.arg,Rational) and self.arg.iszero():
            return Rational(1)
        return self.hold()

class tan(Function):
    
    def getname(self):
        return "tan"
        
    def derivative(self):
        return Rational(1)/cos(self.arg)**Rational(2)
        
    def eval(self):
        return (sin(self.arg)/cos(self.arg)).eval()
