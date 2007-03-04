from sym import Rational,Function

class sin(Function):
    """Return the sine of x (measured in radians)
    """
    
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
    """Return the cosine of x (measured in radians)
    """
    
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
    """Return the tangent of x (measured in radians)
    """
    
    def getname(self):
        return "tan"
        
    def derivative(self):
        return Rational(1) / (cos(self.arg)**2)
        
    def eval(self):
        return (sin(self.arg) / cos(self.arg)).eval()

    
