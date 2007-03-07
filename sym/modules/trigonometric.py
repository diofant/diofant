from sym import Rational,Function
from sym import pi

class sin(Function):
    """Return the sine of x (measured in radians)
    """
    
    def getname(self):
        return "sin"
        
    def derivative(self):
        return cos(self.arg)
        
    def eval(self):
        if self.arg == 0:
            return Rational(0)
        if self.arg==pi or self.arg==2*pi:
            return Rational(0)
        return self

class cos(Function):
    """Return the cosine of x (measured in radians)
    """
    
    def getname(self):
        return "cos"
        
    def derivative(self):
        return -sin(self.arg)
    
    def eval(self):
        if self.arg == 0:
            return Rational(1)
        if self.arg==pi:
            return -Rational(1)
        if self.arg==2*pi:
            return Rational(1)
        return self

class tan(Function):
    """Return the tangent of x (measured in radians)
    """
    
    def getname(self):
        return "tan"
        
    def derivative(self):
        return Rational(1) / (cos(self.arg)**2)
        
    def eval(self):
        return sin(self.arg) / cos(self.arg)

    
