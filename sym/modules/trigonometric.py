from sym import rational,function

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
