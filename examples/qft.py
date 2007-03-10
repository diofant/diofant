import sys
sys.path.append(".")
sys.path.append("..")

from sympy import Basic,exp,Symbol,sin,Rational,I,Mul,NCSymbol
from sympy import hashing

def delta(i,j):
    if i==j:
        return 1
    else:
        return 0

def epsilon(i,j,k):
    if (i,j,k) in [(1,2,3), (2,3,1), (3,1,2)]:
        return 1
    elif (i,j,k) in [(1,3,2), (3,2,1), (2,1,3)]:
        return -1
    else:
        return 0

class Matrix(NCSymbol):

    def __init__(self,mat):
        Basic.__init__(self)
        self.lines=len(mat)
        self.cols=len(mat[0])
        self.mat=[]
        for j in range(self.lines):
            assert len(mat[j])==self.cols
            a=[]
            for i in range(self.cols):
                x=mat[j][i]
                if isinstance(x,int):
                    x=Rational(x)
                assert isinstance(x,Basic)
                a.append(x)
            self.mat.append(a)

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        for line in self.mat:
            for x in line:
                self.mhash.add(x.hash())
        return self.mhash.value

    @staticmethod
    def muleval(x, y):
        print "DDDD"
        if isinstance(x, Matrix) and isinstance(y, Basic) \
            and not isinstance(y, NCSymbol):
                mat=[[1,2,1,1],[1,2,1,1],[1,2,1,1],[1,2,1,1]]
                return Matrix(mat)
        return None

    def multiply(self,b):
        """ return self*b """

        def dotprod(a,b,j,i):
            r=0
            for x in range(4):
                r+=a[j][x]*b[x][i]
            return r

        assert self.cols == b.lines
        r=[
                [0,0,0,0],
                [0,0,0,0],
                [0,0,0,0],
                [0,0,0,0]
            ]
        for j in range(4):
            for i in range(4):
                r[j][i] = dotprod(self.mat,b.mat,j,i)
        return Matrix(r)

    def print_sympy(self):
        s="";
        for j in self.mat:
            for i in j:
                s+="%s "%repr(i);
            s+="\n"
        return s

class Dirac(Matrix):

    def __init__(self,mu):
        if not mu in [0,1,2,3,5]:
            raise "Invalid Dirac index"
        self.mu=mu
        if mu == 0:
            mat = (
                    (1,0,0,0),
                    (0,1,0,0),
                    (0,0,-1,0),
                    (0,0,0,-1)
                    )
        elif mu == 1:
            mat = (
                    (0,0,0,1),
                    (0,0,1,0),
                    (0,-1,0,0),
                    (-1,0,0,0)
                    )
        elif mu == 2:
            mat = (
                    (0,0,0,-I),
                    (0,0,I,0),
                    (0,I,0,0),
                    (-I,0,0,0)
                    )
        elif mu == 3:
            mat = (
                    (0,0,1,0),
                    (0,0,0,-1),
                    (-1,0,0,0),
                    (0,1,0,0)
                    )
        elif mu == 5:
            mat = (
                    (0,0,1,0),
                    (0,0,0,1),
                    (1,0,0,0),
                    (0,1,0,0)
                    )

        Matrix.__init__(self, mat)

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.mu)
        return self.mhash.value

    @staticmethod
    def one():
        return Matrix( ( 
            (1,0,0,0),
            (0,1,0,0),
            (0,0,1,0),
            (0,0,0,1)
            ))

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Dirac) and isinstance(y, Dirac):
            mu=x.mu
            nu=y.mu
            if mu==nu:
                if mu in [0,5]:
                    return Rational(1)
                else:
                    return -Rational(1)
            if mu == 5:
                return I*Dirac(0)*Dirac(1)*Dirac(2)*Dirac(3)*Dirac(nu)
            if nu == 5:
                return I*Dirac(mu)*Dirac(0)*Dirac(1)*Dirac(2)*Dirac(3)
        return None

    def print_sympy(self):
        return "gamma%d"%self.mu

class Pauli(Matrix):

    def __init__(self,i):
        if i==0:
            mat=( (
                (1, 0),
                (0, 1)
                ) )
        elif i==1:
            mat=( (
                (0, 1),
                (1, 0)
                ) )
        elif i==2:
            mat=( (
                (0, -I),
                (I, 0)
                ) )
        elif i==3:
            mat=( (
                (1, 0),
                (0, -1)
                ) )
        else:
            raise "Invalid Pauli index"
        self.i=i
        Matrix.__init__(self, mat)

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Pauli) and isinstance(y, Pauli):
            j=x.i
            k=y.i
            if j == 0: return x
            if k == 0: return y
            return Pauli(0)*delta(j,k) \
                +I*epsilon(j,k,1)*Pauli(1) \
                +I*epsilon(j,k,2)*Pauli(2) \
                +I*epsilon(j,k,3)*Pauli(3)
        return None

    def print_sympy(self):
        if self.i == 0:
            return "one"
        return "sigma%d"%self.i

def doit(e):
    assert isinstance(e, Mul)
    i=0
    while not isinstance(e.args[i],Matrix):
        i+=1
    r=Dirac.one()
    for x in e.args[i:]:
        r = r.multiply(x)
    return Mul([Rational(1)]+e.args[:i])*r


one=Pauli(0)
sigma1=Pauli(1)
sigma2=Pauli(2)
sigma3=Pauli(3)

assert sigma1 == sigma1
assert sigma1 != sigma2

assert sigma1*sigma2 == I*sigma3
assert sigma3*sigma1 == I*sigma2
assert sigma2*sigma3 == I*sigma1

assert sigma1*sigma1 == one
assert sigma2*sigma2 == one
assert sigma3*sigma3 == one

assert sigma1*2*sigma1 == 2*one
assert sigma1*sigma3*sigma1 == -sigma3

gamma0=Dirac(0)
gamma1=Dirac(1)
gamma2=Dirac(2)
gamma3=Dirac(3)
gamma5=Dirac(5)

print gamma0 * gamma1
print gamma0 * gamma0
print gamma2 * gamma5
print gamma5 * gamma2

#print doit(gamma5 * gamma2)

print doit(I*Dirac(0)*Dirac(1)*Dirac(2)*Dirac(3))
