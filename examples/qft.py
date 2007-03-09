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

class DiracMul(Mul):
    """Implements simplifications for Pauli and Dirac matrices"""

    @staticmethod
    def try_to_coerce(x,xbase,xexp,  y):
        if isinstance(x, Pauli) and isinstance(y, Pauli):
            j=y.i
            k=x.i
            return Pauli(0)*delta(j,k) \
                +I*epsilon(j,k,1)*Pauli(1) \
                +I*epsilon(j,k,2)*Pauli(2) \
                +I*epsilon(j,k,3)*Pauli(3), True
        return Mul.try_to_coerce(x,xbase,xexp,y)

#    @staticmethod
#    def _domul(a, b):
#        return DiracMul(Basic.sympify(a), Basic.sympify(b))


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

#    @staticmethod
#    def _domul(a, b):
#        return DiracMul(Basic.sympify(a), Basic.sympify(b))

    def print_sympy(self):
        s="";
        for j in self.mat:
            for i in j:
                s+="%s "%repr(i);
            s+="\n"
        return s

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
            j=y.i
            k=x.i
            return Pauli(0)*delta(j,k) \
                +I*epsilon(j,k,1)*Pauli(1) \
                +I*epsilon(j,k,2)*Pauli(2) \
                +I*epsilon(j,k,3)*Pauli(3), True
        return None, False

    def print_sympy(self):
        if self.i == 0:
            return "one"
        return "sigma%d"%self.i

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
