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
            for i in range(self.cols):
                self.mat.append(self.sympify(mat[j][i]))

    def key2ij(self,key):
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.lines and j>=0 and j < self.cols):
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

    def __getitem__(self,key):
        i,j=self.key2ij(key)
        return self.mat[i*self.cols+j]

    def __setitem__(self,key,value):
        i,j=self.key2ij(key)
        self.mat[i*self.cols+j] = value

    def hash(self):
        if self.mhash: 
            return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.lines)
        self.mhash.addint(self.cols)
        for x in self.mat:
            self.mhash.add(x.hash())
        return self.mhash.value

    @staticmethod
    def muleval(x, y):
        if isinstance(y, Matrix) and not isinstance(x, NCSymbol):
            r=zeronm(y.lines,y.cols)
            for i in range(y.lines):
                for j in range(y.cols):
                    r[i,j]=y[i,j]*x
            return r
        return None

    def multiply(self,b):
        """Returns self*b """

        def dotprod(a,b,i,j):
            assert a.cols == b.lines
            r=0
            for x in range(a.cols):
                r+=a[i,x]*b[x,j]
            return r

        r=zeronm(self.lines,b.cols)
        for i in range(self.lines):
            for j in range(b.cols):
                r[i,j] = dotprod(self,b,i,j)
        return r

    def print_sympy(self):
        s="";
        for i in range(self.lines):
            for j in range(self.cols):
                s+="%s "%repr(self[i,j]);
            s+="\n"
        return s

def zero(n):
    zeronm(n,m)

def zeronm(n,m):
    assert n>0
    assert m>0
    mat = ( [[0]*m]*n )
    return Matrix(mat)

def one(n):
    m = zero(n)
    for i in range(n):
        m[i,i]=1
    return m

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


one2=Pauli(0)
sigma1=Pauli(1)
sigma2=Pauli(2)
sigma3=Pauli(3)

assert sigma1 == sigma1
assert sigma1 != sigma2

assert sigma1*sigma2 == I*sigma3
assert sigma3*sigma1 == I*sigma2
assert sigma2*sigma3 == I*sigma1

assert sigma1*sigma1 == one2
assert sigma2*sigma2 == one2
assert sigma3*sigma3 == one2

assert sigma1*2*sigma1 == 2*one2
assert sigma1*sigma3*sigma1 == -sigma3

a=Matrix((
    (1, 2),
    (3, 1),
    (0, 6),
    ))

b = Matrix ((
    (1, 2),
    (3, 0),
    ))

c= a.multiply(b)
assert c[0,0]==7
assert c[0,1]==2
assert c[1,0]==6
assert c[1,1]==6
assert c[2,0]==18
assert c[2,1]==0

x = Symbol("x")

c = b * Symbol("x")
assert isinstance(c,Matrix)
assert c[0,0] == x
assert c[0,1] == 2*x
assert c[1,0] == 3*x
assert c[1,1] == 0

c = 5 * b
assert isinstance(c,Matrix)
assert c[0,0] == 5
assert c[0,1] == 2*5
assert c[1,0] == 3*5
assert c[1,1] == 0

gamma0=Dirac(0)
gamma1=Dirac(1)
gamma2=Dirac(2)
gamma3=Dirac(3)
gamma5=Dirac(5)

#print doit(gamma5 * gamma2)

#assert doit(I*Dirac(0)*Dirac(1)*Dirac(2)*Dirac(3)) == gamma5
print doit(I*Dirac(0)*Dirac(1)*Dirac(2)*Dirac(3)) 
