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

class Matrix(object):

    def __init__(self,mat):
        self.lines=len(mat)
        self.cols=len(mat[0])
        self.mat=[]
        for j in range(self.lines):
            assert len(mat[j])==self.cols
            for i in range(self.cols):
                self.mat.append(Basic.sympify(mat[j][i]))

    def key2ij(self,key):
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.lines and j>=0 and j < self.cols):
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

    def __getattr__(self,name):
        if name == "T":
            #transposition
            r=zeronm(self.cols,self.lines)
            for i in range(self.lines):
                for j in range(self.cols):
                    r[j,i]=self[i,j]
            return r
        if name == "C":
            #conjugation
            r=zeronm(self.lines,self.cols)
            for i in range(r.lines):
                for j in range(r.cols):
                    r[i,j]=self[i,j].conjugate()
            return r
        if name == "H":
            #hermite conjugation
            return self.T.C
        if name == "D":
            #dirac conjugation
            return self.H * gamma(0)
        raise AttributeError("'%s' object has no attribute '%s'"%
                (self.__class__.__name__, name))

    def __getitem__(self,key):
        i,j=self.key2ij(key)
        return self.mat[i*self.cols+j]

    def __setitem__(self,key,value):
        i,j=self.key2ij(key)
        self.mat[i*self.cols+j] = value

    def hash(self):
        #if self.mhash: 
        #    return self.mhash.value
        self.mhash = hashing.mhash()
        self.mhash.addstr(str(type(self)))
        self.mhash.addint(self.lines)
        self.mhash.addint(self.cols)
        for x in self.mat:
            self.mhash.add(x.hash())
        return self.mhash.value

    def __rmul__(self,a):
        assert not isinstance(a,Matrix)
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=a*self[i,j]
        return r

    def expand(self):
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=self[i,j].expand()
        return r

    def subs(self,a,b):
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=self[i,j].subs(a,b)
        return r

    def __sub__(self,a):
        return self + (-a)

    def __mul__(self,a):
        if isinstance(a,Matrix):
            return self.multiply(a)
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=self[i,j]*a
        return r

    def __add__(self,a):
        return self.addeval(self,a)

    def __div__(self,a):
        return self * (Rational(1)/a)

    @staticmethod
    def addeval(x, y):
        if isinstance(x, Matrix) and isinstance(y, Matrix):
            return x.add(y)
        if isinstance(x, Matrix) and not isinstance(y, NCSymbol):
            assert x.lines == x.cols
            r=zeronm(x.lines,x.cols)
            for i in range(x.lines):
                for j in range(x.cols):
                    if i==j:
                        r[i,j]=x[i,j]+y
                    else:
                        r[i,j]=x[i,j]
            return r
        raise "unimplemented"

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
        if r.lines == 1 and r.cols ==1: 
            return r[0,0]
        return r

    def add(self,b):
        """Returns self+b """

        assert self.lines == b.lines
        assert self.cols == b.cols
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j] = self[i,j]+b[i,j]
        return r

    def __neg__(self):
        return -1*self

    def __eq__(self,a):
        if not isinstance(a, (Matrix, Basic)):
            a = Basic.sympify(a)
        return self.hash() == a.hash()

    def __str__(self):
        return self.print_sympy()

    def print_sympy(self):
        s="";
        for i in range(self.lines):
            for j in range(self.cols):
                s+="%s "%repr(self[i,j]);
            s+="\n"
        return s

def zero(n):
    return zeronm(n,n)

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

#one2=Pauli(0)
#sigma1=Pauli(1)
#sigma2=Pauli(2)
#sigma3=Pauli(3)

#assert sigma1 == sigma1
#assert sigma1 != sigma2

#assert sigma1*sigma2 == I*sigma3
#assert sigma3*sigma1 == I*sigma2
#assert sigma2*sigma3 == I*sigma1

#assert sigma1*sigma1 == one2
#assert sigma2*sigma2 == one2
#assert sigma3*sigma3 == one2

#assert sigma1*2*sigma1 == 2*one2
#assert sigma1*sigma3*sigma1 == -sigma3

a=Matrix((
    (1, 2),
    (3, 1),
    (0, 6),
    ))

b = Matrix ((
    (1, 2),
    (3, 0),
    ))

c= a*b
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

#gamma0=Dirac(0)
#gamma1=Dirac(1)
#gamma2=Dirac(2)
#gamma3=Dirac(3)
#gamma5=Dirac(5)

def sigma(i):
    if i==1:
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
    return Matrix(mat)

def gamma(mu,lower=False):
    if not mu in [0,1,2,3,5]:
        raise "Invalid Dirac index"
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
    m= Matrix(mat)
    if lower:
        g= one(4) 
        g[1,1] = -1
        g[2,2] = -1
        g[3,3] = -1
        m = g * m
    return m

gamma0=gamma(0)
gamma1=gamma(1)
gamma2=gamma(2)
gamma3=gamma(3)
gamma5=gamma(5)

assert I * gamma0 * gamma1 * gamma2 * gamma3 == gamma5

sigma1=sigma(1)
sigma2=sigma(2)
sigma3=sigma(3)

assert sigma1 == sigma1
assert sigma1 != sigma2

assert sigma1*sigma2 == I*sigma3
assert sigma3*sigma1 == I*sigma2
assert sigma2*sigma3 == I*sigma1

a=Symbol("a")
b=Symbol("b")
c=Symbol("c")

#print a*sigma1+b*sigma2+c*sigma3

E = Symbol("E")
m = Symbol("m")

def u(p,r):
    """ p = (p1, p2, p3); r = 0,1 """
    assert r in [1,2]
    p1,p2,p3 = p
    if r == 1:
        ksi = Matrix([ [1],[0] ])
    else:
        ksi = Matrix([ [0],[1] ])
    a = (sigma1*p1 + sigma2*p2 + sigma3*p3) / (E+m) * ksi
    if a ==0:
        a = zeronm(2,1)
    return (E+m).sqrt() * Matrix([ [ksi[0,0]], [ksi[1,0]], [a[0,0]], [a[1,0]] ])

def v(p,r):
    """ p = (p1, p2, p3); r = 0,1 """
    assert r in [1,2]
    p1,p2,p3 = p
    if r == 1:
        ksi = Matrix([ [1],[0] ])
    else:
        ksi = -Matrix([ [0],[1] ])
    a = (sigma1*p1 + sigma2*p2 + sigma3*p3) / (E+m) * ksi
    if a ==0:
        a = zeronm(2,1)
    return (E+m).sqrt() * Matrix([ [a[0,0]], [a[1,0]], [ksi[0,0]], [ksi[1,0]] ])

def pslash(p):
    p1,p2,p3 = p
    p0 = (m**2+p1**2+p2**2+p3**2).sqrt()
    return gamma0*p0-gamma1*p1-gamma2*p2-gamma3*p3

p = (a,b,c)

assert u(p, 1).D * u(p, 2) == 0
assert u(p, 2).D * u(p, 1) == 0

#e= v(p, 2).D * v(p, 2)
#print e.expand().subs(a, (E**2-m**2-b**2-c**2).sqrt()).expand()
#print
#e= u(p, 1) * u(p, 1).D+u(p, 2) * u(p, 2).D
#print e.expand()

#print
#f=pslash(p)+m

#print f.subs(a, (E**2-m**2-b**2-c**2).sqrt())

#print "-"*20
#print e.expand().subs(E, (m**2+a**2+b**2+c**2).sqrt())
#print f.subs(a, (E**2-m**2-b**2-c**2).sqrt()).subs(E,
#        (m**2+a**2+b**2+c**2).sqrt())

p1,p2,p3 =[Symbol(x) for x in ["p1","p2","p3"]]
pp1,pp2,pp3 =[Symbol(x) for x in ["pp1","pp2","pp3"]]
k1,k2,k3 =[Symbol(x) for x in ["k1","k2","k3"]]
kp1,kp2,kp3 =[Symbol(x) for x in ["kp1","kp2","kp3"]]

p = (p1,p2,p3)
pp = (pp1,pp2,pp3)

k = (k1,k2,k3)
kp = (kp1,kp2,kp3)

mu = Symbol("mu")

M0 = [ ( v(pp, 1).D * gamma(mu) * u(p, 1) ) * ( u(k, 1).D * gamma(mu,True) * \
        v(kp, 1) ) for mu in range(4)]
M = M0[0]+M0[1]+M0[2]+M0[3]
assert isinstance(M, Basic)

d=Symbol("d",True) #d=E+m

print M
print "-"*40
M = ((M.subs(E,d-m)).expand() * d**2 ).expand()
print "1/(E+m)**2 * ",M
#assert not M.has(d)
x= (M * M.conjugate())
print x
print x.expand()
