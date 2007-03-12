from sympy.core import Basic,exp,Symbol,Rational,I,Mul,NCSymbol
from sympy.core import hashing

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
        self.mat[i*self.cols+j] = Basic.sympify(value)

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

    def __repr__(self):
        return str(self)

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
