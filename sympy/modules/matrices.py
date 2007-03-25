from sympy.core import Basic,exp,Symbol,Rational,I,Mul,NCSymbol
from sympy.core import hashing

class Matrix(object):

    def __init__(self, *args):
        """
        Matrix can be constructed with values or a rule.
        >>> from sympy import *
        >>> Matrix( (1,2+I), (3,4) )  #doctest: +NORMALIZE_WHITESPACE
        1 2+I
        3 4
        >>> Matrix(2, 2, lambda i,j: i*j )  #doctest: NORMALIZE_WHITESPACE
        1 2
        2 4
        """
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.lines = args[0]
            self.cols = args[1]
            self.mat = []
            for i in range(self.lines):
                for j in range(self.cols):
                    self.mat.append(operation(i+1, j+1))
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            if not isinstance(mat[0], (list, tuple)):
                # make each elements singleton
                mat = [ [element] for element in mat ]
            self.lines=len(mat)
            self.cols=len(mat[0])
            self.mat=[]
            for j in range(self.lines):
                assert len(mat[j])==self.cols
                for i in range(self.cols):
                    self.mat.append(Basic.sympify(mat[j][i]))

    def key2ij(self,key):
        """Converts key=(4,6) to 4,6 and ensures the key is correct."""
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.lines and j>=0 and j < self.cols):
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

    def __getattr__(self,name):
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2+I
        3 4
        >>> m.T #doctest: +NORMALIZE_WHITESPACE
        1 3
        2+I 4
        >>> m.H #doctest: +NORMALIZE_WHITESPACE
        1 3
        2-I 4

        """
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
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2+I
        3 4
        >>> m[1,0]
        3
        >>> m.H[1,0]
        2-I

        """
        i,j=self.key2ij(key)
        return self.mat[i*self.cols+j]

    def __setitem__(self,key,value):
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2+I
        3 4
        >>> m[1,0]=9 
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2+I
        9 4

        """
        i,j=self.key2ij(key)
        self.mat[i*self.cols+j] = Basic.sympify(value)

    def hash(self):
        """Compute a hash every time, because the matrix elements
        could change."""
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
        return self.add(a)

    def __div__(self,a):
        return self * (Rational(1)/a)

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

    def inv(self):
        assert self.cols==self.lines
        m=zero(self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                x=self[i,j]
                if i==j: 
                    x=1/x
                else:
                    if x!=0:
                        raise NotImplementedError("Matrix inversion is \
                            currently only implemented for diagonal matrices")
                m[i,j] = x
        return m

    def print_sympy(self):
        s="";
        for i in range(self.lines):
            for j in range(self.cols):
                s+="%s "%repr(self[i,j]);
            s+="\n"
        return s

    @property
    def mathml(self):
        mml = ""
        for i in range(self.lines):
            mml += "<matrixrow>"
            for j in range(self.cols):
                mml += self[i,j].mathml
            mml += "</matrixrow>"
        return "<matrix>" + mml + "</matrix>"

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
    """Returns a Pauli matrix sigma_i. i=1,2,3 

    See also:

    http://en.wikipedia.org/wiki/Pauli_matrices

    """
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
    """Returns a Dirac gamma matrix gamma^mu in the standard 
    (Dirac) representation.

    If you want gamma_mu, use gamma(mu, True).
    
    We use a convention:

    gamma^5 = I * gamma^0 * gamma^1 * gamma^2 * gamma^3 
    gamma_5 = I * gamma_0 * gamma_1 * gamma_2 * gamma_3 = - gamma^5

    See also:

    http://en.wikipedia.org/wiki/Gamma_matrices

    """
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
        if mu in [1,2,3,5]:
            m = - m
    return m

#Minkowski tensor using the convention (+,-,-,-) used in the Quantum Field
#Theory
minkowski_tensor = Matrix( (
    (1,0,0,0),
    (0,-1,0,0),
    (0,0,-1,0),
    (0,0,0,-1)
    ))
