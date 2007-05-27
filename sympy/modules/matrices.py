from sympy.core import Basic,exp,Symbol,Rational,I,Mul
from sympy.core import hashing

class NonSquareMatrixException(Exception):
    pass

class Matrix(object):

    def __init__(self, *args):
        """
        Matrix can be constructed with values or a rule.

        >>> from sympy import *
        >>> Matrix( (1,2+I), (3,4) ) #doctest:+NORMALIZE_WHITESPACE
        1 2+I
        3 4
        >>> Matrix(2, 2, lambda i,j: (i+1)*j ) #doctest:+NORMALIZE_WHITESPACE
        0 1
        0 2

        Note: in SymPy we count indices from 0. The rule however counts from 1.
        """
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.lines = args[0]
            self.cols = args[1]
            self.mat = []
            for i in range(self.lines):
                for j in range(self.cols):
                    self.mat.append(Basic.sympify(operation(i, j)))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            self.lines=args[0]
            self.cols=args[1]
            mat = args[2]
            self.mat=[]
            for j in range(self.lines):
                for i in range(self.cols):
                    self.mat.append(Basic.sympify(mat[j*self.cols+i]))
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            if not isinstance(mat[0], (list, tuple)):
                # make each element a singleton
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
            return Matrix(self.lines,self.cols, lambda i,j: self[j,i])
        if name == "C":
            #conjugation
            return Matrix(self.lines,self.cols,
                    lambda i,j: self[i,j].conjugate())
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
        if isinstance(key, slice):
            return self.mat[key]
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
        return hash(self.__str__() )

    def __rmul__(self,a):
        assert not isinstance(a,Matrix)
        r=zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=a*self[i,j]
        return r

    def expand(self):
        return Matrix(self.lines,self.cols, lambda i,j: self[i,j].expand())

    def combine(self):
        return Matrix(self.lines,self.cols, lambda i,j: self[i,j].combine())

    def subs(self,a,b):
        return Matrix(self.lines,self.cols, lambda i,j: self[i,j].subs(a,b))

    def __sub__(self,a):
        return self + (-a)

    def __mul__(self,a):
        if isinstance(a,Matrix):
            return self.multiply(a)
        return Matrix(self.lines,self.cols, lambda i,j: self[i,j]*a)

    def __pow__(self, num):
        if isinstance(num, int) or (isinstance(num, Rational) and num.isinteger()):
            if num < 0:
                a = self.inv() # A**-2 = (A**-1)**2
                num = -num
            else:
                a = self
            for i in range(1, num):
                a *= self
            return a
        raise NotImplementedError('Can only rise to the power of an integer for now')

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

        r = Matrix(self.lines,self.cols, lambda i,j: dotprod(self,b,i,j))
        if r.lines == 1 and r.cols ==1:
            return r[0,0]
        return r

    def add(self,b):
        """Returns self+b """

        assert self.lines == b.lines
        assert self.cols == b.cols
        return Matrix(self.lines,self.cols, lambda i,j:
                self[i,j]+b[i,j])

    def __neg__(self):
        return -1*self

    def __eq__(self,a):
        if not isinstance(a, (Matrix, Basic)):
            a = Basic.sympify(a)
        return self.hash() == a.hash()

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

    def __str__(self):
        s="";
        for i in range(self.lines):
            for j in range(self.cols):
                s+="%s "%repr(self[i,j]);
            s+="\n"
        return s

    def __mathml__(self):
        mml = ""
        for i in range(self.lines):
            mml += "<matrixrow>"
            for j in range(self.cols):
                mml += self[i,j].__mathml__()
            mml += "</matrixrow>"
        return "<matrix>" + mml + "</matrix>"

    def row(self, i, f):
        """Elementary row operation using functor"""
        for j in range(0, self.cols):
            self[i, j] = f(self[i, j], j)

    def col(self, j, f):
        """Elementary column operation using functor"""
        for i in range(0, self.lines):
            self[i, j] = f(self[i, j], j)

    def row_swap(self, i, j):
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]

    def col_swap(self, i, j):
        for k in range(0, self.lines):
            self[k, i], self[k, j] = self[k, j], self[k, i]

    def row_del(self, i):
        self.mat = self.mat[:i*self.cols] + self.mat[(i+1)*self.cols:]
        self.lines -= 1

    @property
    def is_square(self):
        return self.lines == self.cols

    def clone(self):
        return Matrix(self.lines, self.cols, lambda i, j: self[i, j])

    def det(self):
        """Compute matrix determinant using Bareis' fraction-free
           algorithm which is an extension of the well known Gaussian
           elimination method. This approach is best suited for dense
           symbolic matrices and will result in a determinant with
           minimal numer of fractions. It means that less term
           rewriting is needed on resulting formulae.

           TODO: Implement algorithm for sparse matrices (SFF).
        """

        if not self.is_square:
            raise NonSquareMatrixException()

        M, n = self.clone(), self.lines

        if n == 1:
            det = M[0, 0]
        elif n == 2:
            det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
        else:
            sign = 1 # track current sign in case of column swap

            for k in range(n-1):
                # look for a pivot in the current column
                # and assume det == 0 if none is found
                if M[k, k] == 0:
                    for i in range(k+1, n):
                        if M[i, k] != 0:
                            M.row_swap(i, k)
                            sign *= -1
                            break
                    else:
                        return Rational(0)

                # proceed with Bareis' fraction-free (FF)
                # form of Gaussian elimination algorithm
                for i in range(k+1, n):
                    for j in range(k+1, n):
                        D = M[k, k]*M[i, j] - M[i, k]*M[k, j]

                        if k > 0:
                            M[i, j] = D / M[k-1, k-1]
                        else:
                            M[i, j] = D

            det = sign * M[n-1, n-1]

        return det

def zero(n):
    return zeronm(n,n)

def zeronm(n,m):
    assert n>0
    assert m>0
    return Matrix(n,m, lambda i,j: 0)

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
