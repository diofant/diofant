import sys
sys.path.append(".")

from sympy.modules.matrices import sigma, gamma, zero, one, I, Matrix, \
                                    minkowski_tensor, eye, randMatrix, \
                                    jacobian
from sympy import Symbol

def test_multiplication():
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

def test_power():
    A = Matrix([[2,3],[4,5]])
    assert (A**5)[:] == [6140, 8097, 10796, 14237]
    A = Matrix([[2, 1, 3],[4,2, 4], [6,12, 1]])
    assert (A**3)[:] == [290, 262, 251, 448, 440, 368, 702, 954, 433]

def test_Pauli():
    #this and the following test are testing both Pauli and Dirac matrices
    #and also that the general Matrix class works correctly in a real world
    #situation
    sigma1=sigma(1)
    sigma2=sigma(2)
    sigma3=sigma(3)

    assert sigma1 == sigma1
    assert sigma1 != sigma2

    assert sigma1*sigma2 == I*sigma3
    assert sigma3*sigma1 == I*sigma2
    assert sigma2*sigma3 == I*sigma1

    assert sigma1*sigma1 == one(2)
    assert sigma2*sigma2 == one(2)
    assert sigma3*sigma3 == one(2)

    assert sigma1*2*sigma1 == 2*one(2)
    assert sigma1*sigma3*sigma1 == -sigma3

def test_Dirac():
    gamma0=gamma(0)
    gamma1=gamma(1)
    gamma2=gamma(2)
    gamma3=gamma(3)
    gamma5=gamma(5)

    assert gamma5 == I * gamma0 * gamma1 * gamma2 * gamma3
    assert gamma1 * gamma2 + gamma2 * gamma1 == zero(4)
    assert gamma0 * gamma0 == one(4) * minkowski_tensor[0,0]
    assert gamma2 * gamma2 != one(4) * minkowski_tensor[0,0]
    assert gamma2 * gamma2 == one(4) * minkowski_tensor[2,2]

    assert gamma(5,True) == \
        I*gamma(0,True)*gamma(1,True)*gamma(2,True)*gamma(3,True)

def test_creation():
    x = Symbol("x")
    a = Matrix([x, 0], [0, 0])
    m = a
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]
    b = Matrix(2,2, [x, 0, 0, 0])
    m = b
    assert m.cols == m.lines
    assert m.cols == 2
    assert m[:] == [x,0,0,0]

    assert a == b

def test_determinant():
    x, y = Symbol('x'), Symbol('y')

    assert Matrix([ [1] ]).det() == 1

    assert Matrix(( (-3,  2),
                    ( 8, -5) )).det() == -1

    assert Matrix(( (x,   1),
                    (y, 2*y) )).det() == 2*x*y-y

    assert Matrix(( (1, 1, 1),
                    (1, 2, 3),
                    (1, 3, 6) )).det() == 1

    assert Matrix(( ( 3, -2,  0, 5),
                    (-2,  1, -2, 2),
                    ( 0, -2,  5, 0),
                    ( 5,  0,  3, 4) )).det() == -289

    assert Matrix(( ( 1,  2,  3,  4),
                    ( 5,  6,  7,  8),
                    ( 9, 10, 11, 12),
                    (13, 14, 15, 16) )).det() == 0

    assert Matrix(( (3, 2, 0, 0, 0),
                    (0, 3, 2, 0, 0),
                    (0, 0, 3, 2, 0),
                    (0, 0, 0, 3, 2),
                    (2, 0, 0, 0, 3) )).det() == 275

    assert Matrix(( (1, 0,  1,  2, 12),
                    (2, 0,  1,  1,  4),
                    (2, 1,  1, -1,  3),
                    (3, 2, -1,  1,  8),
                    (1, 1,  1,  0,  6) )).det() == -55

    assert Matrix(( (-5,  2,  3,  4,  5),
                    ( 1, -4,  3,  4,  5),
                    ( 1,  2, -3,  4,  5),
                    ( 1,  2,  3, -2,  5),
                    ( 1,  2,  3,  4, -1) )).det() == 11664

    assert Matrix(( ( 2,  7, -1, 3, 2),
                    ( 0,  0,  1, 0, 1),
                    (-2,  0,  7, 0, 2),
                    (-3, -2,  4, 5, 3),
                    ( 1,  0,  0, 0, 1) )).det() == 123

def test_submatrix():
    m0 = eye(4)
    assert m0[0:3, 0:3] == eye(3)
    assert m0[2:4, 0:2] == zero(2)

    m1 = Matrix(3,3, lambda i,j: i+j)
    assert m1[0,:] == Matrix(1,3,(0,1,2))
    assert m1[1:3, 1] == Matrix(2,1,(2,3))
    
    m2 = Matrix([0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15])
    assert m2[:,-1] == Matrix(4,1,[3,7,11,15])
    assert m2[-2:,:] == Matrix([[8,9,10,11],[12,13,14,15]])

def test_submatrix_assignment():
    m = zero(4)
    m[2:4, 2:4] = eye(2)
    assert m == Matrix((0,0,0,0), 
                        (0,0,0,0),
                        (0,0,1,0),
                        (0,0,0,1))
    m[0:2, 0:2] = eye(2)
    assert m == eye(4)
    m[:,0] = Matrix(4,1,(1,2,3,4))
    assert m == Matrix((1,0,0,0), 
                        (2,1,0,0),
                        (3,0,1,0),
                        (4,0,0,1))
    m[:,:] = zero(4)
    assert m == zero(4)
    m[:,:] = ((1,2,3,4),(5,6,7,8),(9, 10, 11, 12),(13,14,15,16))
    assert m == Matrix(((1,2,3,4),
                        (5,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))
    m[0:2, 0] = [0,0]
    assert m == Matrix(((0,2,3,4),
                        (0,6,7,8),
                        (9, 10, 11, 12),
                        (13,14,15,16)))

def test_reshape():
    m0 = eye(3)
    assert m0.reshape(1,9) == Matrix(1,9,(1,0,0,0,1,0,0,0,1))
    m1 = Matrix(3,4, lambda i,j: i+j)
    assert m1.reshape(4,3) == Matrix((0,1,2), (3,1,2), (3,4,2), (3,4,5))
    assert m1.reshape(2,6) == Matrix((0,1,2,3,1,2), (3,4,2,3,4,5))

def test_applyfunc():
    m0 = eye(3)
    assert m0.applyfunc(lambda x:2*x) == eye(3)*2
    assert m0.applyfunc(lambda x: 0 ) == zero(3)

def test_random():
    M = randMatrix(3,3)
    M = randMatrix(3,3,seed=3)
    M = randMatrix(3,4,0,150)

def test_LUdecomp():
    testmat = Matrix([[0,2,5,3],
                      [3,3,7,4],
                      [8,4,0,2],
                      [-2,6,3,4]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)
    
    testmat = Matrix([[6,-2,7,4],
                      [0,3,6,7],
                      [1,-2,7,4],
                      [-9,2,6,3]])
    L,U,p = testmat.LUdecomposition()
    assert L.is_lower()
    assert U.is_upper()
    assert (L*U).permuteBkwd(p)-testmat == zero(4)

def test_LUsolve():
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    x = Matrix(3,1,[3,7,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x
    A = Matrix([[0,-1,2],
                [5,10,7],
                [8,3,4]])
    x = Matrix(3,1,[-1,2,5])
    b = A*x
    soln = A.LUsolve(b)
    assert soln == x

def test_inverse():
    A = eye(4)
    assert A.inv() == eye(4)
    A = Matrix([[2,3,5],
                [3,6,2],
                [8,3,6]])
    Ainv = A.inv()
    assert A*Ainv == eye(3)
    Ainv2 = A.inv("LU")
    assert Ainv == Ainv2

def test_cross():
    v1 = Matrix(1,3,[1,2,3])
    v2 = Matrix(1,3,[3,4,5])
    assert v1.cross(v2) == Matrix(1,3,[-2,4,-2])

def test_cofactor():
    assert eye(3) == eye(3).cofactorMatrix()
    test = Matrix([[1,3,2],[2,6,3],[2,3,6]])
    assert test.cofactorMatrix() == Matrix([[27,-6,-6],[-12,2,3],[-3,1,0]])
    test = Matrix([[1,2,3],[4,5,6],[7,8,9]])
    assert test.cofactorMatrix() == Matrix([[-3,6,-3],[6,-12,6],[-3,6,-3]])

def test_jacobian():
    x = Symbol('x')
    y = Symbol('y')
    L = [x**2*y, 2*y**2 + x*y]
    syms = [x,y]
    assert jacobian(L, syms) == Matrix([[2*x*y, x**2],[y, 4*y+x]])

    L = [x, x**2*y**3]
    assert jacobian(L, syms) == Matrix([[1, 0], [2*x*y**3, x**2*3*y**2]])

def test_charpoly():
    x = Symbol('x')
    y = Symbol('y')
    eye3 = eye(3)
    assert eye3.charpoly(x) == (1-x)**3
    assert eye3.charpoly(y) == (1-y)**3

