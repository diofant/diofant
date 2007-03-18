import sys
sys.path.append(".")
sys.path.append("..")

from sympy import Basic,exp,Symbol,sin,Rational,I,Mul,NCSymbol, Matrix, \
    gamma, sigma

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

#implement this using the NCSymbols, in paulialgebra.py in modules
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


def get_re_im(a):
    x = a.subs(I,0)
    y = (a+(-x).expand()).subs(I,1)
    return x,y

#gamma^mu
gamma0=gamma(0)
gamma1=gamma(1)
gamma2=gamma(2)
gamma3=gamma(3)
gamma5=gamma(5)

#sigma_i
sigma1=sigma(1)
sigma2=sigma(2)
sigma3=sigma(3)

a=Symbol("a")
b=Symbol("b")
c=Symbol("c")

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
print "-"*40
x,y= get_re_im(M)
print x,y
e = x**2+y**2
print e
