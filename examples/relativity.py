import sys
sys.path.append(".")
sys.path.append("..")

from sympy import Basic,exp,Symbol,sin,Rational, Matrix, Function, \
    Derivative

def grad(f,X):
    a=[]
    for x in X:
        a.append( f.diff(x) )
    return a

def d(m,x):
    return grad(m[0,0],x)

class MT(object):
    def __init__(self,m):
        self.gdd=m
        self.guu=m.inv()

    def __str__(self):
        return "g_dd =\n" + str(self.gdd)

    def dd(self,i,j):
        return self.gdd[i,j]

    def uu(self,i,j):
        return self.guu[i,j]

class G(object):
    def __init__(self,g,x):
        self.g = g
        self.x = x

    def udd(self,i,k,l):
        g=self.g
        x=self.x
        r=0
        for m in [0,1,2,3]:
            r+=g.uu(i,m)/2 * (g.dd(m,k).diff(x[l])+g.dd(m,l).diff(x[k]) \
                    - g.dd(k,l).diff(x[m]))
        return r

class Riemann(object):
    def __init__(self,G,x):
        self.G = G
        self.x = x

    def uddd(self,rho,sigma,mu,nu):
        G=self.G
        x=self.x
        r=G.udd(rho,nu,sigma).diff(x[mu])-G.udd(rho,mu,sigma).diff(x[nu])
        for lam in [0,1,2,3]:
            r+=G.udd(sigma,mu,lam)*G.udd(lam,nu,sigma) \
                -G.udd(sigma,nu,lam)*G.udd(lam,mu,sigma)
        return r

class Ricci(object):
    def __init__(self,R,x):
        self.R = R
        self.x = x

    def dd(self,mu,nu):
        R=self.R
        x=self.x
        r=0
        for lam in [0,1,2,3]:
            r+=R.uddd(lam,mu,lam,nu)
        return r

class nu(Function):
    def getname(self):
        return "nu"

class lam(Function):
    def getname(self):
        return "lambda"

t=Symbol("t")
r=Symbol("r")
theta=Symbol("theta")
phi=Symbol("phi")

gdd=Matrix(( 
    (-exp(nu(r)),0,0,0), 
    (0, exp(lam(r)), 0, 0),
    (0, 0, r**2, 0),
    (0, 0, 0, r**2*sin(theta)**2)
    ))
#gdd=Matrix(( 
#    (1, 0, 0, 0), 
#    (0, 1, 0, 0),
#    (0, 0, r**2, 0),
#    (0, 0, 0, r**2*sin(theta)**2)
#    ))
#gdd=Matrix(( 
#    (1, 0, 0, 0), 
#    (0, 1, 0, 0),
#    (0, 0, 1, 0),
#    (0, 0, 0, r**2)
#    ))
g=MT(gdd)
X=(t,r,theta,phi)
Gamma=G(g,X)
Rmn=Ricci(Riemann(Gamma,X),X)

#print g
#print Gamma.udd(0,1,0)
#print Gamma.udd(0,0,1)
#print
#print Gamma.udd(1,0,0)
#print Gamma.udd(1,1,1)
#print Gamma.udd(1,2,2)
#print Gamma.udd(1,3,3)
#print
#print Gamma.udd(2,2,1)
#print Gamma.udd(2,1,2)
#print Gamma.udd(2,3,3)
#print
#print Gamma.udd(3,2,3)
#print Gamma.udd(3,3,2)
#print Gamma.udd(3,1,3)
#print Gamma.udd(3,3,1)
#print "-"*40
print Rmn.dd(0,0)
print "scalar curvature:"
print Rmn.dd(0,0)+Rmn.dd(1,1)+Rmn.dd(2,2)+Rmn.dd(3,3)

"""
#this is the input to eigenmath:

# This script calculates the Einstein tensor for a static spherically symmetric
# metric.
#
# Cf. \"A first course in general relativity,\" Bernard F. Schutz, p. 255.
#
# The book tells us exactly what the Einstein tensor components should be.
# If we get the right answer then we can be reasonably sure that the script is
# correct. Once that is known then we can use the functions defined here in
# other scripts.
#
# This is the line element for the metric (Equation 10.7)
#
#   2     2 Phi   2    2 Lambda   2    2        2
# ds  = -e      dt  + e         dr  + r  d Omega
#
# where
#
#  2        2    2         2      2            2
# r  d Omega  = r  (d theta  + sin  theta d phi )
#
# Note: Phi and Lambda are both functions of r.

# Given the line element we can write the metric tensor by inspection:

gdd = ((-exp(2 Phi(r)),                0,   0,                0),
       (             0, exp(2 Lambda(r)),   0,                0),
       (             0,                0, r^2,                0),
       (             0,                0,   0, r^2 sin(theta)^2))

# Note: \"dd\" stands for two \"down\" indices, \"uu\" stands for two \"up\" indices.

# X is our coordinate system. We need this for computing gradients.

X = (t,r,theta,phi)

# Step 1: Calculate guu.

guu = inv(gdd)

# Step 2: Calculate the connection coefficients. Cf. Gravitation, p. 210.
#
# Gamma    = 1/2 (g     + g     - g    )
#      abc         ab,c    ac,b    bc,a
#
# Note: The comma means gradient which increases the rank of gdd by 1.

gddd = d(gdd,X)

# Note: We transpose indices so they match up with Gamma, i.e., we put them in
# alphabetical order.

GAMDDD = 1/2 (gddd +                # indices are already in correct order
transpose(gddd,2,3) -               # transpose c and b
transpose(transpose(gddd,2,3),1,2)) # transpose c and a, then b and a

# Raise first index.
#
#      a      au
# Gamma    = g   Gamma
#       bc            ubc
#
# Note: Sum over index u means contraction.

GAMUDD = contract(outer(guu,GAMDDD),2,3)

# Step 3. Calculate the Riemann tensor. Cf. Gravitation, p. 219.
#
# a is alpha
# b is beta
# c is gamma
# d is delta
# u is mu
#
#  a           a            a            a        u          a        u
# R     = Gamma      - Gamma      + Gamma    Gamma    - Gamma    Gamma
#   bcd         bd,c         bc,d         uc       bd         ud       bc
#
# Do the gradient once and save in a temporary variable.

tmp1 = d(GAMUDD,X)

# The Gamma Gamma product is a rank 6 tensor with dim 4 per rank.
# That works out to 4 to the 6th or 4,096 elements.
# Of course, we'll do the outer product and contract over u just once and save
# the result in a second temporary variable.

tmp2 = contract(outer(GAMUDD,GAMUDD),2,4)

# Now put it all together. Do the transpositions so the indices get matched up
# with R on the left, i.e., put them in alphabetical order.

RUDDD = transpose(tmp1,3,4) -             # transpose d and c
  tmp1 +                                  # already in correct order
  transpose(tmp2,2,3) -                   # transpose c and b
  transpose(transpose(tmp2,2,3),3,4)      # transpose d and b, then d and c

# Step 4: Calculate the Ricci tensor. Cf. Gravitation, p. 343.
#
#        a
# R   = R
#  uv     uav
#
# Contract over \"a\" (1st and 3rd indices).

RDD = contract(RUDDD,1,3)

# Step 5: Calculate the Ricci scalar. Cf. Gravitation, p. 343.
#
#      uv
# R = g   R
#          vu  ...the book has uv, does it give the same result?
#              Yes because the metric tensor is symmetric so it's ok to
#              transpose.
#              I prefer vu because it looks like raising an index.

R = contract(contract(outer(guu,RDD),2,3),1,2)

# Step 6: Finally, calculate the Einstein tensor. Cf. Gravitation, p. 343.
#
# G   = R   - 1/2 g   R
#  uv    uv        uv

GDD = RDD - 1/2 gdd R

# Next we compare this result with Schutz' book. Schutz p. 255 gives the
# following Einstein tensor components (all other components are zero):
#
#        1                d
# G   = ----  exp(2 Phi) ---- [r (1 - exp(-2 Lambda))]
#  tt     2               dr
#        r
#
#          1                                         2
# G   = - ---- exp(2 Lambda) (1 - exp(-2 Lambda)) + --- Phi'
#  rr       2                                        r
#          r
#
#                 2                               2
# G            = r  exp(-2 Lambda) [Phi'' + (Phi')  + Phi'/r
#  theta theta
#
#                                                   - Phi' Lambda' - Lamda'/r]
#
#               2
# G        = sin  theta G
#  phi phi               theta theta

Gtt = 1/r^2 exp(2 Phi(r)) d(r (1 - exp(-2 Lambda(r))),r)

Grr = -1/r^2 exp(2 Lambda(r)) (1 - exp(-2 Lambda(r))) + 2/r d(Phi(r),r)

Gthetatheta = r^2 exp(-2 Lambda(r)) (
  d(d(Phi(r),r),r) +
  d(Phi(r),r)^2 +
  d(Phi(r),r) / r -
  d(Phi(r),r) d(Lambda(r),r) -
  d(Lambda(r),r) / r)

Gphiphi = sin(theta)^2 Gthetatheta

# Put together the expected tensor:

expect = ((Gtt,   0,           0,       0),
          (  0, Grr,           0,       0),
          (  0,   0, Gthetatheta,       0),
          (  0,   0,           0, Gphiphi))

# Check that GDD is correct.

check(GDD = expect)

# Display the non-zero components of GDD.

display(Gtt)
display(Grr)
display(Gthetatheta)
display(Gphiphi)
"""
