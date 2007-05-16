"""
This module contain solvers for all kinds of equations,
both algebraic (solve) and differential (dsolve).
"""

from sympy import Basic, Symbol, Number, Mul, log, Add, \
        sin, cos, integrate, sqrt, exp, Rational
        
from sympy.core.functions import Derivative, diff
from sympy.modules.polynomials import collect
from sympy.modules.matrices import zeronm

class Equation(Basic):

    def __init__(self, lhs, rhs):
        Basic.__init__(self)

        self._args = [Basic.sympify(lhs),
                      Basic.sympify(rhs)]

    @property
    def lhs(self):
        return self._args[0]

    @property
    def rhs(self):
        return self._args[1]

    @staticmethod
    def addeval(x, y):
        if isinstance(x, Equation) and isinstance(y, Equation):
            return Equation(x.lhs + y.lhs, x.rhs + y.rhs)
        else:
            return None

    @staticmethod
    def muleval(x, y):
        if isinstance(x, Equation) and isinstance(y, Equation):
            return Equation(x.lhs * y.lhs, x.rhs * y.rhs)
        else:
            return None
            
    def __eq__(self, other):
        if isinstance(other, bool):
            if other == True:
                return self.lhs.hash() == self.rhs.hash()
            else:
                return self.lhs.hash() != self.rhs.hash()
        else:
            raise TypeError
        
    def __nonzero__(self):
        return self.lhs.hash() == self.rhs.hash()

def solve(eq, vars):
    """
    Solves any (supported) kind of equation (not differential).

    Examples
    ========
      >>> from sympy import Symbol
      >>> x, y = Symbol('x'), Symbol('y')
      >>> solve(2*x-3, [x])
      3/2

    """
    
    if isinstance(vars, Basic):
        vars = [vars]

    if isinstance(vars, Symbol) or len(vars) == 1:
        x = vars[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]

        r = eq.match(a*x + b, [a,b]) # linear equation
        if r and _wo(r,x): return solve_linear(r[a], r[b])

        r = eq.match(a*x**2 + c, [a,c]) # quadratic equation
        if r and _wo(r,x): return solve_quadratic(r[a], 0, r[c])

        r = eq.match(a*x**2 + b*x + c, [a,b,c]) # quadratic equation
        if r and _wo(r,x): return solve_quadratic(r[a], r[b], r[c])
        
        d = Symbol('d', dummy=True)        
        r = eq.match(a*x**3 + b*x**2 + c*x + d, [a,b,c,d])
        if r and _wo(r, x): return solve_cubic(r[a], r[b], r[c], r[d])
        
        r = eq.match(a*x**3 - b*x + c)
        if r and _wo(r, x): return solve_cubic(r[a], 0, r[b], r[c])
    else:
        # augmented matrix
        n, m = len(eq), len(vars)
        matrix = zeronm(n, m+1)
        
        index = {}
        
        for i in range(0, len(vars)):
            index[vars[i]] = i
        
        for i in range(0, n):
            poly = collect(Basic.sympify(eq[i]).expand(), vars)
        
            if poly is not None:
                content, tail = poly
            else:
                break # at this point we will need Groebner basis 
                      # and multivariate polynomial factorization

            for (sym, coeff) in content.iteritems():
                matrix[i, index[sym]] = coeff
                
            matrix[i, m] = -tail
        else:
            return solve_linear_system(matrix, vars)

    raise "Sorry, can't solve it (yet)."

def solve_linear_system(matrix, syms):
    i, m = 0, matrix.cols-1 # don't count augmentation

    while i < matrix.lines:
        if matrix [i, i] == 0:
            # there is no pivot in current column
            # so try to find one in other colums
            for k in range(i+1, m):
                if matrix[i, k] != 0:
                    break
            else:
                if matrix[i, m] != 0:
                    return {}   # no solutions
                else:
                    # zero row or was a linear combination of 
                    # other rows so now we can safely skip it
                    matrix.row_del(i)
                    continue
    
            # we want to change the order of colums so
            # the order of variables must also change
            syms[i], syms[k] = syms[k], syms[i]            
            matrix.col_swap(i, k)
            
        pivot = matrix [i, i] 

        # divide all elements in the current row by the pivot
        matrix.row(i, lambda x, _: x / pivot)
            
        for k in range(i+1, matrix.lines):
            if matrix[k, i] != 0:
                coeff = matrix[k, i]
                
                # subtract from the current row the row containing 
                # pivot and multiplied by extracted coefficient
                matrix.row(k, lambda x, j: x - matrix[i, j]*coeff)
            
        i += 1
       
    # if there weren't any problmes, augmented matrix is now
    # in row-echelon form so we can check how many solutions
    # there are and extract them using back substitution

    if len(syms) == matrix.lines:
        # this system is Cramer equivalent so there is
        # finite number of solutions, exactly len(syms)
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # make back substitution for variables
            for j in range(k+1, m):
                content -= matrix[k, j]*solutions[syms[j]]
    
            solutions[syms[k]] = content.expand()
            
            k -= 1
            
        return solutions
    elif len(syms) > matrix.lines:
        # this system will have infinite number of solutions
        # dependent on exactly len(syms) - i parameters
        k, solutions = i-1, {}

        while k >= 0:
            content = matrix[k, m]

            # make back substitution for variables
            for j in range(k+1, i):
                content -= matrix[k, j]*solutions[syms[j]]
    
            # make back substitution for parameters
            for j in range(i, m):
                content -= matrix[k, j]*syms[j]
                
            solutions[syms[k]] = content.expand()
            
            k -= 1
            
        return solutions
    else:
        return {}   # no solutions

def solve_linear(a, b):
    """Solve a*x + b == 0"""
    return -b/a

def solve_quadratic(a, b, c):
    """Solve the cuadratic a*x**2 + b*x + c == 0"""
    D = b**2-4*a*c
    if D == 0:
        return [-b/(2*a)]
    else:
        return [
                (-b+sqrt(D))/(2*a),
                (-b-sqrt(D))/(2*a)
               ]
def solve_cubic(a, b, c, d):
    """Solve the cubic a*x**3 + b*x**2 + c*x + d == 0
    
    arguments are supposed to be sympy objects (so no python float's, int's, etc.)
    
    Cardano's method: http://en.wikipedia.org/wiki/Cubic_equation#Cardano.27s_method
    """
    # we calculate the depressed cubic t**3 + p*t + q
    
    #normalize
    a_1 = b / a
    b_1 = c / a
    c_1 = c / a
    
    del a, b, c
    
    p = b_1 - (a_1**2)/3
    q = c_1 + (2*a_1**3 - 9*a_1*b_1)/27
    
    u_1 = ( (q/2) + sqrt((q**2)/4 + (p**3)/27) )**Rational(1,3)
    u_2 = ( (q/2) - sqrt((q**2)/4 + (p**3)/27) )**Rational(1,3)
    # todo: this irnores
    
    x_1 = p/(3*u_1) - u_1 - a_1/3
    x_2 = p/(3*u_2) - u_2 - a_1/3
    
    return (x_1, x_2)
    

def dsolve(eq, funcs):
    """
    Solves any (supported) kind of differential equation.

    Usage
    =====
        dsolve(f, y(x)) -> Solve a differential equation f for the function y

        
    Details
    =======
        @param f: ordinary differential equation
        
        @param y: indeterminate function of one variable
    
        - you can declare the derivative of an unknown function this way:
        >>> from sympy import *
        >>> x = Symbol('x') # x is the independent variable
        >>> f = Function(x) # f is a function of f
        >>> f_ = Derivative(f, x) # f_ will be the derivative of f with respect to x
        
        - This function just parses the equation "eq" and determines the type of
        differential equation, then it determines all the coefficients and then
        calls the particular solver, which just accepts the coefficients.
        
    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> f = Function(x)
        >>> dsolve(Derivative(Derivative(f,x),x)+9*f, f)
        C1*sin(3*x)+C2*cos(3*x)

        #this is probably returned on amd64
        sin(3*x)*C1+cos(3*x)*C2

    """

    #currently only solve for one function
    if isinstance(funcs, Basic) or len(funcs) == 1:
        if isinstance(funcs, (list, tuple)): # normalize args
            f = funcs[0]
        else:
            f = funcs
            
        x = f[0]
        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]

        r = eq.match(a*Derivative(f,x) + b, [a,b])
        if r and _wo(r,f): return solve_ODE_first_order(r[a], r[b], f, x)

        r = eq.match(a*Derivative(Derivative(f,x),x) + b*f, [a,b])
        if r and _wo(r,f): return solve_ODE_second_order(r[a], 0, r[b], f, x)

        #special equations, that we know how to solve
        t = x*exp(-f)
        tt = (a*diff(t, x, 2)/t).expand()
        r = eq.match(tt, [a])
        if r:
            #check, that we've rewritten the equation correctly:
            #assert ( r[a]*diff(t, x,2)/t ) == eq.subs(f, t)
            return solve_ODE_1(f, x)
        eq = (eq*exp(f)/exp(-f)).expand()
        r = eq.match(tt, [a])
        if r:
            #check, that we've rewritten the equation correctly:
            #assert ( diff(t, x,2)*r[a]/t ).expand() == eq
            return solve_ODE_1(f, x)

    raise NotImplementedError("Sorry, can't solve it (yet)")

def solve_ODE_first_order(a, b, f, x):
    """ a*f'(x)+b = 0 """
    return integrate(-b/a, x) + Symbol("C1")

def solve_ODE_second_order(a, b, c, f, x):
    """ a*f''(x) + b*f'(x) + c = 0 """
    #a very special case, for b=0 and a,c not depending on x:
    return Symbol("C1")*sin(sqrt(c/a)*x)+Symbol("C2")*cos(sqrt(c/a)*x)

def solve_ODE_1(f, x):
    """ (x*exp(-f(x)))'' = 0 """
    C1 = Symbol("C1")
    C2 = Symbol("C2")
    return -log(C1+C2/x)


def _wo(di, x):
    """Are all items in the dictionary "di" without "x"?"""
    for d in di:
        if di[d].has(x):
            return False
    return True
