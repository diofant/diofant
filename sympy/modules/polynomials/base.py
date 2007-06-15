"""Module providing the class Polynomial and low-level functions"""

from sympy import Add, Basic, Mul, Number, Pow, Rational, Real, Symbol
from sympy.modules.polynomials.common import *

coeff_rings = ['int', 'rat', 'real', 'cplx', 'sym']

class PolynomialException(Exception):
    pass

class Polynomial:
    """Polynomial representation in coefficient list form.

    This brings higher efficiency, but also more readable code, since
    there is no more need of two representation side-by-side that need
    to be in sync, in all the internal algorithms.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> f = Polynomial(x)
    >>> f
    Polynomial(x, [x], 'grevlex', 'int')
    >>> g = Polynomial(y+1)
    >>> f + g
    Polynomial(1+x+y, [x, y], 'grevlex', 'int')
    >>> f * g
    Polynomial(x+x*y, [x, y], 'grevlex', 'int')
    """
    def __init__(self, p, var=None, order='grevlex', coeff='int'):
        p = Basic.sympify(p)
        if var == None:
            var = p.atoms(type=Symbol)
            var.sort(key=str)
        self.var = var
        self.order = order
        self.p = coeff_list(p, self.var, self.order)
        # Try to find the lowest coefficient ring for this polynomial,
        # (they are sorted by inclusion, from left to right.)
        if coeff == None:
            coeff = 'int'
            for term in self.p:
                if not isinstance(term[0], Number):
                    coeff = 'sym'
                    break
                elif isinstance(term[0], Rational):
                    if term[0].is_integer:
                        coeff = max(coeff, 'int', key=coeff_rings.index)
                    else:
                        coeff = max(coeff, 'rat', key=coeff_rings.index)
                elif isinstance(term[0], Real):
                    coeff = max('real', coeff, key=coeff_rings.index)
                else:
                    coeff = 'cplx'
        self.coeff = coeff
        
    def __str__(self):
        return str(poly(self.p, self.var))

    def __repr__(self):
        return 'Polynomial(%s, %s, %s, %s)' % (repr(poly(self.p, self.var)),
               repr(self.var), repr(self.order), repr(self.coeff))

    def copy(self):
        p = Polynomial(0, self.var, self.order, self.coeff)
        # deep copy nested list
        p.p = []
        for term in self.p:
            p.p.append(term[:])
        return p

    def set(self, p=None, var=None, order=None, coeff=None):
        if p == None and (var != None or order != None):
            p = poly(self.p, self.var)
        if var != None:
            self.var = var
        if order != None:
            self.order = order
        if coeff != None:
            self.coeff = coeff
        if p != None:
            self.p = coeff_list(p, self.var, self.order)
                
    def __add__(self, other):
        if not isinstance(other, Polynomial):
            return self.__add__(Polynomial(other))

        if self.coeff == other.coeff:
            coeff = self.coeff
        else:
            coeff = max(self.coeff, other.coeff, key=coeff_rings.index)
            
        if self.var != other.var or self.order != other.order:
            if self.var != other.var:
                var = self.var + filter(lambda x: not x in self.var, other.var)
                var.sort(key=str)
            else:
                var = self.var
            if self.order == other.order:
                order = self.order
            else:
                order = 'grevlex'
            s = self.copy()
            s.set(var=var, order=order, coeff=coeff)
            o = other.copy()
            o.set(var=var, order=order, coeff=coeff)
            return s+o

        # Finally, the actual addition can begin!
        r = Polynomial(0, self.var, self.order, coeff)
        r.p = []
        # Merge the terms of self and other:
        i, j = 0, 0
        while i < len(self.p) and j < len(other.p):
            if (self.p[i][1:] == other.p[j][1:]):
                c = self.p[i][0]+other.p[j][0]
                if c != 0:
                    r.p.append([c] + self.p[i][1:])
                i += 1
                j += 1
            elif term_cmp(self.p[i], other.p[j], r.order) > 0:
                r.p.append(self.p[i])
                i += 1
            else:
                r.p.append(other.p[j])
                j += 1
        r.p += self.p[i:]
        r.p += other.p[j:]
        # Check if something was appended to the (empty) result.
        if len(r.p) == 0:
            r.set(p=0)
        # Check for remaining zero terms (coming directly from self or other?).
        if len(r.p) > 1:
            r.p = filter(lambda t:t[0] != 0, r.p)
        return r

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + (-1)*other

    def __rsub__(self, other):
        return other + (-1)*self

    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            return self.__mul__(Polynomial(other))

        if self.coeff == other.coeff:
            coeff = self.coeff
        else:
            coeff = max(self.coeff, other.coeff, key=coeff_rings.index)
            
        if self.var != other.var or self.order != other.order:
            if self.var != other.var:
                var = self.var + filter(lambda x: not x in self.var, other.var)
                var.sort(key=str)
            else:
                var = self.var
            if self.order == other.order:
                order = self.order
            else:
                order = 'grevlex'
            s = self.copy()
            s.set(var=var, order=order, coeff=coeff)
            o = other.copy()
            o.set(var=var, order=order, coeff=coeff)
            return s*o

        # Finally, the actual multiplication can begin!
        r = Polynomial(0, self.var, self.order, coeff)
        # Distribute the multiplication
        for self_term in self.p:
            temp = other.copy()
            for i in range(0, len(temp.p)):
                temp.p[i][0] *= self_term[0]
                for j in range(1, len(self_term)):
                    temp.p[i][j] += self_term[j]
            r += temp
        return r
        
    def __rmul__(self, other):
        return self.__mul__(other)

    def poly(self):
        return poly(self.p, self.var)

    def diff(self, v):
        if not v in self.var:
            return Polynomial(p=Rational(0), order=self.order)
        else:
            r = self.copy()
            r.p = []
            i = self.var.index(v) + 1
            for term in self.p: # iterate over non-changed list
                if term[i] > 0:
                    copy = term[:]
                    copy[0] *= copy[i]
                    copy[i] -= 1
                    r.p.append(copy)
            if len(r.p) == 0:
                r.var = []
                r.p = [[0]]
            return r


def coeff_list(p, var=None, order='grevlex'):
    """Return the list of coeffs and exponents.

    Currently, lexicographic ('lex'), graded lex ('grlex'), graded
    reverse lex ('grevlex') and 1-elimination ('1-el')orders are implemented.
    The list of variables determines the order of the variables.
    The coefficients are assumed to be non-zero real numbers, that is,
    you can divide by them.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> from sympy.modules.trigonometric import sin
    >>> coeff_list(x+sin(y)*x**2, [x])
    [[sin(y), 2], [1, 1]]

    >>> coeff_list(6*y**3+7*y*x**2)
    [[7, 2, 1], [6, 0, 3]]

    >>> coeff_list(6*y**3+7*y*x**2, [y,x])
    [[6, 3, 0], [7, 1, 2]]

    """

    p = Basic.sympify(p)
    p = p.expand()

    if isinstance(var, Symbol):
        var = [var]
    if var == None:
        var = p.atoms(type=Symbol)
        var.sort(Basic.cmphash)

    res = []

    if isinstance(p, Add):
        for a in p._args:
            res.append(*coeff_list(a, var, order))
    else:
        if not isinstance(p, Mul):
            p = Mul(Rational(1), p, evaluate=False)
        item = [Rational(1)] + len(var)*[Rational(0)]
        for factor in p._args:
            # check if any of the var appear
            if filter(lambda x:x in var, factor.atoms(type=Symbol)):
                if isinstance(factor, Pow) \
                   and (factor.base in var) \
                   and isinstance(factor.exp, Number) \
                   and factor.exp.is_integer \
                   and factor.exp > 0:
                    item[var.index(factor.base)+1] += factor.exp
                elif isinstance(factor, Symbol):
                    item[var.index(factor)+1] += 1
                else:
                    raise PolynomialException('Not a polynomial!')
            else: # the factor is relativly constant
                item[0] *= factor
        res = [item]

    # sort list according to monomial ordering
    if order == 'lex':
        res.sort(key=lambda x: x[1:], reverse=True)
    elif order == 'grlex':
        res.sort(key=lambda x: [sum(x[1:])] + x[1:], reverse=True)
    elif order == 'grevlex':
        res.sort(key=lambda x: [sum(x[1:])]
                 + reverse(map(lambda l:-l, x[1:])), reverse=True)
    elif order == '1-el':
        res.sort(key=lambda x: [x[1]] + [sum(x[2:])]
                 + reverse(map(lambda l:-l, x[2:])), reverse=True)
    else:
        raise PolynomialException(str(order) + 'is not an implemented order.')

    # unify monomials
    result = []
    for item in res:
        filt = filter(lambda x: x[1:] == item[1:], result)
        if filt:
            result[ result.index(filt[0]) ][0] += item[0]
        else:
            result.append(item)

    return result


def ispoly(p, var=None):
    """
    Usage:
      ispoly(p, var) -> Returns True if p is a polynomial in variable var.
                        Returns False otherwise.

    Notes:
        You can check wether it's a polynomial in several variables at
        once giving a tuple of symbols second argument
        (like ispoly(x**2 + y + 1, (x,y)) ).

    Examples:
        >>> from sympy import *
        >>> from sympy.modules.polynomials import *
        >>> x = Symbol('x')
        >>> ispoly(x**2+x+1, x)
        True
        >>> y = Symbol('y')
        >>> ispoly(x**2 + y + 1, (x,y)) #polynomial in variables x and y
        True
        >>> ispoly(x**2 + exp(y) + 1, (x,y))
        False

    See also:
       L{coeff_list}, L{coeff}

    """
    p = Basic.sympify(p)

    if var == None:
        var = p.atoms(type=Symbol)
    elif isinstance(var, Basic):
        # if the var argument is not a tuple or list
        var = [var] # so we can iterate over it

    if len(var) == 0:
        return True # constant is polynomial.
    elif len(var) > 1:
        return ispoly(p, var[0]) and ispoly(p, var[1:])

    if not var[0] in p.atoms(type=Symbol):
        return True # constant is polynomial.

    # Now we look for one variable, that is in the expression.
    if isinstance(p, Pow):
        if isinstance(p.exp, Number) \
           and p.exp.is_integer \
           and p.exp > 0:
            return ispoly(p.base, var[0])
        else:
            return False
    elif isinstance(p, (Add, Mul)):
        a, b = p.getab()
        return ispoly(a, var[0]) and ispoly(b, var[0])
    elif isinstance(p, Number):
        return True
    elif isinstance(p, Symbol):
        return True
    else:
        return False


def poly(p, var):
    """Returns a sympy polynomial from the representation "p" returned by
    coeff_list().
    """

    if isinstance(var, Symbol):
        var = [var]

    if len(p) == 0:
        return Rational(0)
    elif len(p[0]) != len(var) + 1:
        raise PolynomialException('Wrong number of variables given.')

    r = 0
    for item in p:
        c = item[0]
        for v in var:
            c *= v**item[var.index(v)+1]
        r += c
    return r
