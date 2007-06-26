"""Polynomial ideals with some arithmetic."""

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials.wrapper import div, groebner

class Ideal:
    """Describes a polynomial ideal over a field, with several variables.

    Try to avoid different variables and orders between the ideals, give
    them explicitly; the automatic handlers don't always act as expected.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> I = Ideal([x, y**2])
    >>> J = Ideal(x*y)
    >>> I == J
    False
    >>> I*J == I.intersect(J)
    False
    >>> I**2 == I*I
    True
    """
    def __init__(self, f=[Rational(0)], var=None, order='grevlex',
                 is_groebner=None):
        if not isinstance(f, list):
            f = [f]
        self.f = map(Basic.sympify, f)
        if var == None:
            var = []
            for p in self.f:
                var += filter(lambda x: not x in var, p.atoms(type=Symbol))
            var.sort(key=str)
        self.var = var
        self.order = order
        self.is_groebner = is_groebner
        for p in self.f:
            if not ispoly(p, self.var):
                raise PolynomialException('Non-polynomial as generator.')

    def __len__(self):
        return len(self.f)

    def __iter__(self):
        return self.f.__iter__()

    def __str__(self):
        return self.f.__str__()

    def __repr__(self):
        return 'Ideal(%s, %s, %s, %s)' % (repr(self.f), repr(self.var), 
               self.order, repr(self.is_groebner) )

    def __add__(self, other):
        """f is in I + J iff there are g in I, h in J such that f = g+h.
        """
        if not isinstance(other, Ideal):
            other = Ideal(other)
    
        var = self.var + filter(lambda x: not x in self.var, other.var)
        var.sort(key=str)

        if self.order == other.order:
            order = self.order
        else:
            order = None

        f = self.f + filter(lambda x: not x in self.f, other.f)
        
        return Ideal(f, var, order)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        """f is in I*J iff f is a sum of products of elements of I and J.
        """
        if not isinstance(other, Ideal):
            other = Ideal(other)

        var = self.var + filter(lambda x: not x in self.var, other.var)
        var.sort(key=str)

        if self.order == other.order:
            order = self.order
        else:
            order = None

        f = []
        for p in self.f:
            for q in other.f:
                f.append(p*q)

        return Ideal(f, var, order)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        """Repeated product of self.
        """
        other = Basic.sympify(other)
        if not (isinstance(other, Number)
                and other.is_integer
                and other >= 0):
            raise PolynomialException('Illegal power of ideal.')

        if other == 0:
            return Ideal(Rational(1), self.var, self.order, True)
        elif other == 1:
            return self

        # avoid repeated elements
        f = self.f[:]
        I = Ideal(Rational(0), self.var, self.order)

        for p in self.f:
            I += p*Ideal(f,self.var, self.order)**(other-1)
            f.remove(p)

        return I

    def __mod__(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other)
        other.groebner()

        f = []
        for p in self.f:
            f.append(div(p, other.f, other.var, other.order)[-1])
        f = filter(lambda x: x!=0, f)
        return Ideal(f, other.var, other.order)

    def __contains__(self, other):
        other = Basic.sympify(other)
        if not ispoly(other, self.var):
            return False
        self.groebner()
        rem = div(other, self.f, self.var, self.order)[-1]
        if rem == 0:
            return True
        else:
            return False

    def __eq__(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other, self.var, self.order)
        # Try to save Groebner base computations.       
        if self.var != other.var or self.order != other.order:
            if self.is_groebner:
                s = self
                o = Ideal(other.f, self.var, self.order)
                o.groebner()
            elif other.is_groebner:
                o = other
                s = Ideal(self.f, other.var, other.order)
                s.groebner()
            else:
                s = Ideal(self.f)
                o = Ideal(other.f)
                s.groebner()
                o.groebner()
            return s.f == o.f
        else:
            self.groebner()
            other.groebner()
            return self.f == other.f

    def groebner(self, reduced=True):
        if not self.is_groebner:
            self.f = groebner(self.f, self.var, self.order, None, reduced)
            self.is_groebner = True

    def intersect(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other)
        t = Symbol('t', dummy=True)
        I = t*self + (1-t)*other
        I.order = '1-el'
        I.groebner()
        f = filter(lambda p: not t in p.atoms(type=Symbol), I.f)
        if self.var == other.var:
            var = self.var
        else:
            var = None
        if self.order == other.order:
            order = self.order
        else:
            order = None
        return Ideal(f, var, order)
