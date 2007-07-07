"""
This module is experimental. There are likely bugs,
and the interface is likely to change.
------------------------------------------------------------------

This module implements arbitrary-precision binary floating-point
and fixed-point arithmetic. Due to performing all internal
arithmetic on long Python integers instead of lists of
single-digit integers, the arithmetic implemented here is
typically 10-100 faster than decimal.py.

The two main features are the classes BinaryReal and
BinaryComplex. Precision is controlled via the global 'context'
object (similar to decimal.py's contexts).

Test:

from sympy import *
from sympy.modules.numerics import *

context.prec = 100   # bits
pif()
expf(1)
cosf(1)
sinf(1)
logf(2)
expf(logf(3))
powerf(2,'0.5')
expf(BinaryComplex((0,1)))
logf(-1)
abs(expf(BinaryComplex((0,1)) * 4.56))

"""


from sympy import Basic, Rational, Add, Mul, Pow, pi
import math
import decimal

isi = isinstance
Q = Rational

#------------------------------------------------------------------
#
# First some utilities
#


def bitcount(n):
    """Give position of highest set bit in an integer"""
    if n == 0: return 0
    if n < 0: n = -n
    # math.log gives a good estimate, and never overflows, but not
    # is not always exact. Subtract 2 to understimate, then
    # count remaining bits "manually"
    bc = max(0, int(math.log(n, 2)) - 2)
    n >>= bc
    while n:
        n >>= 1
        bc += 1
    return bc

def trailing_zeros(n):
    """Count trailing zero bits in an integer"""
    if n == 0: return 0
    if n < 0: n = -n
    B=64; M=(1<<B)-1; t = 0
    while not n & M:
        n >>= B; t += B
    while not n & 1:
        n >>= 1; t += 1
    return t

# These functions are used to ensure symmetric rounding
# toward zero for both positive and negative integers
# and to allow negative shifts
def rshift(x, n):
    if n == 0: return x
    if n < 0:  return lshift(x, -n)
    if x >= 0: return x >> n
    else:      return -((-x) >> n)

def lshift(x, n):
    if n == 0: return x
    if n < 0:  return rshift(x, -n)
    if x >= 0: return x << n
    else:      return -((-x) << n)

def div(x, y):
    if y < 0 and x > 0: return -(x // -y)
    if y > 0 and x < 0: return -(-x // y)
    return x // y


# Fixed-point arithmetic helper functions

def fixed_normalize(a, abits, target_bits):
    if abits == target_bits: return a
    if abits > target_bits:  return a >> (abits - target_bits)
    if abits < target_bits:  return a << (target_bits - abits)

def fixed_add(a, abits, b, bbits, target_bits):
    x = fixed_normalize(a, abits, target_bits)
    y = fixed_normalize(b, bbits, target_bits)
    return x+y

def fixed_mul(a, abits, b, bbits, target_bits):
    if a == 0 or b == 0: return 0
    return fixed_normalize(a*b, abits+bbits, target_bits)

def fixed_div(a, abits, b, bbits, target_bits):
    if b == 0: raise ZeroDivisionError
    return fixed_normalize((a << bbits) // b, abits, target_bits)

def decimals_to_bits(n):
    return int(n*math.log(10,2)+2)

def bits_to_decimals(n):
    return int((n-2)/math.log(10,2))

def isreal(x):
    return isi(x, (int, long, float, Q, BinaryReal))


class NumericalContext:
    stack = []
    def enter(self, **kwargs):
        self.stack.append(self.__dict__.copy())
        self.__dict__.update(kwargs)
    def revert(self):
        self.__dict__.update(self.stack.pop())

FLOAT = 0
FIXED = 1

context = NumericalContext()
context.prec = 64
context.prec_mode = FLOAT


class Numerical:

    def __init__(s, x, prec=None, prec_mode=None):
        if prec == None:
            prec = context.prec
        if prec_mode == None:
            prec_mode = context.prec_mode
        s.init(x, prec, prec_mode)

    def do_op(s, method, t):
        prec = context.prec
        prec_mode = context.prec_mode

        if   method == '+': f = s.add
        elif method == '*': f = s.mul
        elif method == '/': f = s.div
        elif method == '<<': f = s.lshift
        elif method == '>>': f = s.rshift

        try:
            return f(s, t, prec, prec_mode)

        except (TypeError, NotImplementedError):
            try:
                if isi(t, BinaryComplex):
                    s = s.complex()
                t = s.__class__(t, prec, prec_mode)
            except (TypeError, NotImplementedError):
                raise Exception, ":("
            return s.do_op(method, t)

    def __pos__(s):
        prec = context.prec
        prec_mode = context.prec_mode
        return s.__class__(s, prec, prec_mode)

    def __add__(s, t): return s.do_op('+', t)
    def __radd__(s, t): return s + t
    def __sub__(s, t): return s + (-t)
    def __rsub__(s, t): return (-s) + t
    def __mul__(s, t): return s.do_op('*', t)
    def __rmul__(s, t): return s * t
    def __div__(s, t): return s.do_op('/', t)
    def __lshift__(s, t): return s.do_op('<<', t)
    def __rshift__(s, t): return s.do_op('>>', t)


class BinaryReal(Numerical):

    def init(s, x, prec, prec_mode):
        s.man = s.exp = 0

        if prec_mode == FLOAT:
            if   isi(x, tuple): s.man, s.exp = x
            elif isi(x, BinaryReal): s.man, s.exp = x.man, x.exp
            elif isi(x, (int, long)): s.man = x; s.exp = 0
            elif isi(x, float):
                m, e = math.frexp(x)
                s.man = int(m * 2**53)
                s.exp = e - 53
            elif isi(x, Q):
                n = prec + bitcount(x.q)
                s.man = div(lshift(x.p, n), x.q)
                s.exp = -n
            elif isi(x, str):
                s.man, s.exp = BinaryReal(Q(x), prec).man_exp()
            else:
                raise TypeError
            # Normalize
            bc = bitcount(s.man)
            if bc > prec:
                s.man = rshift(s.man, bc-prec)
                s.exp += (bc-prec)
            tr = trailing_zeros(s.man)
            if tr:
                s.man = rshift(s.man, tr)
                s.exp += tr
            if s.man == 0:
                s.exp = 0

        elif prec_mode == FIXED:
            if isi(x, tuple):
                s.man, s.exp = x
                assert isi(s.man, (int, long))
            elif isi(x, (int, long)):
                s.man = lshift(x, prec);
                s.exp = -prec
            elif isi(x, BinaryReal):
                s.man = lshift(x.man, x.exp+prec)
                s.exp = -prec
            elif isi(x, float):
                # XXX: remove this constructor
                t = BinaryReal(BinaryReal(x), prec=prec, prec_mode=FIXED)
                s.man, s.exp = t.man_exp()
            elif isi(x, Q):
                s.man = div(lshift(x.p, prec), x.q)
                s.exp = -prec
            else:
                raise TypeError

    def decimal(s, n):
        """Represent as a decimal string with at most n digits"""
        prec_ = decimal.getcontext().prec
        decimal.getcontext().prec = n
        if s.exp >= 0:
            d = decimal.Decimal(s.man) * (1<<s.exp)
        else:
            d = decimal.Decimal(s.man) / (1<<-s.exp)
        a = str(d)
        decimal.getcontext().prec = prec_
        return a

    def __str__(s):
        if context.prec_mode == FLOAT:
            return s.decimal(int(context.prec/math.log(10,2)))
        if context.prec_mode == FIXED:
            return s.decimal(int(bitcount(s.man)/math.log(10,2)))

    __repr__ = __str__

    def man_exp(s):
        return s.man, s.exp

    def complex(s):
        return BinaryComplex((s, 0))

    def rational(s):
        return s.man * Q(2)**(s.exp)

    def __float__(s):
        try:
            return math.ldexp(s.man, s.exp)
        except OverflowError:
            n = bitcount(s.man) - 64
            m = s.man >> n
            return math.ldexp(m, s.exp + n)

    def __nonzero__(s):
        return s.man != 0

    def __cmp__(s, t):
        if t == 0:
            return cmp(s.man, 0)
        elif s.man == 0:
            return cmp(0, t)
        elif isi(t, BinaryReal):
            if s.man < 0 and t.man > 0: return -1
            if s.man > 0 and t.man < 0: return 1
            if s.exp == t.exp: return cmp(s.man, t.man)
            a = cmp((s - t).man, 0)
            return a
        else:
            return cmp(s, BinaryReal(t))

    def __abs__(s):
        if s.man < 0:
            return -s
        return s

    def __neg__(s):
        return BinaryReal((-s.man, s.exp))

    @staticmethod
    def add(s, t, prec, prec_mode):
        if prec_mode == FLOAT:
            if isi(t, BinaryReal):
                if t.exp > s.exp:
                    s, t = t, s
                m = t.man + (s.man << (s.exp - t.exp))
                return BinaryReal((m, t.exp), prec)
        if prec_mode == FIXED:
            if isi(t, BinaryReal):
                r = fixed_add(s.man, -s.exp, t.man, -t.exp, prec)
                return BinaryReal((r, -prec))
            if isi(t, (int, long)):
                r = fixed_add(s.man, -s.exp, t, 0, prec)
                return BinaryReal((r, -prec))
        raise NotImplementedError

    @staticmethod
    def mul(s, t, prec, prec_mode):
        if prec_mode == FLOAT:
            if isi(t, BinaryReal):
                return BinaryReal((s.man*t.man, s.exp+t.exp))
            elif isi(t, (int, long)):
                return BinaryReal((s.man*t, s.exp))
        if prec_mode == FIXED:
            if isi(t, BinaryReal):
                r = fixed_mul(s.man, -s.exp, t.man, -t.exp, prec)
                return BinaryReal((r, -prec))
            if isi(t, (int, long)):
                return BinaryReal((s.man*t, s.exp))
        raise NotImplementedError

    @staticmethod
    def div(s, t, prec, prec_mode):
        if t == 0:
            raise ZeroDivisionError
        if prec_mode == FLOAT:
            if isi(t, BinaryReal):
                # XXX: do this faster
                extra = prec + bitcount(t.man)
                man = div(lshift(s.man, extra), t.man)
                exp = s.exp - t.exp - extra
                return BinaryReal((man, exp))
            elif isi(t, (int, long)):
                extra = prec + bitcount(t)
                man = div(lshift(s.man, extra), t)
                exp = s.exp - extra
                return BinaryReal((man, exp))
        if prec_mode == FIXED:
            if isi(t, BinaryReal):
                a = fixed_div(s.man, -s.exp, t.man, -t.exp, prec)
                return BinaryReal((a, -prec))
            if isi(t, (int, long)):
                return BinaryReal((div(s.man, t), -prec))
        raise NotImplementedError

    @staticmethod
    def lshift(s, n, prec, prec_mode):
        if prec_mode == FLOAT:
            return BinaryReal((s.man, s.exp+int(n)))
        if prec_mode == FIXED:
            return BinaryReal((lshift(s.man,n), s.exp))

    @staticmethod
    def rshift(s, n, prec, prec_mode):
        return s.lshift(s, -n, prec, prec_mode)



#------------------------------------------------------------------
#
# BinaryComplex
#

class BinaryComplex(Numerical):

    def init(s, val, prec, prec_mode):
        if isi(val, (BinaryComplex, complex)):
            real = val.real
            imag = val.imag
        elif isi(val, tuple):
            real, imag = val
        else:
            real = val
            imag = 0
        s.real = BinaryReal(real)
        s.imag = BinaryReal(imag)

    def __str__(s):
        return "%s + %s*I" % (s.real, s.imag)

    __repr__ = __str__

    def __complex__(s):
        return complex(float(s.real), float(s.imag))

    def conjugate(s):
        return BinaryComplex(s.real, -s.imag)

    def __eq__(s, t):
        t = BinaryComplex(t)
        return s.real == t.real and s.imag == t.imag

    def __nonzero__(s):
        return s.real != 0 or s.imag != 0

    def __neg__(s):
        return BinaryComplex((-s.real, -s.imag))

    def __abs__(s):
        return powerf(s.real*s.real + s.imag*s.imag, '0.5')

    @staticmethod
    def add(s, t, prec, prec_mode):
        if isi(t, BinaryComplex):
            return BinaryComplex((s.real+t.real, s.imag+t.imag))
        raise NotImplementedError

    @staticmethod
    def mul(s, t, prec, prec_mode):
        if isi(t, BinaryComplex):
            a = s.real; b = s.imag; c = t.real; d = t.imag
            return BinaryComplex((a*c-b*d, a*d+b*c))
        if isreal(t):
            return BinaryComplex((s.real*t, s.imag*t))
        raise NotImplementedError

    @staticmethod
    def div(s, t, prec, prec_mode):
        if isi(t, BinaryComplex):
            a = s.real; b = s.imag; c = t.real; d = t.imag
            mag = c*c + d*d
            return BinaryComplex(((a*c+b*d)/mag, (b*c-a*d)/mag))
        if isreal(t):
            return BinaryComplex((s.real/t, s.imag/t))
        raise NotImplementedError

    @staticmethod
    def lshift(s, n, prec, prec_mode):
        return BinaryComplex((s.real<<n, s.imag<<n))

    @staticmethod
    def rshift(s, n, prec, prec_mode):
        return BinaryComplex((s.real>>n, s.imag>>n))




#------------------------------------------------------------------
#
# Transcendental functions
#

def pif():
    context.enter(prec=context.prec+10, prec_mode=FIXED)
    def acot(x):
        s = w = BinaryReal(1)/x
        x = x**2
        n = 3
        while 1:
            w /= x
            term = w / n
            if not term:
                break
            if n & 2: s -= term
            else:     s += term
            n += 2
        return s
    p = 4*(4*acot(5)-acot(239))
    context.revert()
    return +p

def expf(z):
    r = 32
    context.enter(prec=context.prec+r+10, prec_mode=FIXED)
    z = BinaryComplex(z) >> r
    a = w = BinaryComplex(1)
    k = 1
    while a:
        a = a * z / k
        w += a
        k += 1
    for i in range(r):
        w = w*w
    if z.imag == 0:
        w = w.real
    context.revert()
    return +w

def cosf(z):
    z = BinaryComplex(z)
    eiz = expf(z * BinaryComplex((0, 1)))
    w = (eiz + BinaryComplex(1)/eiz)/2
    if z.imag == 0:
        w = w.real
    return +w

def sinf(z):
    z = BinaryComplex(z)
    eiz = expf(z * BinaryComplex((0, 1)))
    # XXX: this sum formula is bad for z ~= 0
    w = (eiz - BinaryComplex(1)/eiz) / BinaryComplex((0, 2))
    if z.imag == 0:
        w = w.real
    return +w

def newton_polish(f, r0, prec, start_prec):
    def quadratic_steps(start, target):
        L = [target + 2]
        while L[-1] > start*2:
            L = L + [L[-1]//2 + 1]
        return L[::-1]
    r = r0
    for p in quadratic_steps(start_prec, prec):
        context.enter(prec=p)
        r = f(r, p)
        context.revert()
    return r

def logf(z):
    import cmath
    prec = context.prec
    z = BinaryComplex(z, prec=prec+5)
    r = BinaryComplex(cmath.log(complex(z)), prec=53)
    start_prec = 50
    def f(r, prc):
        return r + z/expf(r) - 1
    y = newton_polish(f, r, prec, start_prec)
    y = +y
    if (isi(z, BinaryReal) and z > 0) or (z.imag == 0 and z.real > 0):
        return y.real
    else:
        return y

def powerf(x, y):
    x = BinaryComplex(x)
    y = BinaryComplex(y)
    context.enter(prec=context.prec+15)
    z = expf(logf(x) * y)
    context.revert()
    return +z
