"""
This module is experimental. There are likely bugs,
and the interface is likely to change.
------------------------------------------------------------------

This module implements arbitrary-precision binary floating-point
and fixed-point arithmetic. Due to performing all internal
arithmetic on long Python integers instead of lists of
single-digit integers, the arithmetic implemented here is
typically 10-100 times faster than decimal.py.

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


######################################################################
#
#                     Low-level arithmetic types
#
######################################################################


# Utilities

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


#----------------------------------------------------------------------
# BinaryReal
#

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
                raise NotImplementedError
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

    def __int__(s):
        try:
            return int(float(s))
        except OverflowError:
            return s.man << s.exp

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
                # TODO: do this faster
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


#----------------------------------------------------------------------
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
            # TODO: faster addition would make this redundant
            if b == d == 0:
                return BinaryComplex((a*c, 0))
            else:
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



######################################################################
#
#                     Transcendental functions
#
######################################################################


# Utilities

def make_binary(x):
    if isreal(x) or (isinstance(x, str) and 'j' not in x):
        return BinaryReal(x)
    else:
        return BinaryComplex(x)

def make_arg_binary(f):
    def g(x):
        return f(make_binary(x))
    return g

def real_to_real(f):
    def g(x):
        y = f(x)
        if isinstance(y, BinaryComplex) and isinstance(x, BinaryReal):
            return y.real
        return y
    return g

def constmemo(f):
    f.memo_prec = -1
    f.memo_val = None
    def calc():
        if context.prec <= f.memo_prec:
            return +f.memo_val
        f.memo_val = f()
        f.memo_prec = context.prec
        return +f.memo_val
    return calc

def extraprec(n):
    def decorator(f):
        def g(*args, **kwargs):
            context.enter(prec=context.prec+n)
            y = f(*args, **kwargs)
            context.revert()
            return +y
        return g
    return decorator

_i = BinaryComplex((0,1))

def _divmod(x, c):
    if isi(x, BinaryReal):
        n = int(x / c)
        t = x - n*c
        return n, t
    if isi(x, BinaryComplex):
        n = int(x.real / c)
        t = x - n*c
        return n, t

#----------------------------------------------------------------------
# Mathematical constants
#

def machin(coefs, hyperbolic=False):
    def acot(x):
        s = w = BinaryReal(1)/x
        x = x**2
        n = 3
        while 1:
            w /= x
            term = w / n
            if not term:
                break
            if hyperbolic or n & 2 == 0: s += term
            else: s -= term
            n += 2
        return s
    context.enter(prec=context.prec+15, prec_mode=FIXED)
    s = 0
    for a, b in coefs:
        s += a*acot(b)
    context.revert()
    return +s

@constmemo
def pif():
    return machin([(16, 5), (-4, 239)])

@constmemo
def log2f():
    return machin([(18, 26), (-2, 4801), (8, 8749)], True)

@constmemo
def log10f():
    return machin([(46, 31), (34, 49), (20, 161)], True)

@constmemo
def sqrt2f():
    context.enter(prec=context.prec+10, prec_mode=FIXED)
    # Newton's method
    x = half = BinaryReal(1)>>1
    d = 1
    while d:
        d = x * (half - x*x)
        x += d
    x <<= 1
    context.revert()
    return +x

@constmemo
def eulergammaf():
    """
    Compute a numerical approximation of Euler's constant ~= 0.577216

    We use the Brent-McMillan formula, g ~= A(n)/B(n) - log(n), where
        A(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
        B(n) = sum_{k=0,1,2,...} (n**k / k!)**2
        H(k) = 1 + 1/2 + 1/3 + ... + 1/k

    The error is bounded by O(exp(-4n)). Choosing n to be a power
    of two, 2**p, the logarithm becomes particularly easy to calculate.

    Reference:
    Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
    http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
    """
    # TODO: may need even more extra precision
    context.enter(prec=context.prec+30, prec_mode=FIXED)
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((context.prec/4) * math.log(2), 2)) + 1
    n = 1 << p
    one = BinaryReal(1)
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    r = one
    while r:
        A += r * H
        B += r
        r = r * (n*n) / (k*k)
        H += one / k
        k += 1
    S = (A / B) - p*log2f()
    context.revert()
    return +S


#----------------------------------------------------------------------
# Exponential function
#

@make_arg_binary
@extraprec(3)
def expf(x):
    """
    Calculate exp of a BinaryReal or BinaryComplex.

    If x is a BinaryReal, we first rewrite x as t + n*log(2) and then
    calculate exp(x) as exp(t)*(2**n). With |t| <= log(2) ~= 0.7,
    exp(t) can be computed very quickly from the Maclaurin series
    of exp, and multiplying a BinaryReal by 2**n costs nothing.

    If x is a BinaryComplex, we use Euler's formula. It would be
    possible to use Maclaurin series directly for complex numbers as
    well, but that results in loss of precision near a root for the
    real or imaginary part. In the future, such a method could be
    used for prec_mode=FIXED.
    """
    if isi(x, BinaryReal):
        n, t = _divmod(x, log2f())
        return exp_near_0(t) << n
    if isi(x, BinaryComplex):
        mag = expf(x.real)
        re = mag * cosf(x.imag)
        im = mag * sinf(x.imag)
        return BinaryComplex((re, im))

# Helper for exp

def exp_near_0(x, r=None):
    """
    Calculate exp(x) for a BinaryReal x that should be close to 0.
    The method will work for large x and BinaryComplex's as well,
    but less efficiently and with loss of relative precision
    around a root of the real or imaginary part.

    The basic algorithm is to sum the Maclaurin series

        exp(x) = 1 + x + x**2/2 + x**3/6 + ...

    using fixed-point arithmetic.

    To improve the rate of convergence, we choose an integer r
    and instead calculate exp(x/2**r)**(2**r) = exp(x). The optimal
    value for r depends on the Python platform, the magnitude
    of x, and the target precision, and has to be estimated
    from experimental timings. One test with x ~= 0.3 showed that
    r = 2.2*prec**0.42 gave a good fit to the optimal values for r
    for prec between 1 and 10000 bits, on one particular machine.
    This is used as a default value for r (it could be tweaked
    in the future.)

    This optimization makes the summation about twice as fast at
    low precision levels and much faster at high precision
    (roughly five times faster at 1000 decimal digits).
    """
    if r == None:
        r = int(2.2*context.prec**0.42)
    newprec = context.prec + r+8
    context.enter(prec=newprec, prec_mode=FIXED)
    a = s = BinaryReal(1)
    # TODO: should be able to do x >> r . fix __rshift__
    x = (+x) >> r
    k = 1
    while a:
        a = a*x / k
        s += a
        k += 1
    for j in range(r):
        s = s*s
    context.revert()
    return s

#----------------------------------------------------------------------
# Sine
#

@real_to_real
@make_arg_binary
def sinf(x):
    """
    To evaluate sin(x) for a real x or a complex x that lies close to
    the real axis, we first use the elementary translation and
    reflection identities for trigonometric functions to obtain
    an argument [with real part] between 0 and pi/4. Then we sum
    the Maclaurin series for either cos or sin.

    If x has a large imaginary part, we use the formula

      sin(x+i*y) = sin(x)*cosh(y) + i*cos(x)*sinh(y)

    with cosh and sinh computed from the real exp.

    The trickery is needed to preserve relative precision near a
    root. For prec_mode=FIXED, it would be possible to compute
    cos and sin more directly from the complex exponential (instead of
    vice versa).

    The other trigonometric functions are computed from sin (and exp).
    """
    if isi(x, BinaryComplex):
        re, im = x.real, x.imag
    else:
        re, im = x, 0
    if re < 0:
        return -sinf(-x)
    if abs(im) < 1:
        w = sin_near_real(x)
    else:
        context.enter(prec=context.prec+5)
        a = sin_near_real(re)
        b = sin_near_real(pif()/2 - re)
        ey = expf(im)
        eyinv = BinaryReal(1) / ey
        cosh = (ey + eyinv)/2
        sinh = (ey - eyinv)/2
        w = BinaryComplex((a*cosh, b*sinh))
        context.revert()
    return +w

# Sine helpers

@extraprec(3)
def sin_near_real(x):
    """Compute an approximation of sin(x) where x is a BinaryReal
    or a BinaryComplex. The real part of x can be arbitrarily
    large, but the imaginary part should be close to 0, or this
    function will be slow and cause loss of precision.
    """
    # Reduce argument mod pi/4 to obtain a base case
    pi4 = pif() >> 2
    n, t = _divmod(x, pi4)
    if n%2 == 1: t = pi4 - t
    if n%8 in (0, 3, 4, 7): r = sin_near_0(t)
    else: r = cos_near_0(t)
    if n%8 > 3: r = -r
    return r

def cos_near_0(x):
    """Maclaurin series approximation for cos."""
    context.enter(prec=context.prec+10, prec_mode=FIXED)
    a = s = BinaryReal(1)
    x2 = x*x
    k = 2
    while a:
        a = a*x2 / ((k-1) * k)
        if (k>>1)&1: s -= a
        else:        s += a
        k += 2
    context.revert()
    return s

def sin_near_0(x):
    """Maclaurin series approximation for sin. Unlike cos_near_0,
    this function increases the precision for x close to 0 to
    maintain a high relative precision."""
    if isi(x, BinaryComplex):
        extraprec = max(0, bitcount(x.real.man)-x.real.exp)
    else:
        extraprec = max(0, bitcount(x.man)-x.exp)
    context.enter(prec=context.prec+10+extraprec, prec_mode=FIXED)
    a = s = +x # x?
    x2 = x*x
    k = 3
    while a:
        a = a * x2 / ((k-1) * k)
        if (k>>1)&1: s -= a
        else:        s += a
        k += 2
    context.revert()
    return s


#----------------------------------------------------------------------
# Additional trigonometric and hyperbolic functions
#

@make_arg_binary
@extraprec(3)
def cosf(x):
    """Calculate cos of a BinaryReal or BinaryComplex"""
    return sinf(x + pif()/2)

@make_arg_binary
@extraprec(3)
def tanf(x):
    """Calculate tan of a BinaryReal or BinaryComplex"""
    return sinf(x)/cosf(x)

@make_arg_binary
def sinhf(x):
    """Calculate sinh of a BinaryReal or BinaryComplex"""
    return -_i*sinf(_i*x)

@make_arg_binary
def coshf(x):
    """Calculate cosh of a BinaryReal or BinaryComplex"""
    return cosf(_i*x)

@make_arg_binary
def tanhf(x):
    """Calculate tanh of a BinaryReal or BinaryComplex"""
    return -_i*tanf(_i*x)


#----------------------------------------------------------------------
# Inverse functions etc (much more work needed here)
#

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
    # TODO: argument reduction
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

@extraprec(15)
def powerf(x, y):
    return expf(logf(x) * y)
