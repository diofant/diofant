"""
Adaptive numerical evaluation of Diofant expressions, using mpmath
for mathematical functions.
"""

import math
import numbers

import mpmath.libmp as libmp
from mpmath import inf as mpmath_inf
from mpmath import (make_mpc, make_mpf, mp, mpc, mpf, nsum, quadosc, quadts,
                    workprec)
from mpmath.libmp import bitcount as mpmath_bitcount
from mpmath.libmp import (fhalf, fnan, fnone, fone, from_int, from_man_exp,
                          from_rational, fzero, mpf_abs, mpf_add, mpf_atan,
                          mpf_atan2, mpf_cmp, mpf_cos, mpf_e, mpf_exp, mpf_log,
                          mpf_lt, mpf_mul, mpf_neg, mpf_pi, mpf_pow,
                          mpf_pow_int, mpf_shift, mpf_sin, mpf_sqrt, normalize,
                          round_nearest)
from mpmath.libmp.backend import MPZ
from mpmath.libmp.gammazeta import mpf_bernoulli
from mpmath.libmp.libmpc import _infs_nan
from mpmath.libmp.libmpf import dps_to_prec, prec_to_dps

from .compatibility import is_sequence
from .singleton import S
from .sympify import sympify


LG10 = math.log(10, 2)
rnd = round_nearest


def bitcount(n):
    return mpmath_bitcount(int(n))


# Used in a few places as placeholder values to denote exponents and
# precision levels, e.g. of exact numbers. Must be careful to avoid
# passing these to mpmath functions or returning them in final results.
INF = float(mpmath_inf)
MINUS_INF = float(-mpmath_inf)

# ~= 100 digits. Real men set this to INF.
DEFAULT_MAXPREC = int(110*LG10)  # keep in sync with maxn kwarg of evalf


class PrecisionExhausted(ArithmeticError):
    """Raised when precision is exhausted."""


############################################################################
#                                                                          #
#              Helper functions for arithmetic and complex parts           #
#                                                                          #
############################################################################


"""
An mpf value tuple is a tuple of integers (sign, man, exp, bc)
representing a floating-point number: [1, -1][sign]*man*2**exp where
sign is 0 or 1 and bc should correspond to the number of bits used to
represent the mantissa (man) in binary notation, e.g.

>>> sign, man, exp, bc = 0, 5, 1, 3
>>> n = [1, -1][sign]*man*2**exp
>>> n, bitcount(man)
(10, 3)

A temporary result is a tuple (re, im, re_acc, im_acc) where
re and im are nonzero mpf value tuples representing approximate
numbers, or None to denote exact zeros.

re_acc, im_acc are integers denoting log2(e) where e is the estimated
relative accuracy of the respective complex part, but may be anything
if the corresponding complex part is None.

"""


def fastlog(x):
    """Fast approximation of log2(x) for an mpf value tuple x.

    Notes: Calculated as exponent + width of mantissa. This is an
    approximation for two reasons: 1) it gives the ceil(log2(abs(x)))
    value and 2) it is too high by 1 in the case that x is an exact
    power of 2. Although this is easy to remedy by testing to see if
    the odd mpf mantissa is 1 (indicating that one was dealing with
    an exact power of 2) that would decrease the speed and is not
    necessary as this is only being used as an approximation for the
    number of bits in x. The correct return value could be written as
    "x[2] + (x[3] if x[1] != 1 else 0)".
        Since mpf tuples always have an odd mantissa, no check is done
    to see if the mantissa is a multiple of 2 (in which case the
    result would be too large by 1).

    Examples
    ========

    >>> s, m, e = 0, 5, 1
    >>> bc = bitcount(m)
    >>> n = [1, -1][s]*m*2**e
    >>> n, (log(n)/log(2)).evalf(2), fastlog((s, m, e, bc))
    (10, 3.3, 4)

    """
    if not x or x == fzero:
        return MINUS_INF
    return x[2] + x[3]


def pure_complex(v):
    """Return a and b if v matches a + I*b where b is not zero and
    a and b are Numbers, else None.

    >>> a, b = Tuple(2, 3)
    >>> pure_complex(a)
    >>> pure_complex(a + b*I)
    (2, 3)
    >>> pure_complex(I)
    (0, 1)

    """
    from .numbers import I
    h, t = v.as_coeff_Add()
    c, i = t.as_coeff_Mul()
    if i is I:
        return h, c


def scaled_zero(mag, sign=1):
    """Return an mpf representing a power of two with magnitude ``mag``
    and -1 for precision. Or, if ``mag`` is a scaled_zero tuple, then just
    remove the sign from within the list that it was initially wrapped
    in.

    Examples
    ========

    >>> z, p = scaled_zero(100)
    >>> z, p
    (([0], 1, 100, 1), -1)
    >>> ok = scaled_zero(z)
    >>> ok
    (0, 1, 100, 1)
    >>> Float(ok)
    1.26765060022823e+30
    >>> Float(ok, p)
    0.e+30
    >>> ok, p = scaled_zero(100, -1)
    >>> Float(scaled_zero(ok), p)
    -0.e+30

    """
    if type(mag) is tuple and len(mag) == 4 and iszero(mag, scaled=True):
        return (mag[0][0],) + mag[1:]
    elif isinstance(mag, numbers.Integral):
        if sign not in [-1, 1]:
            raise ValueError('sign must be +/-1')
        rv, p = mpf_shift(fone, mag), -1
        s = 0 if sign == 1 else 1
        rv = ([s],) + rv[1:]
        return rv, p
    else:
        raise ValueError('scaled zero expects int or scaled_zero tuple.')


def iszero(mpf, scaled=False):
    if not scaled:
        return not mpf or not mpf[1] and not mpf[-1]
    return mpf and type(mpf[0]) is list and mpf[1] == mpf[-1] == 1


def complex_accuracy(result):
    """
    Returns relative accuracy of a complex number with given accuracies
    for the real and imaginary parts. The relative accuracy is defined
    in the complex norm sense as ||z|+|error|| / |z| where error
    is equal to (real absolute error) + (imag absolute error)*i.

    The full expression for the (logarithmic) error can be approximated
    easily by using the max norm to approximate the complex norm.

    In the worst case (re and im equal), this is wrong by a factor
    sqrt(2), or by log2(sqrt(2)) = 0.5 bit.

    """
    re, im, re_acc, im_acc = result
    if not im:
        if not re:
            return INF
        return re_acc
    if not re:
        return im_acc
    re_size = fastlog(re)
    im_size = fastlog(im)
    absolute_error = max(re_size - re_acc, im_size - im_acc)
    relative_error = absolute_error - max(re_size, im_size)
    return -relative_error


def get_abs(expr, prec, options):
    re, im, re_acc, im_acc = evalf(expr, prec + 2, options)
    if not re:
        re, re_acc, im, im_acc = im, im_acc, re, re_acc
    if im:
        return libmp.mpc_abs((re, im), prec), None, re_acc, None
    elif re:
        return mpf_abs(re), None, re_acc, None
    else:
        return None, None, None, None


def get_complex_part(expr, no, prec, options):
    """Selector no = 0 for real part, no = 1 for imaginary part."""
    workprec = prec
    i = 0
    while 1:
        res = evalf(expr, workprec, options)
        value, accuracy = res[no::2]
        if (not value) or accuracy >= prec or expr.is_Float:
            return value, None, accuracy, None
        workprec += max(30, 2**i)
        i += 1


def evalf_abs(expr, prec, options):
    return get_abs(expr.args[0], prec, options)


def evalf_re(expr, prec, options):
    return get_complex_part(expr.args[0], 0, prec, options)


def evalf_im(expr, prec, options):
    return get_complex_part(expr.args[0], 1, prec, options)


def finalize_complex(re, im, prec):
    assert re != fzero or im != fzero

    if re == fzero:
        return None, im, None, prec
    elif im == fzero:
        return re, None, prec, None

    size_re = fastlog(re)
    size_im = fastlog(im)
    if size_re > size_im:
        re_acc = prec
        im_acc = prec + min(-(size_re - size_im), 0)
    else:
        im_acc = prec
        re_acc = prec + min(-(size_im - size_re), 0)
    return re, im, re_acc, im_acc


def chop_parts(value, prec):
    """Chop off tiny real or complex parts."""
    re, im, re_acc, im_acc = value
    # chop based on absolute value
    if re and re not in _infs_nan and (fastlog(re) < -prec + 4):
        re, re_acc = None, None
    if im and im not in _infs_nan and (fastlog(im) < -prec + 4):
        im, im_acc = None, None
    return re, im, re_acc, im_acc


def check_target(expr, result, prec):
    a = complex_accuracy(result)
    if a < prec:
        raise PrecisionExhausted('Failed to distinguish the expression: \n\n%s\n\n'
                                 'from zero. Try simplifying the input, using chop=True, or providing '
                                 'a higher maxn for evalf' % expr)


############################################################################
#                                                                          #
#                            Arithmetic operations                         #
#                                                                          #
############################################################################


def add_terms(terms, prec, target_prec):
    """
    Helper for evalf_add. Adds a list of (mpfval, accuracy) terms.

    Returns
    =======

    - None, None if there are no non-zero terms;
    - terms[0] if there is only 1 term;
    - scaled_zero if the sum of the terms produces a zero by cancellation
      e.g. mpfs representing 1 and -1 would produce a scaled zero which need
      special handling since they are not actually zero and they are purposely
      malformed to ensure that they can't be used in anything but accuracy
      calculations;
    - a tuple that is scaled to target_prec that corresponds to the
      sum of the terms.

    The returned mpf tuple will be normalized to target_prec; the input
    prec is used to define the working precision.

    XXX explain why this is needed and why one can't just loop using mpf_add

    """
    terms = [t for t in terms if not iszero(t)]
    if not terms:
        return None, None
    elif len(terms) == 1:
        return terms[0]

    # see if any argument is NaN or oo and thus warrants a special return
    special = []
    from .numbers import Float, nan
    for t in terms:
        arg = Float._new(t[0], 1)
        if arg is nan or arg.is_infinite:
            special.append(arg)
    if special:
        from .add import Add
        rv = evalf(Add(*special), prec + 4, {})
        return rv[0], rv[2]

    working_prec = 2*prec
    sum_man, sum_exp, absolute_error = 0, 0, MINUS_INF

    for x, accuracy in terms:
        sign, man, exp, bc = x
        if sign:
            man = -man
        absolute_error = max(absolute_error, bc + exp - accuracy)
        delta = exp - sum_exp
        if exp >= sum_exp:
            # x much larger than existing sum?
            # first: quick test
            if ((delta > working_prec) and
                ((not sum_man) or
                 delta - bitcount(abs(sum_man)) > working_prec)):
                sum_man = man
                sum_exp = exp
            else:
                sum_man += (man << delta)
        else:
            delta = -delta
            # x much smaller than existing sum?
            if delta - bc > working_prec:
                if not sum_man:
                    sum_man, sum_exp = man, exp
            else:
                sum_man = (sum_man << delta) + man
                sum_exp = exp
    if not sum_man:
        return scaled_zero(absolute_error)
    if sum_man < 0:
        sum_sign = 1
        sum_man = -sum_man
    else:
        sum_sign = 0
    sum_bc = bitcount(sum_man)
    sum_accuracy = sum_exp + sum_bc - absolute_error
    r = normalize(sum_sign, sum_man, sum_exp, sum_bc, target_prec,
                  rnd), sum_accuracy
    return r


def evalf_add(v, prec, options):
    res = pure_complex(v)
    if res:
        h, c = res
        re, _, re_acc, _ = evalf(h, prec, options)
        im, _, im_acc, _ = evalf(c, prec, options)
        return re, im, re_acc, im_acc

    oldmaxprec = options['maxprec']

    i = 0
    target_prec = prec
    while 1:
        options['maxprec'] = min(oldmaxprec, 2*prec)

        terms = [evalf(arg, prec + 10, options) for arg in v.args]
        re, re_acc = add_terms(
            [a[0::2] for a in terms if a[0]], prec, target_prec)
        im, im_acc = add_terms(
            [a[1::2] for a in terms if a[1]], prec, target_prec)
        acc = complex_accuracy((re, im, re_acc, im_acc))
        if acc >= target_prec:
            break
        else:
            if (prec - target_prec) > options['maxprec']:
                break

            prec = prec + max(10 + 2**i, target_prec - acc)
            i += 1

    options['maxprec'] = oldmaxprec
    if iszero(re, scaled=True):
        re = scaled_zero(re)
    if iszero(im, scaled=True):
        im = scaled_zero(im)
    return re, im, re_acc, im_acc


def evalf_mul(v, prec, options):
    res = pure_complex(v)
    if res:
        # the only pure complex that is a mul is h*I
        _, h = res
        im, _, im_acc, _ = evalf(h, prec, options)
        return None, im, None, im_acc
    args = list(v.args)

    # see if any argument is NaN or oo and thus warrants a special return
    special, other = [], []
    from .numbers import Float, nan
    for arg in args:
        arg = evalf(arg, prec, options)
        if arg[0] is None:
            continue
        arg = Float._new(arg[0], 1)
        if arg is nan or arg.is_infinite:
            special.append(arg)
        else:
            other.append(arg)
    if special:
        from .mul import Mul
        other = Mul(*other)
        special = Mul(*special)
        return evalf(special*other, prec + 4, {})

    # With guard digits, multiplication in the real case does not destroy
    # accuracy. This is also true in the complex case when considering the
    # total accuracy; however accuracy for the real or imaginary parts
    # separately may be lower.
    acc = prec

    # XXX: big overestimate
    working_prec = prec + len(args) + 5

    # Empty product is 1
    start = man, exp, bc = MPZ(1), 0, 1

    # First, we multiply all pure real or pure imaginary numbers.
    # direction tells us that the result should be multiplied by
    # I**direction; all other numbers get put into complex_factors
    # to be multiplied out after the first phase.
    last = len(args)
    direction = 0
    args.append(S.One)
    complex_factors = []

    for i, arg in enumerate(args):
        if i != last and pure_complex(arg):
            args[-1] = (args[-1]*arg).expand()
            continue
        elif i == last and arg is S.One:
            continue
        re, im, re_acc, im_acc = evalf(arg, working_prec, options)
        if re and im:
            complex_factors.append((re, im, re_acc, im_acc))
            continue
        elif re:
            (s, m, e, b), w_acc = re, re_acc
        elif im:
            (s, m, e, b), w_acc = im, im_acc
            direction += 1
        else:
            return None, None, None, None
        direction += 2*s
        man *= m
        exp += e
        bc += b
        if bc > 3*working_prec:
            man >>= working_prec
            exp += working_prec
        acc = min(acc, w_acc)
    sign = (direction & 2) >> 1
    if not complex_factors:
        v = normalize(sign, man, exp, bitcount(man), prec, rnd)
        # multiply by i
        if direction & 1:
            return None, v, None, acc
        else:
            return v, None, acc, None
    else:
        # initialize with the first term
        if (man, exp, bc) != start:
            # there was a real part; give it an imaginary part
            re, im = (sign, man, exp, bitcount(man)), (0, MPZ(0), 0, 0)
            i0 = 0
        else:
            # there is no real part to start (other than the starting 1)
            wre, wim, wre_acc, wim_acc = complex_factors[0]
            acc = min(acc,
                      complex_accuracy((wre, wim, wre_acc, wim_acc)))
            re = wre
            im = wim
            i0 = 1

        for wre, wim, wre_acc, wim_acc in complex_factors[i0:]:
            # acc is the overall accuracy of the product; we aren't
            # computing exact accuracies of the product.
            acc = min(acc,
                      complex_accuracy((wre, wim, wre_acc, wim_acc)))

            use_prec = working_prec
            A = mpf_mul(re, wre, use_prec)
            B = mpf_mul(mpf_neg(im), wim, use_prec)
            C = mpf_mul(re, wim, use_prec)
            D = mpf_mul(im, wre, use_prec)
            re = mpf_add(A, B, use_prec)
            im = mpf_add(C, D, use_prec)
        # multiply by I
        if direction & 1:
            re, im = mpf_neg(im), re
        return re, im, acc, acc


def evalf_pow(v, prec, options):
    from .numbers import E

    target_prec = prec
    base, exp = v.args

    # We handle x**n separately. This has two purposes: 1) it is much
    # faster, because we avoid calling evalf on the exponent, and 2) it
    # allows better handling of real/imaginary parts that are exactly zero
    if exp.is_Integer:
        p = exp.numerator
        # Exact
        if not p:
            return fone, None, prec, None
        # Exponentiation by p magnifies relative error by |p|, so the
        # base must be evaluated with increased precision if p is large
        prec += int(math.log(abs(p), 2))
        re, im, re_acc, im_acc = evalf(base, prec + 5, options)
        # Real to integer power
        if re and not im:
            return mpf_pow_int(re, p, target_prec), None, target_prec, None
        # (x*I)**n = I**n * x**n
        if im and not re:
            z = mpf_pow_int(im, p, target_prec)
            case = p % 4
            if case == 0:
                return z, None, target_prec, None
            elif case == 1:
                return None, z, None, target_prec
            elif case == 2:
                return mpf_neg(z), None, target_prec, None
            else:
                return None, mpf_neg(z), None, target_prec
        # Zero raised to an integer power
        if not re:
            return None, None, None, None
        # General complex number to arbitrary integer power
        re, im = libmp.mpc_pow_int((re, im), p, prec)
        # Assumes full accuracy in input
        return finalize_complex(re, im, target_prec)

    # Pure square root
    if exp is S.Half:
        xre, xim, _, _ = evalf(base, prec + 5, options)
        # General complex square root
        if xim:
            re, im = libmp.mpc_sqrt((xre or fzero, xim), prec)
            return finalize_complex(re, im, prec)
        if not xre:
            return None, None, None, None
        # Square root of a negative real number
        if mpf_lt(xre, fzero):
            return None, mpf_sqrt(mpf_neg(xre), prec), None, prec
        # Positive square root
        return mpf_sqrt(xre, prec), None, prec, None

    # We first evaluate the exponent to find its magnitude
    # This determines the working precision that must be used
    prec += 10
    yre, yim, _, _ = evalf(exp, prec, options)
    # Special cases: x**0
    if not (yre or yim):
        return fone, None, prec, None

    ysize = fastlog(yre)
    # Restart if too big
    # XXX: prec + ysize might exceed maxprec
    if ysize > 5:
        prec += ysize
        yre, yim, _, _ = evalf(exp, prec, options)

    # Pure exponential function; no need to evalf the base
    if base is E:
        if yim:
            re, im = libmp.mpc_exp((yre or fzero, yim), prec)
            return finalize_complex(re, im, target_prec)
        return mpf_exp(yre, target_prec), None, target_prec, None

    xre, xim, _, _ = evalf(base, prec + 5, options)
    # 0**y
    if not (xre or xim):
        return None, None, None, None

    # (real ** complex) or (complex ** complex)
    if yim:
        re, im = libmp.mpc_pow(
            (xre or fzero, xim or fzero), (yre or fzero, yim),
            target_prec)
        return finalize_complex(re, im, target_prec)
    # complex ** real
    if xim:
        re, im = libmp.mpc_pow_mpf((xre or fzero, xim), yre, target_prec)
        return finalize_complex(re, im, target_prec)
    # negative ** real
    elif mpf_lt(xre, fzero):
        re, im = libmp.mpc_pow_mpf((xre, fzero), yre, target_prec)
        return finalize_complex(re, im, target_prec)
    # positive ** real
    else:
        return mpf_pow(xre, yre, target_prec), None, target_prec, None


############################################################################
#                                                                          #
#                            Special functions                             #
#                                                                          #
############################################################################
def evalf_trig(v, prec, options):
    """
    This function handles sin and cos of complex arguments.

    TODO: should also handle tan of complex arguments.

    """
    from ..functions import cos, sin
    if isinstance(v, cos):
        func = mpf_cos
    elif isinstance(v, sin):
        func = mpf_sin
    else:
        raise NotImplementedError
    arg = v.args[0]
    # 20 extra bits is possibly overkill. It does make the need
    # to restart very unlikely
    xprec = prec + 20
    re, im, re_acc, im_acc = evalf(arg, xprec, options)
    if im:
        if 'subs' in options:
            v = v.subs(options['subs'])
        return evalf(v._eval_evalf(prec), prec, options)
    if not re:
        if isinstance(v, cos):
            return fone, None, prec, None
        elif isinstance(v, sin):
            return None, None, None, None
        else:
            raise NotImplementedError
    # For trigonometric functions, we are interested in the
    # fixed-point (absolute) accuracy of the argument.
    xsize = fastlog(re)
    # Magnitude <= 1.0. OK to compute directly, because there is no
    # danger of hitting the first root of cos (with sin, magnitude
    # <= 2.0 would actually be ok)
    if xsize < 1:
        return func(re, prec, rnd), None, prec, None
    # Very large
    if xsize >= 10:
        xprec = prec + xsize
        re, im, re_acc, im_acc = evalf(arg, xprec, options)
    # Need to repeat in case the argument is very close to a
    # multiple of pi (or pi/2), hitting close to a root
    while 1:
        y = func(re, prec, rnd)
        ysize = fastlog(y)
        gap = -ysize
        accuracy = (xprec - xsize) - gap
        if accuracy < prec:
            if xprec > options['maxprec']:
                return y, None, accuracy, None
            xprec += gap
            re, im, re_acc, im_acc = evalf(arg, xprec, options)
            continue
        else:
            return y, None, prec, None


def evalf_log(expr, prec, options):
    from ..functions import Abs, log
    from .add import Add

    if len(expr.args) > 1:
        expr = expr.doit()
        return evalf(expr, prec, options)

    arg = expr.args[0]
    workprec = prec + 10
    xre, xim, xacc, _ = evalf(arg, workprec, options)

    if xim:
        # XXX: use get_abs etc instead
        re = evalf_log(
            log(Abs(arg, evaluate=False), evaluate=False), prec, options)
        im = mpf_atan2(xim, xre or fzero, prec)
        return re[0], im, re[2], prec

    imaginary_term = (mpf_cmp(xre, fzero) < 0)

    re = mpf_log(mpf_abs(xre), prec, rnd)
    size = fastlog(re)
    if prec - size > workprec:
        # We actually need to compute 1+x accurately, not x
        arg = Add(S.NegativeOne, arg, evaluate=False)
        xre, xim, _, _ = evalf_add(arg, prec, options)
        prec2 = workprec - fastlog(xre)
        # xre is now x - 1 so we add 1 back here to calculate x
        re = mpf_log(mpf_abs(mpf_add(xre, fone, prec2)), prec, rnd)

    re_acc = prec

    if imaginary_term:
        return re, mpf_pi(prec), re_acc, prec
    else:
        return re, None, re_acc, None


def evalf_atan(v, prec, options):
    arg = v.args[0]
    xre, xim, reacc, imacc = evalf(arg, prec + 5, options)
    if xre is xim is None:
        return (None,)*4
    if xim:
        raise NotImplementedError
    return mpf_atan(xre, prec, rnd), None, prec, None


def evalf_subs(prec, subs):
    """Change all Float entries in `subs` to have precision prec."""
    newsubs = {}
    for a, b in subs.items():
        b = sympify(b)
        if b.is_Float:
            b = b._eval_evalf(prec)
        newsubs[a] = b
    return newsubs


def evalf_piecewise(expr, prec, options):
    if 'subs' in options:
        expr = expr.subs(evalf_subs(prec, options['subs']))
        newopts = options.copy()
        del newopts['subs']
        return evalf(expr, prec, newopts)

    # We still have undefined symbols
    raise NotImplementedError


def evalf_bernoulli(expr, prec, options):
    arg = expr.args[0]
    if not arg.is_Integer:
        raise ValueError('Bernoulli number index must be an integer')
    n = int(arg)
    b = mpf_bernoulli(n, prec, rnd)
    if b == fzero:
        return None, None, None, None
    return b, None, prec, None

############################################################################
#                                                                          #
#                            High-level operations                         #
#                                                                          #
############################################################################


def as_mpmath(x, prec, options):
    from .numbers import oo
    x = sympify(x)
    if x == 0:
        return mpf(0)
    if x == oo:
        return mpf('inf')
    if x == -oo:
        return mpf('-inf')
    # XXX
    re, im, _, _ = evalf(x, prec, options)
    if im:
        return mpc(re or fzero, im)
    return mpf(re)


def do_integral(expr, prec, options):
    func = expr.args[0]
    x, xlow, xhigh = expr.args[1]
    if xlow == xhigh:
        xlow = xhigh = 0
    elif x not in func.free_symbols:
        # only the difference in limits matters in this case
        # so if there is a symbol in common that will cancel
        # out when taking the difference, then use that
        # difference
        if xhigh.free_symbols & xlow.free_symbols:
            diff = xhigh - xlow
            if not diff.free_symbols:
                xlow, xhigh = 0, diff

    oldmaxprec = options['maxprec']
    options['maxprec'] = min(oldmaxprec, 2*prec)

    with workprec(prec + 5):
        xlow = as_mpmath(xlow, prec + 15, options)
        xhigh = as_mpmath(xhigh, prec + 15, options)

        # Integration is like summation, and we can phone home from
        # the integrand function to update accuracy summation style
        # Note that this accuracy is inaccurate, since it fails
        # to account for the variable quadrature weights,
        # but it is better than nothing

        from ..functions import cos, sin
        from .numbers import pi
        from .symbol import Wild

        have_part = [False, False]
        max_real_term = [MINUS_INF]
        max_imag_term = [MINUS_INF]

        def f(t):
            re, im, re_acc, im_acc = evalf(func, mp.prec, {'subs': {x: t}, 'maxprec': DEFAULT_MAXPREC})

            have_part[0] = re or have_part[0]
            have_part[1] = im or have_part[1]

            max_real_term[0] = max(max_real_term[0], fastlog(re))
            max_imag_term[0] = max(max_imag_term[0], fastlog(im))

            if im:
                return mpc(re or fzero, im)
            return mpf(re or fzero)

        if options.get('quad') == 'osc':
            A = Wild('A', exclude=[x])
            B = Wild('B', exclude=[x])
            D = Wild('D')
            m = func.match(cos(A*x + B)*D)
            if not m:
                m = func.match(sin(A*x + B)*D)
            if not m:
                raise ValueError('An integrand of the form sin(A*x+B)*f(x) '
                                 'or cos(A*x+B)*f(x) is required for oscillatory quadrature')
            period = as_mpmath(2*pi/m[A], prec + 15, options)
            result = quadosc(f, [xlow, xhigh], period=period)
            # XXX: quadosc does not do error detection yet
            quadrature_error = MINUS_INF
        else:
            result, quadrature_error = quadts(f, [xlow, xhigh], error=1)
            quadrature_error = fastlog(quadrature_error._mpf_)

    options['maxprec'] = oldmaxprec

    if have_part[0]:
        re = result.real._mpf_
        if re == fzero:
            re, re_acc = scaled_zero(
                min(-prec, -max_real_term[0], -quadrature_error))
            re = scaled_zero(re)  # handled ok in evalf_integral
        else:
            re_acc = -max(max_real_term[0] - fastlog(re) -
                          prec, quadrature_error)
    else:
        re, re_acc = None, None

    if have_part[1]:
        im = result.imag._mpf_
        if im == fzero:
            im, im_acc = scaled_zero(
                min(-prec, -max_imag_term[0], -quadrature_error))
            im = scaled_zero(im)  # handled ok in evalf_integral
        else:
            im_acc = -max(max_imag_term[0] - fastlog(im) -
                          prec, quadrature_error)
    else:
        im, im_acc = None, None

    result = re, im, re_acc, im_acc
    return result


def evalf_integral(expr, prec, options):
    limits = expr.limits
    if len(limits) != 1 or len(limits[0]) != 3:
        raise NotImplementedError
    workprec = prec
    i = 0
    maxprec = options.get('maxprec', INF)
    while 1:
        result = do_integral(expr, workprec, options)
        accuracy = complex_accuracy(result)
        if accuracy >= prec:  # achieved desired precision
            break
        if workprec >= maxprec:  # can't increase accuracy any more
            break
        if accuracy == -1:
            # maybe the answer really is zero and maybe we just haven't increased
            # the precision enough. So increase by doubling to not take too long
            # to get to maxprec.
            workprec *= 2
        else:
            workprec += max(prec, 2**i)
        workprec = min(workprec, maxprec)
        i += 1
    return result


def check_convergence(numer, denom, n):
    """
    Returns (h, g, p) where
    -- h is:
        > 0 for convergence of rate 1/factorial(n)**h
        < 0 for divergence of rate factorial(n)**(-h)
        = 0 for geometric or polynomial convergence or divergence

    -- abs(g) is:
        > 1 for geometric convergence of rate 1/h**n
        < 1 for geometric divergence of rate h**n
        = 1 for polynomial convergence or divergence

        (g < 0 indicates an alternating series)

    -- p is:
        > 1 for polynomial convergence of rate 1/n**h
        <= 1 for polynomial divergence of rate n**(-h)

    """
    npol = numer.as_poly(n)
    dpol = denom.as_poly(n)
    p = npol.degree()
    q = dpol.degree()
    rate = q - p
    if rate:
        return rate, None, None
    constant = dpol.LC() / npol.LC()
    if abs(constant) != 1:
        return rate, constant, None
    if npol.degree() == dpol.degree() == 0:
        return rate, constant, 0
    pc = npol.all_coeffs()[-2]
    qc = dpol.all_coeffs()[-2]
    return rate, constant, (qc - pc)/dpol.LC()


def hypsum(expr, n, start, prec):
    """
    Sum a rapidly convergent infinite hypergeometric series with
    given general term, e.g. e = hypsum(1/factorial(n), n). The
    quotient between successive terms must be a quotient of integer
    polynomials.

    """
    from ..simplify import hypersimp
    from ..utilities import lambdify
    from .numbers import Float

    if prec == float('inf'):
        raise NotImplementedError('does not support inf prec')

    if start:
        expr = expr.subs({n: n + start})
    hs = hypersimp(expr, n)
    if hs is None:
        raise NotImplementedError('a hypergeometric series is required')
    num, den = hs.as_numer_denom()

    func1 = lambdify(n, num)
    func2 = lambdify(n, den)

    h, g, p = check_convergence(num, den, n)

    if h < 0:
        raise ValueError('Sum diverges like (n!)^%i' % (-h))

    term = expr.subs({n: 0})
    if not term.is_Rational:
        raise NotImplementedError('Non rational term functionality is not implemented.')

    # Direct summation if geometric or faster
    if h > 0 or (h == 0 and abs(g) > 1):
        term = (MPZ(term.numerator) << prec) // term.denominator
        s = term
        k = 1
        while abs(term) > 5:
            term *= MPZ(func1(k - 1))
            term //= MPZ(func2(k - 1))
            s += term
            k += 1
        return from_man_exp(s, -prec)
    else:
        alt = g < 0
        if abs(g) < 1:
            raise ValueError('Sum diverges like (%i)^n' % abs(1/g))
        if p < 1 or (p == 1 and not alt):
            raise ValueError('Sum diverges like n^%i' % (-p))
        # We have polynomial convergence: use Richardson extrapolation
        vold = None
        ndig = prec_to_dps(prec)
        while True:
            # Need to use at least quad precision because a lot of cancellation
            # might occur in the extrapolation process; we check the answer to
            # make sure that the desired precision has been reached, too.
            prec2 = 4*prec
            term0 = (MPZ(term.numerator) << prec2) // term.denominator

            def summand(k, _term=[term0]):
                if k:
                    k = int(k)
                    _term[0] *= MPZ(func1(k - 1))
                    _term[0] //= MPZ(func2(k - 1))
                return make_mpf(from_man_exp(_term[0], -prec2))

            with workprec(prec):
                v = nsum(summand, [0, mpmath_inf], method='richardson')
            vf = Float(v, ndig)
            if vold is not None and vold == vf:
                break
            prec += prec  # double precision each time
            vold = vf

        return v._mpf_


def evalf_prod(expr, prec, options):
    from ..concrete import Sum
    if all((l[1] - l[2]).is_Integer for l in expr.limits):
        re, im, re_acc, im_acc = evalf(expr.doit(), prec=prec, options=options)
    else:
        re, im, re_acc, im_acc = evalf(expr.rewrite(Sum), prec=prec, options=options)
    return re, im, re_acc, im_acc


def evalf_sum(expr, prec, options):
    from .numbers import Float, oo
    if 'subs' in options:
        expr = expr.subs(options['subs'])
    func = expr.function
    limits = expr.limits
    if len(limits) != 1 or len(limits[0]) != 3:
        raise NotImplementedError
    if func is S.Zero:
        return mpf(0), None, None, None
    prec2 = prec + 10
    try:
        n, a, b = limits[0]
        if b != oo or a != int(a):
            raise NotImplementedError
        # Use fast hypergeometric summation if possible
        v = hypsum(func, n, int(a), prec2)
        delta = prec - fastlog(v)
        if fastlog(v) < -10:
            v = hypsum(func, n, int(a), delta)
        return v, None, min(prec, delta), None
    except NotImplementedError:
        # Euler-Maclaurin summation for general series
        m, err, eps = prec, oo, Float(2.0)**(-prec)
        while err > eps:
            m <<= 1
            s, err = expr.euler_maclaurin(m=m, n=m, eps=eps,
                                          eval_integral=False)
            err = err.evalf(strict=False)
        err = fastlog(evalf(abs(err), 20, options)[0])
        re, im, re_acc, im_acc = evalf(s, prec2, options)
        if re_acc is None:
            re_acc = -err
        if im_acc is None:
            im_acc = -err
        return re, im, re_acc, im_acc


############################################################################
#                                                                          #
#                            Symbolic interface                            #
#                                                                          #
############################################################################

def evalf_symbol(x, prec, options):
    val = options['subs'][x]
    if isinstance(val, mpf):
        if not val:
            return None, None, None, None
        return val._mpf_, None, prec, None
    else:
        if '_cache' not in options:
            options['_cache'] = {}
        cache = options['_cache']
        cached, cached_prec = cache.get(x, (None, MINUS_INF))
        if cached_prec >= prec:
            return cached
        v = evalf(sympify(val), prec, options)
        cache[x] = (v, prec)
        return v


evalf_table = None


def _create_evalf_table():
    global evalf_table
    from ..concrete.products import Product
    from ..concrete.summations import Sum
    from ..functions.combinatorial.numbers import bernoulli
    from ..functions.elementary.complexes import Abs, im, re
    from ..functions.elementary.exponential import log
    from ..functions.elementary.piecewise import Piecewise
    from ..functions.elementary.trigonometric import atan, cos, sin
    from ..integrals.integrals import Integral
    from .add import Add
    from .mul import Mul
    from .numbers import (Exp1, Float, Half, ImaginaryUnit, Integer, NaN,
                          NegativeOne, One, Pi, Rational, Zero)
    from .power import Pow
    from .symbol import Dummy, Symbol
    evalf_table = {
        Symbol: evalf_symbol,
        Dummy: evalf_symbol,
        Float: lambda x, prec, options: (x._mpf_, None, prec if prec <= x._prec else x._prec, None),
        Rational: lambda x, prec, options: (from_rational(x.numerator, x.denominator, prec),
                                            None, prec, None),
        Integer: lambda x, prec, options: (from_int(x.numerator, prec),
                                           None, prec, None),
        Zero: lambda x, prec, options: (None, None, prec, None),
        One: lambda x, prec, options: (fone, None, prec, None),
        Half: lambda x, prec, options: (fhalf, None, prec, None),
        Pi: lambda x, prec, options: (mpf_pi(prec), None, prec, None),
        Exp1: lambda x, prec, options: (mpf_e(prec), None, prec, None),
        ImaginaryUnit: lambda x, prec, options: (None, fone, None, prec),
        NegativeOne: lambda x, prec, options: (fnone, None, prec, None),
        NaN: lambda x, prec, options: (fnan, None, prec, None),

        cos: evalf_trig,
        sin: evalf_trig,

        Add: evalf_add,
        Mul: evalf_mul,
        Pow: evalf_pow,

        log: evalf_log,
        atan: evalf_atan,
        Abs: evalf_abs,

        re: evalf_re,
        im: evalf_im,

        Integral: evalf_integral,
        Sum: evalf_sum,
        Product: evalf_prod,
        Piecewise: evalf_piecewise,

        bernoulli: evalf_bernoulli,
    }


def evalf(x, prec, options):
    from ..functions import im as im_
    from ..functions import re as re_
    try:
        rf = evalf_table[x.func]
        r = rf(x, prec, options)
    except KeyError:
        try:
            # Fall back to ordinary evalf if possible
            if 'subs' in options:
                x = x.subs(evalf_subs(prec, options['subs']))
            re, im = x._eval_evalf(prec).as_real_imag()
            if re.has(re_) or im.has(im_):
                raise NotImplementedError
            if re == 0:
                re = None
                reprec = None
            else:
                re = re._to_mpmath(prec)._mpf_
                reprec = prec
            if im == 0:
                im = None
                imprec = None
            else:
                im = im._to_mpmath(prec)._mpf_
                imprec = prec
            r = re, im, reprec, imprec
        except AttributeError:
            raise NotImplementedError
    chop = options.get('chop', False)
    if chop:
        if chop is True:
            chop_prec = prec
        else:
            # convert (approximately) from given tolerance;
            # the formula here will will make 1e-i rounds to 0 for
            # i in the range +/-27 while 2e-i will not be chopped
            chop_prec = round(-3.321*math.log10(chop) + 2.5)
        r = chop_parts(r, chop_prec)
    if options.get('strict'):
        check_target(x, r, prec)
    return r


class EvalfMixin:
    """Mixin class adding evalf capability."""

    def evalf(self, dps=15, subs=None, maxn=110, chop=False, strict=True, quad=None):
        """
        Evaluate the given formula to an accuracy of dps decimal digits.
        Optional keyword arguments:

            subs=<dict>
                Substitute numerical values for symbols, e.g.
                subs={x:3, y:1+pi}. The substitutions must be given as a
                dictionary.

            maxn=<integer>
                Allow a maximum temporary working precision of maxn digits
                (default=110)

            chop=<bool>
                Replace tiny real or imaginary parts in subresults
                by exact zeros (default=False)

            strict=<bool>
                Raise PrecisionExhausted if any subresult fails to evaluate
                to full accuracy, given the available maxprec
                (default=True)

            quad=<str>
                Choose algorithm for numerical quadrature. By default,
                tanh-sinh quadrature is used. For oscillatory
                integrals on an infinite interval, try quad='osc'.

        """
        from .numbers import Float, I

        if subs and is_sequence(subs):
            raise TypeError('subs must be given as a dictionary')

        if not evalf_table:
            _create_evalf_table()
        prec = dps_to_prec(dps)
        options = {'maxprec': max(prec, int(maxn*LG10)), 'chop': chop,
                   'strict': strict}
        if subs is not None:
            options['subs'] = subs
        if quad is not None:
            options['quad'] = quad
        try:
            result = evalf(self, prec + 4, options)
        except PrecisionExhausted:
            if self.is_Float and self._prec >= prec:
                return Float._new(self._mpf_, prec)
            else:
                raise
        except NotImplementedError:
            # Fall back to the ordinary evalf
            v = self._eval_evalf(prec)
            if v is None:
                return self
            else:
                # Normalize result
                return v.subs({_: _.evalf(dps, strict=strict)
                               for _ in v.atoms(Float)})
        re, im, re_acc, im_acc = result
        if re:
            p = max(min(prec, re_acc), 1)
            re = Float._new(re, p)
        else:
            re = S.Zero
        if im:
            p = max(min(prec, im_acc), 1)
            im = Float._new(im, p)
            return re + im*I
        else:
            return re

    def _evalf(self, prec):
        """Helper for evalf. Does the same thing but takes binary precision."""
        r = self._eval_evalf(prec)
        if r is None:
            r = self
        return r

    def _eval_evalf(self, prec):
        return

    def _to_mpmath(self, prec):
        # mpmath functions accept ints as input
        errmsg = 'cannot convert to mpmath number'
        if hasattr(self, '_as_mpf_val'):
            return make_mpf(self._as_mpf_val(prec))
        try:
            re, im, _, _ = evalf(self, prec, {'maxprec': DEFAULT_MAXPREC})
            if im:
                if not re:
                    re = fzero
                return make_mpc((re, im))
            elif re:
                return make_mpf(re)
            else:
                return make_mpf(fzero)
        except NotImplementedError:
            v = self._eval_evalf(prec)
            if v is None:
                raise ValueError(errmsg)
            re, im = v.as_real_imag()
            if re.is_Float:
                re = re._mpf_
            else:
                raise ValueError(errmsg)
            if im.is_Float:
                im = im._mpf_
            else:
                raise ValueError(errmsg)
            return make_mpc((re, im))


def N(x, dps=15, **options):
    r"""
    Calls x.evalf(dps, \*\*options).

    Examples
    ========

    >>> Sum(1/k**k, (k, 1, oo))
    Sum(k**(-k), (k, 1, oo))
    >>> N(_, 4)
    1.291

    See Also
    ========

    diofant.core.evalf.EvalfMixin.evalf

    """
    return sympify(x).evalf(dps, **options)
