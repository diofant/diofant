"""
Integer factorization
"""

import math
import numbers
import random

from ..core import (Function, Integer, Mul, Pow, Rational, integer_nthroot,
                    sympify)
from ..core.compatibility import as_int
from ..core.evalf import bitcount
from .generate import nextprime, primerange, sieve
from .primetest import isprime


small_trailing = [i and max(int(not i % 2**j) and j for j in range(1, 8))
                  for i in range(256)]


def smoothness(n):
    """Return the B-smooth and B-power smooth values of n.

    The smoothness of n is the largest prime factor of n; the power-
    smoothness is the largest divisor raised to its multiplicity.

    >>> smoothness(2**7*3**2)
    (3, 128)
    >>> smoothness(2**4*13)
    (13, 16)
    >>> smoothness(2)
    (2, 2)

    See Also
    ========

    factorint, smoothness_p

    """

    if n == 1:
        return (1, 1)  # not prime, but otherwise this causes headaches
    facs = factorint(n)
    return max(facs), max(m**facs[m] for m in facs)


def smoothness_p(n, m=-1, power=0, visual=None):
    r"""
    Return a list of [m, (p, (M, sm(p + m), psm(p + m)))...]
    where:

    1. p**M is the base-p divisor of n
    2. sm(p + m) is the smoothness of p + m (m = -1 by default)
    3. psm(p + m) is the power smoothness of p + m

    The list is sorted according to smoothness (default) or by power smoothness
    if power=1.

    The smoothness of the numbers to the left (m = -1) or right (m = 1) of a
    factor govern the results that are obtained from the p +/- 1 type factoring
    methods.

        >>> smoothness_p(10431, m=1)
        (1, [(3, (2, 2, 4)), (19, (1, 5, 5)), (61, (1, 31, 31))])
        >>> smoothness_p(10431)
        (-1, [(3, (2, 2, 2)), (19, (1, 3, 9)), (61, (1, 5, 5))])
        >>> smoothness_p(10431, power=1)
        (-1, [(3, (2, 2, 2)), (61, (1, 5, 5)), (19, (1, 3, 9))])

    If visual=True then an annotated string will be returned:

        >>> print(smoothness_p(21477639576571, visual=1))
        p**i=4410317**1 has p-1 B=1787, B-pow=1787
        p**i=4869863**1 has p-1 B=2434931, B-pow=2434931

    This string can also be generated directly from a factorization dictionary
    and vice versa:

        >>> f = factorint(17*9)
        >>> f
        {3: 2, 17: 1}
        >>> smoothness_p(f)
        'p**i=3**2 has p-1 B=2, B-pow=2\np**i=17**1 has p-1 B=2, B-pow=16'
        >>> smoothness_p(_)
        {3: 2, 17: 1}

    The table of the output logic is:

        ====== ====== ======= =======
        |              Visual
        ------ ----------------------
        Input  True   False   other
        ====== ====== ======= =======
        dict    str    tuple   str
        str     str    tuple   dict
        tuple   str    tuple   str
        n       str    tuple   tuple
        mul     str    tuple   tuple
        ====== ====== ======= =======

    See Also
    ========

    factorint, smoothness

    """
    from ..utilities import flatten

    # visual must be True, False or other (stored as None)
    if visual in (1, 0):
        visual = bool(visual)
    if visual not in (True, False):
        visual = None

    if type(n) is str:
        if visual:
            return n
        d = {}
        for li in n.splitlines():
            k, v = [int(i) for i in
                    li.split('has')[0].split('=')[1].split('**')]
            d[k] = v
        if visual is not True and visual is not False:
            return d
        return smoothness_p(d, visual=False)
    elif type(n) is not tuple:
        facs = factorint(n, visual=False)

    if power:
        k = -1
    else:
        k = 1
    if type(n) is not tuple:
        rv = (m, sorted(((f, tuple([M] + list(smoothness(f + m))))
                         for f, M in [i for i in facs.items()]),
                        key=lambda x: (x[1][k], x[0])))
    else:
        rv = n

    if visual is False or (visual is not True) and (type(n) in [int, Mul]):
        return rv
    lines = []
    for dat in rv[1]:
        dat = flatten(dat)
        dat.insert(2, m)
        lines.append('p**i=%i**%i has p%+i B=%i, B-pow=%i' % tuple(dat))
    return '\n'.join(lines)


def trailing(n):
    """Count the number of trailing zero digits in the binary
    representation of n, i.e. determine the largest power of 2
    that divides n.

    Examples
    ========

    >>> trailing(128)
    7
    >>> trailing(63)
    0

    """
    n = int(n)
    if not n:
        return 0
    low_byte = n & 0xff
    if low_byte:
        return small_trailing[low_byte]

    # 2**m is quick for z up through 2**30
    z = bitcount(n) - 1
    if n == 1 << z:
        return z

    t = 0
    p = 8
    while not n & 1:
        while not n & ((1 << p) - 1):
            n >>= p
            t += p
            p *= 2
        p //= 2
    return t


def multiplicity(p, n):
    """Find the greatest integer m such that p**m divides n.

    Examples
    ========

    >>> [multiplicity(5, n) for n in [8, 5, 25, 125, 250]]
    [0, 1, 2, 3, 3]
    >>> multiplicity(3, Rational(1, 9))
    -2

    """
    try:
        p, n = as_int(p), as_int(n)
    except ValueError:
        if all(isinstance(i, (numbers.Integral, Rational)) for i in (p, n)):
            p, n = Rational(p), Rational(n)
            if p.denominator == 1:
                if n.numerator == 1:
                    return -multiplicity(p.numerator, n.denominator)
                return Integer(0)
            elif p.numerator == 1:
                return multiplicity(p.denominator, n.denominator)
            else:
                like = min(multiplicity(p.numerator, n.numerator), multiplicity(p.denominator, n.denominator))
                cross = min(multiplicity(p.denominator, n.numerator), multiplicity(p.numerator, n.denominator))
                return like - cross
        raise ValueError('expecting ints or fractions, got %s and %s' % (p, n))

    if n == 0:
        raise ValueError('multiplicity of 0 is not defined')
    if p == 2:
        return trailing(n)
    if p < 2:
        raise ValueError('p must be an integer, 2 or larger, but got %s' % p)
    if p == n:
        return 1

    m = 0
    n, rem = divmod(n, p)
    while not rem:
        m += 1
        if m > 5:
            # The multiplicity could be very large. Better
            # to increment in powers of two
            e = 2
            while 1:
                ppow = p**e
                if ppow < n:
                    nnew, rem = divmod(n, ppow)
                    if not rem:
                        m += e
                        e *= 2
                        n = nnew
                        continue
                return m + multiplicity(p, n)
        n, rem = divmod(n, p)
    return m


def perfect_power(n, candidates=None, big=True, factor=True):
    """
    Return ``(b, e)`` such that ``n`` == ``b**e`` if ``n`` is a
    perfect power; otherwise return ``False``.

    By default, the base is recursively decomposed and the exponents
    collected so the largest possible ``e`` is sought. If ``big=False``
    then the smallest possible ``e`` (thus prime) will be chosen.

    If ``candidates`` for exponents are given, they are assumed to be sorted
    and the first one that is larger than the computed maximum will signal
    failure for the routine.

    If ``factor=True`` then simultaneous factorization of n is attempted
    since finding a factor indicates the only possible root for n. This
    is True by default since only a few small factors will be tested in
    the course of searching for the perfect power.

    Examples
    ========

    >>> perfect_power(16)
    (2, 4)
    >>> perfect_power(16, big = False)
    (4, 2)

    """
    n = int(n)
    if n < 3:
        return False
    logn = math.log(n, 2)
    max_possible = int(logn) + 2  # only check values less than this
    not_square = n % 10 in [2, 3, 7, 8]  # squares cannot end in 2, 3, 7, 8
    if not candidates:
        candidates = primerange(2 + not_square, max_possible)

    afactor = 2 + n % 2
    for e in candidates:
        if e < 3:
            if e == 1 or e == 2 and not_square:
                continue
        if e > max_possible:
            return False

        # see if there is a factor present
        if factor:
            if n % afactor == 0:
                # find what the potential power is
                if afactor == 2:
                    e = trailing(n)
                else:
                    e = multiplicity(afactor, n)
                # if it's a trivial power we are done
                if e == 1:
                    return False

                # maybe the bth root of n is exact
                r, exact = integer_nthroot(n, e)
                if not exact:
                    # then remove this factor and check to see if
                    # any of e's factors are a common exponent; if
                    # not then it's not a perfect power
                    n //= afactor**e
                    m = perfect_power(n, candidates=primefactors(e), big=big)
                    if m is False:
                        return False
                    else:
                        r, m = m
                        # adjust the two exponents so the bases can
                        # be combined
                        g = math.gcd(m, e)
                        if g == 1:
                            return False
                        m //= g
                        e //= g
                        r, e = r**m*afactor**e, g
                if not big:
                    e0 = primefactors(e)
                    if len(e0) > 1 or e0[0] != e:
                        e0 = e0[0]
                        r, e = r**(e//e0), e0
                return r, e
            else:
                # get the next factor ready for the next pass through the loop
                afactor = nextprime(afactor)

        # Weed out downright impossible candidates
        if logn/e < 40:
            b = 2.0**(logn/e)
            if abs(int(b + 0.5) - b) > 0.01:
                continue

        # now see if the plausible e makes a perfect power
        r, exact = integer_nthroot(n, e)
        if exact:
            if big:
                m = perfect_power(r, big=big, factor=factor)
                if m is not False:
                    r, e = m[0], e*m[1]
            return int(r), e
    return False


def pollard_rho(n, s=2, a=1, retries=5, seed=1234, max_steps=None, F=None):
    r"""
    Use Pollard's rho method to try to extract a nontrivial factor
    of ``n``. The returned factor may be a composite number. If no
    factor is found, ``None`` is returned.

    The algorithm generates pseudo-random values of x with a generator
    function, replacing x with F(x). If F is not supplied then the
    function x**2 + ``a`` is used. The first value supplied to F(x) is ``s``.
    Upon failure (if ``retries`` is > 0) a new ``a`` and ``s`` will be
    supplied; the ``a`` will be ignored if F was supplied.

    The sequence of numbers generated by such functions generally have a
    a lead-up to some number and then loop around back to that number and
    begin to repeat the sequence, e.g. 1, 2, 3, 4, 5, 3, 4, 5 -- this leader
    and loop look a bit like the Greek letter rho, and thus the name, 'rho'.

    For a given function, very different leader-loop values can be obtained
    so it is a good idea to allow for retries:

    >>> n = 16843009
    >>> F = lambda x: (2048*pow(x, 2, n) + 32767)%n
    >>> for s in range(5):
    ...     print('loop length = %4i; leader length = %3i' % next(cycle_length(F, s)))
    ...
    loop length = 2489; leader length =  42
    loop length =   78; leader length = 120
    loop length = 1482; leader length =  99
    loop length = 1482; leader length = 285
    loop length = 1482; leader length = 100

    Here is an explicit example where there is a two element leadup to
    a sequence of 3 numbers (11, 14, 4) that then repeat:

    >>> x = 2
    >>> for i in range(9):
    ...     x = (x**2 + 12)%17
    ...     print(x)
    ...
    16
    13
    11
    14
    4
    11
    14
    4
    11
    >>> next(cycle_length(lambda x: (x**2+12)%17, 2))
    (3, 2)
    >>> list(cycle_length(lambda x: (x**2+12)%17, 2, values=True))
    [16, 13, 11, 14, 4]

    Instead of checking the differences of all generated values for a gcd
    with n, only the kth and 2*kth numbers are checked, e.g. 1st and 2nd,
    2nd and 4th, 3rd and 6th until it has been detected that the loop has been
    traversed. Loops may be many thousands of steps long before rho finds a
    factor or reports failure. If ``max_steps`` is specified, the iteration
    is cancelled with a failure after the specified number of steps.

    Examples
    ========

    >>> n = 16843009
    >>> F = lambda x: (2048*pow(x, 2, n) + 32767) % n
    >>> pollard_rho(n, F=F)
    257

    Use the default setting with a bad value of ``a`` and no retries:

    >>> pollard_rho(n, a=n-2, retries=0)

    If retries is > 0 then perhaps the problem will correct itself when
    new values are generated for a:

    >>> pollard_rho(n, a=n-2, retries=1)
    257

    References
    ==========

    * Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
      A Computational Perspective", Springer, 2nd edition, 229-231

    """
    n = int(n)
    if n < 5:
        raise ValueError('pollard_rho should receive n > 4')
    prng = random.Random(seed + retries)
    V = s
    for i in range(retries + 1):
        U = V
        if not F:
            def F(x):
                return (pow(x, 2, n) + a) % n
        j = 0
        while 1:
            if max_steps and (j > max_steps):
                break
            j += 1
            U = F(U)
            V = F(F(V))  # V is 2x further along than U
            g = math.gcd(U - V, n)
            if g == 1:
                continue
            if g == n:
                break
            return int(g)
        V = prng.randint(0, n - 1)
        a = prng.randint(1, n - 3)  # for x**2 + a, a%n should not be 0 or -2
        F = None
    return


def pollard_pm1(n, B=10, a=2, retries=0, seed=1234):
    """
    Use Pollard's p-1 method to try to extract a nontrivial factor
    of ``n``. Either a divisor (perhaps composite) or ``None`` is returned.

    The value of ``a`` is the base that is used in the test gcd(a**M - 1, n).
    The default is 2.  If ``retries`` > 0 then if no factor is found after the
    first attempt, a new ``a`` will be generated randomly (using the ``seed``)
    and the process repeated.

    Note: the value of M is lcm(1..B) = reduce(ilcm, range(2, B + 1)).

    A search is made for factors next to even numbers having a power smoothness
    less than ``B``. Choosing a larger B increases the likelihood of finding a
    larger factor but takes longer. Whether a factor of n is found or not
    depends on ``a`` and the power smoothness of the even mumber just less than
    the factor p (hence the name p - 1).

    Although some discussion of what constitutes a good ``a`` some
    descriptions are hard to interpret. At the modular.math site referenced
    below it is stated that if gcd(a**M - 1, n) = N then a**M % q**r is 1
    for every prime power divisor of N. But consider the following:

        >>> n = 257*1009
        >>> smoothness_p(n)
        (-1, [(257, (1, 2, 256)), (1009, (1, 7, 16))])

    So we should (and can) find a root with B=16:

        >>> pollard_pm1(n, B=16, a=3)
        1009

    If we attempt to increase B to 256 we find that it doesn't work:

        >>> pollard_pm1(n, B=256)
        >>>

    But if the value of ``a`` is changed we find that only multiples of
    257 work, e.g.:

        >>> pollard_pm1(n, B=256, a=257)
        1009

    Checking different ``a`` values shows that all the ones that didn't
    work had a gcd value not equal to ``n`` but equal to one of the
    factors:

        >>> M = 1
        >>> for i in range(2, 256):
        ...     M = ilcm(M, i)
        ...
        >>> {math.gcd(pow(a, M, n) - 1, n) for a in range(2, 256) if
        ...  math.gcd(pow(a, M, n) - 1, n) != n}
        {1009}

    But does aM % d for every divisor of n give 1?

        >>> aM = pow(255, M, n)
        >>> [(d, aM%Pow(*d.args)) for d in factorint(n, visual=True).args]
        [(257**1, 1), (1009**1, 1)]

    No, only one of them. So perhaps the principle is that a root will
    be found for a given value of B provided that:

    1) the power smoothness of the p - 1 value next to the root
       does not exceed B
    2) a**M % p != 1 for any of the divisors of n.

    By trying more than one ``a`` it is possible that one of them
    will yield a factor.

    Examples
    ========

    With the default smoothness bound, this number can't be cracked:

        >>> pollard_pm1(21477639576571)

    Increasing the smoothness bound helps:

        >>> pollard_pm1(21477639576571, B=2000)
        4410317

    Looking at the smoothness of the factors of this number we find:

        >>> print(smoothness_p(21477639576571, visual=1))
        p**i=4410317**1 has p-1 B=1787, B-pow=1787
        p**i=4869863**1 has p-1 B=2434931, B-pow=2434931

    The B and B-pow are the same for the p - 1 factorizations of the divisors
    because those factorizations had a very large prime factor:

        >>> factorint(4410317 - 1)
        {2: 2, 617: 1, 1787: 1}
        >>> factorint(4869863-1)
        {2: 1, 2434931: 1}

    Note that until B reaches the B-pow value of 1787, the number is not cracked;

        >>> pollard_pm1(21477639576571, B=1786)
        >>> pollard_pm1(21477639576571, B=1787)
        4410317

    The B value has to do with the factors of the number next to the divisor,
    not the divisors themselves. A worst case scenario is that the number next
    to the factor p has a large prime divisisor or is a perfect power. If these
    conditions apply then the power-smoothness will be about p/2 or p. The more
    realistic is that there will be a large prime factor next to p requiring
    a B value on the order of p/2. Although primes may have been searched for
    up to this level, the p/2 is a factor of p - 1, something that we don't
    know. The modular.math reference below states that 15% of numbers in the
    range of 10**15 to 15**15 + 10**4 are 10**6 power smooth so a B of 10**6
    will fail 85% of the time in that range. From 10**8 to 10**8 + 10**3 the
    percentages are nearly reversed...but in that range the simple trial
    division is quite fast.

    References
    ==========

    * Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
      A Computational Perspective", Springer, 2nd edition, 236-238
    * https://web.archive.org/web/20150716201437/http://modular.math.washington.edu/edu/2007/spring/ent/ent-html/node81.html
    * https://web.archive.org/web/20170830055619/http://www.cs.toronto.edu/~yuvalf/Factorization.pdf

    """

    n = int(n)
    if n < 4 or B < 3:
        raise ValueError('pollard_pm1 should receive n > 3 and B > 2')
    prng = random.Random(seed + B)

    # computing a**lcm(1,2,3,..B) % n for B > 2
    # it looks weird, but it's right: primes run [2, B]
    # and the answer's not right until the loop is done.
    for i in range(retries + 1):
        aM = a
        for p in sieve.primerange(2, B + 1):
            e = int(math.log(B, p))
            aM = pow(aM, pow(p, e), n)
        g = math.gcd(aM - 1, n)
        if 1 < g < n:
            return int(g)

        # get a new a:
        # since the exponent, lcm(1..B), is even, if we allow 'a' to be 'n-1'
        # then (n - 1)**even % n will be 1 which will give a g of 0 and 1 will
        # give a zero, too, so we set the range as [2, n-2]. Some references
        # say 'a' should be coprime to n, but either will detect factors.
        a = prng.randint(2, n - 2)


def _trial(factors, n, candidates, verbose=False):
    """
    Helper function for integer factorization. Trial factors ``n`
    against all integers given in the sequence ``candidates``
    and updates the dict ``factors`` in-place. Returns the reduced
    value of ``n`` and a flag indicating whether any factors were found.

    """
    if verbose:
        factors0 = list(factors)
    nfactors = len(factors)
    for d in candidates:
        if n % d == 0:
            m = multiplicity(d, n)
            n //= d**m
            factors[d] = m
    if verbose:
        for k in sorted(set(factors).difference(set(factors0))):
            print(factor_msg % (k, factors[k]))
    return int(n), len(factors) != nfactors


def _check_termination(factors, n, limitp1, use_trial, use_rho, use_pm1,
                       verbose):
    """
    Helper function for integer factorization. Checks if ``n``
    is a prime or a perfect power, and in those cases updates
    the factorization and raises ``StopIteration``.

    """

    if verbose:
        print('Check for termination')

    # since we've already been factoring there is no need to do
    # simultaneous factoring with the power check
    p = perfect_power(n, factor=False)
    if p is not False:
        base, exp = p
        if limitp1:
            limit = limitp1 - 1
        else:
            limit = limitp1
        facs = factorint(base, limit, use_trial, use_rho, use_pm1,
                         verbose=False)
        for b, e in facs.items():
            if verbose:
                print(factor_msg % (b, e))
            factors[b] = exp*e
        return True

    if isprime(n):
        factors[int(n)] = 1
        return True

    if n == 1:
        return True

    return False


trial_int_msg = "Trial division with ints [%i ... %i] and fail_max=%i"
trial_msg = "Trial division with primes [%i ... %i]"
rho_msg = "Pollard's rho with retries %i, max_steps %i and seed %i"
pm1_msg = "Pollard's p-1 with smoothness bound %i and seed %i"
factor_msg = '\t%i ** %i'
fermat_msg = 'Close factors satisying Fermat condition found.'
complete_msg = 'Factorization is complete.'


def _factorint_small(factors, n, limit, fail_max):
    """
    Return the value of n and either a 0 (indicating that factorization up
    to the limit was complete) or else the next near-prime that would have
    been tested.

    Factoring stops if there are fail_max unsuccessful tests in a row.

    If factors of n were found they will be in the factors dictionary as
    {factor: multiplicity} and the returned value of n will have had those
    factors removed. The factors dictionary is modified in-place.

    """

    def done(n, d):
        """return n, d if the sqrt(n) wasn't reached yet, else
        n, 0 indicating that factoring is done.

        """
        if d*d <= n:
            return n, d
        return n, 0

    d = 2
    m = trailing(n)
    if m:
        factors[d] = m
        n >>= m
    d = 3
    if limit < d:
        if n > 1:
            factors[n] = 1
        return done(n, d)
    # reduce
    m = 0
    while n % d == 0:
        n //= d
        m += 1
        if m == 20:
            mm = multiplicity(d, n)
            m += mm
            n //= d**mm
            break
    if m:
        factors[d] = m

    # when d*d exceeds maxx or n we are done; if limit**2 is greater
    # than n then maxx is set to zero so the value of n will flag the finish
    if limit*limit > n:
        maxx = 0
    else:
        maxx = limit*limit

    dd = maxx or n
    d = 5
    fails = 0
    while fails < fail_max:
        if d*d > dd:
            break
        # d = 6*i - 1
        # reduce
        m = 0
        while n % d == 0:
            n //= d
            m += 1
            if m == 20:
                mm = multiplicity(d, n)
                m += mm
                n //= d**mm
                break
        if m:
            factors[d] = m
            dd = maxx or n
            fails = 0
        else:
            fails += 1
        d += 2
        if d*d > dd:
            break
        # d = 6*i - 1
        # reduce
        m = 0
        while n % d == 0:
            n //= d
            m += 1
            if m == 20:
                mm = multiplicity(d, n)
                m += mm
                n //= d**mm
                break
        if m:
            factors[d] = m
            dd = maxx or n
            fails = 0
        else:
            fails += 1
        # d = 6*(i+1) - 1
        d += 4

    return done(n, d)


def factorint(n, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              verbose=False, visual=None):
    r"""
    Given a positive integer ``n``, ``factorint(n)`` returns a dict containing
    the prime factors of ``n`` as keys and their respective multiplicities
    as values. For example:

    >>> factorint(2000)    # 2000 = (2**4) * (5**3)
    {2: 4, 5: 3}
    >>> factorint(65537)   # This number is prime
    {65537: 1}

    For input less than 2, factorint behaves as follows:

        - ``factorint(1)`` returns the empty factorization, ``{}``
        - ``factorint(0)`` returns ``{0:1}``
        - ``factorint(-n)`` adds ``-1:1`` to the factors and then factors ``n``

    Partial Factorization:

    If ``limit`` (> 3) is specified, the search is stopped after performing
    trial division up to (and including) the limit (or taking a
    corresponding number of rho/p-1 steps). This is useful if one has
    a large number and only is interested in finding small factors (if
    any). Note that setting a limit does not prevent larger factors
    from being found early; it simply means that the largest factor may
    be composite. Since checking for perfect power is relatively cheap, it is
    done regardless of the limit setting.

    This number, for example, has two small factors and a huge
    semi-prime factor that cannot be reduced easily:

    >>> a = 1407633717262338957430697921446883
    >>> f = factorint(a, limit=10000)
    >>> f
    {7: 1, 991: 1, 202916782076162456022877024859: 1}
    >>> isprime(max(f))
    False

    This number has a small factor and a residual perfect power whose
    base is greater than the limit:

    >>> factorint(3*101**7, limit=5)
    {3: 1, 101: 7}

    Visual Factorization:

    If ``visual`` is set to ``True``, then it will return a visual
    factorization of the integer.  For example:

    >>> pprint(factorint(4200, visual=True), use_unicode=False)
     3  1  2  1
    2 *3 *5 *7

    Note that this is achieved by using the evaluate=False flag in Mul
    and Pow. If you do other manipulations with an expression where
    evaluate=False, it may evaluate.  Therefore, you should use the
    visual option only for visualization, and use the normal dictionary
    returned by visual=False if you want to perform operations on the
    factors.

    You can easily switch between the two forms by sending them back to
    factorint:

    >>> regular = factorint(1764); regular
    {2: 2, 3: 2, 7: 2}
    >>> pprint(factorint(regular), use_unicode=False)
     2  2  2
    2 *3 *7

    >>> visual = factorint(1764, visual=True); pprint(visual, use_unicode=False)
     2  2  2
    2 *3 *7
    >>> print(factorint(visual))
    {2: 2, 3: 2, 7: 2}

    If you want to send a number to be factored in a partially factored form
    you can do so with a dictionary or unevaluated expression:

    >>> factorint(factorint({4: 2, 12: 3})) # twice to toggle to dict form
    {2: 10, 3: 3}
    >>> factorint(Mul(4, 12, evaluate=False))
    {2: 4, 3: 1}

    The table of the output logic is:

        ====== ====== ======= =======
                       Visual
        ------ ----------------------
        Input  True   False   other
        ====== ====== ======= =======
        dict    mul    dict    mul
        n       mul    dict    dict
        mul     mul    dict    dict
        ====== ====== ======= =======

    Notes
    =====

    The function switches between multiple algorithms. Trial division
    quickly finds small factors (of the order 1-5 digits), and finds
    all large factors if given enough time. The Pollard rho and p-1
    algorithms are used to find large factors ahead of time; they
    will often find factors of the order of 10 digits within a few
    seconds:

    >>> factors = factorint(12345678910111213141516)
    >>> for base, exp in sorted(factors.items()):
    ...     print('%s %s' % (base, exp))
    ...
    2 2
    2507191691 1
    1231026625769 1

    Any of these methods can optionally be disabled with the following
    boolean parameters:

        - ``use_trial``: Toggle use of trial division
        - ``use_rho``: Toggle use of Pollard's rho method
        - ``use_pm1``: Toggle use of Pollard's p-1 method

    ``factorint`` also periodically checks if the remaining part is
    a prime number or a perfect power, and in those cases stops.


    If ``verbose`` is set to ``True``, detailed progress is printed.

    See Also
    ========

    smoothness, smoothness_p, divisors

    """
    factordict = {}
    if visual and not isinstance(n, Mul) and not isinstance(n, dict):
        factordict = factorint(n, limit=limit, use_trial=use_trial,
                               use_rho=use_rho, use_pm1=use_pm1,
                               verbose=verbose, visual=False)
    elif isinstance(n, Mul):
        factordict = {int(k): int(v) for k, v in
                      list(n.as_powers_dict().items())}
    elif isinstance(n, dict):
        factordict = n
    if factordict and (isinstance(n, Mul) or isinstance(n, dict)):
        # check it
        for k in list(factordict):
            if isprime(k):
                continue
            e = factordict.pop(k)
            d = factorint(k, limit=limit, use_trial=use_trial, use_rho=use_rho,
                          use_pm1=use_pm1, verbose=verbose, visual=False)
            for k, v in d.items():
                if k in factordict:
                    factordict[k] += v*e
                else:
                    factordict[k] = v*e
    if visual or (type(n) is dict and
                  visual is not True and
                  visual is not False):
        if factordict == {}:
            return Integer(1)
        if -1 in factordict:
            factordict.pop(-1)
            args = [Integer(-1)]
        else:
            args = []
        args.extend([Pow(*i, evaluate=False)
                     for i in sorted(factordict.items())])
        return Mul(*args, evaluate=False)
    elif isinstance(n, dict) or isinstance(n, Mul):
        return factordict

    assert use_trial or use_rho or use_pm1

    n = as_int(n)
    if limit:
        limit = int(limit)

    # special cases
    if n < 0:
        factors = factorint(
            -n, limit=limit, use_trial=use_trial, use_rho=use_rho,
            use_pm1=use_pm1, verbose=verbose, visual=False)
        factors[-1] = 1
        return factors

    if limit and limit < 2:
        if n == 1:
            return {}
        return {n: 1}
    elif n < 10:
        # doing this we are assured of getting a limit > 2
        # when we have to compute it later
        return [{0: 1}, {}, {2: 1}, {3: 1}, {2: 2}, {5: 1},
                {2: 1, 3: 1}, {7: 1}, {2: 3}, {3: 2}][n]

    factors = {}

    # do simplistic factorization
    if verbose:
        sn = str(n)
        if len(sn) > 50:
            print('Factoring %s' % sn[:5] +
                  '..(%i other digits)..' % (len(sn) - 10) + sn[-5:])
        else:
            print('Factoring', n)

    if use_trial:
        # this is the preliminary factorization for small factors
        small = 2**15
        fail_max = 600
        small = min(small, limit or small)
        if verbose:
            print(trial_int_msg % (2, small, fail_max))
        n, next_p = _factorint_small(factors, n, small, fail_max)
    else:
        next_p = 2
    if factors and verbose:
        for k in sorted(factors):
            print(factor_msg % (k, factors[k]))
    if next_p == 0:
        if n > 1:
            factors[int(n)] = 1
        if verbose:
            print(complete_msg)
        return factors

    # continue with more advanced factorization methods

    # first check if the simplistic run didn't finish
    # because of the limit and check for a perfect
    # power before exiting
    if limit and next_p > limit:
        if verbose:
            print('Exceeded limit:', limit)

        if _check_termination(factors, n, limit, use_trial,
                              use_rho, use_pm1, verbose):
            if verbose:
                print(complete_msg)
            return factors

        factors[int(n)] = 1
        return factors
    else:
        # Before quitting (or continuing on)...

        # ...do a Fermat test since it's so easy and we need the
        # square root anyway. Finding 2 factors is easy if they are
        # "close enough." This is the big root equivalent of dividing by
        # 2, 3, 5.
        sqrt_n = integer_nthroot(n, 2)[0]
        a = sqrt_n + 1
        a2 = a**2
        b2 = a2 - n
        for i in range(3):
            b, fermat = integer_nthroot(b2, 2)
            if fermat:
                break
            b2 += 2*a + 1  # equiv to (a+1)**2 - n
            a += 1
        if fermat:
            if verbose:
                print(fermat_msg)
            if limit:
                limit -= 1
            for r in [a - b, a + b]:
                facs = factorint(r, limit=limit, use_trial=use_trial,
                                 use_rho=use_rho, use_pm1=use_pm1,
                                 verbose=verbose)
                factors.update(facs)
            if verbose:
                print(complete_msg)
            return factors

        # ...see if factorization can be terminated
        if _check_termination(factors, n, limit, use_trial,
                              use_rho, use_pm1, verbose):
            if verbose:
                print(complete_msg)
            return factors

    # these are the limits for trial division which will
    # be attempted in parallel with pollard methods
    low, high = next_p, 2*next_p

    limit = limit or sqrt_n
    # add 1 to make sure limit is reached in primerange calls
    limit += 1

    while 1:

        high_ = high
        if limit < high_:
            high_ = limit

        # Trial division
        if use_trial:
            if verbose:
                print(trial_msg % (low, high_))
            ps = sieve.primerange(low, high_)
            n, found_trial = _trial(factors, n, ps, verbose)
            if found_trial:
                if _check_termination(factors, n, limit, use_trial, use_rho,
                                      use_pm1, verbose):
                    if verbose:
                        print(complete_msg)
                    return factors
        else:
            found_trial = False

        if high > limit:
            if verbose:
                print('Exceeded limit:', limit)
            factors[int(n)] = 1
            if verbose:
                print(complete_msg)
            return factors

        # Only used advanced methods when no small factors were found
        if not found_trial:
            if use_pm1 or use_rho:
                high_root = max(int(math.log(high_**0.7)), low, 3)

                # Pollard p-1
                if use_pm1:
                    if verbose:
                        print(pm1_msg % (high_root, high_))
                    c = pollard_pm1(n, B=high_root, seed=high_)
                    if c:
                        # factor it and let _trial do the update
                        ps = factorint(c, limit=limit - 1,
                                       use_trial=use_trial,
                                       use_rho=use_rho,
                                       use_pm1=use_pm1,
                                       verbose=verbose)
                        n, _ = _trial(factors, n, ps, verbose=False)
                        if _check_termination(factors, n, limit, use_trial,
                                              use_rho, use_pm1, verbose):
                            if verbose:
                                print(complete_msg)
                            return factors

                # Pollard rho
                if use_rho:
                    max_steps = high_root
                    if verbose:
                        print(rho_msg % (1, max_steps, high_))
                    c = pollard_rho(n, retries=1, max_steps=max_steps,
                                    seed=high_)
                    if c:
                        # factor it and let _trial do the update
                        ps = factorint(c, limit=limit - 1,
                                       use_trial=use_trial,
                                       use_rho=use_rho,
                                       use_pm1=use_pm1,
                                       verbose=verbose)
                        n, _ = _trial(factors, n, ps, verbose=False)
                        if _check_termination(factors, n, limit, use_trial,
                                              use_rho, use_pm1, verbose):
                            if verbose:
                                print(complete_msg)
                            return factors

        low, high = high, high*2


def factorrat(rat, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              verbose=False, visual=None):
    r"""
    Given a Rational ``r``, ``factorrat(r)`` returns a dict containing
    the prime factors of ``r`` as keys and their respective multiplicities
    as values. For example:

    >>> factorrat(Rational(8, 9))  # = (2**3) * (3**-2)
    {2: 3, 3: -2}
    >>> factorrat(Rational(-1, 987))  # = -1 * (3**-1) * (7**-1) * (47**-1)
    {-1: 1, 3: -1, 7: -1, 47: -1}

    Please see the docstring for ``factorint`` for detailed explanations
    and examples of the following keywords:

        - ``limit``: Integer limit up to which trial division is done
        - ``use_trial``: Toggle use of trial division
        - ``use_rho``: Toggle use of Pollard's rho method
        - ``use_pm1``: Toggle use of Pollard's p-1 method
        - ``verbose``: Toggle detailed printing of progress
        - ``visual``: Toggle product form of output

    """
    from collections import defaultdict
    f = factorint(rat.numerator, limit=limit, use_trial=use_trial,
                  use_rho=use_rho, use_pm1=use_pm1,
                  verbose=verbose).copy()
    f = defaultdict(int, f)
    for p, e in factorint(rat.denominator, limit=limit,
                          use_trial=use_trial,
                          use_rho=use_rho,
                          use_pm1=use_pm1,
                          verbose=verbose).items():
        f[p] += -e

    if not visual:
        return dict(f)
    else:
        if -1 in f:
            f.pop(-1)
            args = [Integer(-1)]
        else:
            args = []
        args.extend([Pow(*i, evaluate=False)
                     for i in sorted(f.items())])
        return Mul(*args, evaluate=False)


def primefactors(n, limit=None, verbose=False):
    """Return a sorted list of n's prime factors, ignoring multiplicity
    and any composite factor that remains if the limit was set too low
    for complete factorization. Unlike factorint(), primefactors() does
    not return -1 or 0.

    Examples
    ========

    >>> primefactors(6)
    [2, 3]
    >>> primefactors(-5)
    [5]

    >>> sorted(factorint(123456).items())
    [(2, 6), (3, 1), (643, 1)]
    >>> primefactors(123456)
    [2, 3, 643]

    >>> sorted(factorint(10000000001, limit=200).items())
    [(101, 1), (99009901, 1)]
    >>> isprime(99009901)
    False
    >>> primefactors(10000000001, limit=300)
    [101]

    See Also
    ========

    divisors

    """
    n = int(n)
    factors = sorted(factorint(n, limit=limit, verbose=verbose))
    s = [f for f in factors[:-1:] if f not in [-1, 0, 1]]
    if factors and isprime(factors[-1]):
        s += [factors[-1]]
    return s


def _divisors(n):
    """Helper function for divisors which generates the divisors."""

    factordict = factorint(n)
    ps = sorted(factordict)

    def rec_gen(n=0):
        if n == len(ps):
            yield 1
        else:
            pows = [1]
            for j in range(factordict[ps[n]]):
                pows.append(pows[-1] * ps[n])
            for q in rec_gen(n + 1):
                for p in pows:
                    yield p * q

    for p in rec_gen():
        yield p


def divisors(n, generator=False):
    r"""Return all divisors of n.

    Divisors are sorted from 1..n by default.  If generator is True an
    unordered generator is returned.

    The number of divisors of n can be quite large if there are many
    prime factors (counting repeated factors). If only the number of
    factors is desired use divisor_count(n).

    Examples
    ========

    >>> divisors(24)
    [1, 2, 3, 4, 6, 8, 12, 24]
    >>> divisor_count(24)
    8

    >>> list(divisors(120, generator=True))
    [1, 2, 4, 8, 3, 6, 12, 24, 5, 10, 20, 40, 15, 30, 60, 120]

    See Also
    ========

    primefactors, factorint, divisor_count

    References
    ==========

    * https://stackoverflow.com/questions/1010381/python-factorization

    """

    n = as_int(abs(n))
    if isprime(n):
        return [1, n]
    if n == 1:
        return [1]
    if n == 0:
        return []
    rv = _divisors(n)
    if not generator:
        return sorted(rv)
    return rv


def divisor_count(n, modulus=1):
    """Return the number of divisors of ``n``.

    If ``modulus`` is not 1 then only those that are divisible by
    ``modulus`` are counted.

    References
    ==========

    * https://web.archive.org/web/20130629014824/http://www.mayer.dial.pipex.com:80/maths/formulae.htm

    Examples
    ========

    >>> divisor_count(6)
    4

    See Also
    ========

    factorint, divisors, totient

    """

    if not modulus:
        return 0
    elif modulus != 1:
        n, r = divmod(n, modulus)
        if r:
            return 0
    if n == 0:
        return 0
    return Mul(*[v + 1 for k, v in factorint(n).items() if k > 1])


def _antidivisors(n):
    """Helper function for antidivisors which generates the antidivisors."""

    for d in _divisors(n):
        y = 2*d
        if n > y and n % y:
            yield y
    for d in _divisors(2*n-1):
        if n > d >= 2 and n % d:
            yield d
    for d in _divisors(2*n+1):
        if n > d >= 2 and n % d:
            yield d


def antidivisors(n, generator=False):
    r"""Return all antidivisors of n sorted from 1..n by default.

    Antidivisors of n are numbers that do not divide n by the largest
    possible margin.  If generator is True an unordered generator is returned.

    References
    ==========

    * definition is described in http://oeis.org/A066272/a066272a.html

    Examples
    ========

    >>> antidivisors(24)
    [7, 16]

    >>> sorted(antidivisors(128, generator=True))
    [3, 5, 15, 17, 51, 85]

    See Also
    ========

    primefactors, factorint, divisors, divisor_count, antidivisor_count

    """

    n = as_int(abs(n))
    if n <= 2:
        return []
    rv = _antidivisors(n)
    if not generator:
        return sorted(rv)
    return rv


def antidivisor_count(n):
    """Return the number of antidivisors of ``n``.

    References
    ==========

    * formula from https://oeis.org/A066272

    Examples
    ========

    >>> antidivisor_count(13)
    4
    >>> antidivisor_count(27)
    5

    See Also
    ========

    factorint, divisors, antidivisors, divisor_count, totient

    """

    n = as_int(abs(n))
    if n <= 2:
        return 0
    return divisor_count(2*n-1) + divisor_count(2*n+1) + \
        divisor_count(n) - divisor_count(n, 2) - 5


class totient(Function):
    """Calculate the Euler totient function phi(n)

    >>> totient(1)
    1
    >>> totient(25)
    20

    See Also
    ========

    divisor_count

    """

    @classmethod
    def eval(cls, n):
        n = sympify(n)
        if n.is_Integer:
            if n < 1:
                raise ValueError("n must be a positive integer")
            factors = factorint(n)
            t = 1
            for p, k in factors.items():
                t *= (p - 1) * p**(k - 1)
            return t

    def _eval_is_integer(self):
        if self.args[0].is_integer and self.args[0].is_positive:
            return True


class divisor_sigma(Function):
    r"""Calculate the divisor function `\sigma_k(n)` for positive integer n.

    ``divisor_sigma(n, k)`` is equal to ``sum(x**k for x in divisors(n))``

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^\omega p_i^{m_i},

    then

    .. math ::
        \sigma_k(n) = \prod_{i=1}^\omega (1+p_i^k+p_i^{2k}+\cdots
        + p_i^{m_ik}).

    Parameters
    ==========

    k : power of divisors in the sum

        for k = 0, 1:
        ``divisor_sigma(n, 0)`` is equal to ``divisor_count(n)``
        ``divisor_sigma(n, 1)`` is equal to ``sum(divisors(n))``

        Default for k is 1.

    References
    ==========

    * https://en.wikipedia.org/wiki/Divisor_function

    Examples
    ========

    >>> divisor_sigma(18, 0)
    6
    >>> divisor_sigma(39, 1)
    56
    >>> divisor_sigma(12, 2)
    210
    >>> divisor_sigma(37)
    38

    See Also
    ========

    divisor_count, totient, divisors, factorint

    """

    @classmethod
    def eval(cls, n, k=1):
        n = sympify(n)
        k = sympify(k)
        if n.is_prime:
            return 1 + n**k
        if n.is_Integer:
            if n <= 0:
                raise ValueError("n must be a positive integer")
            else:
                return Mul(*[(p**(k*(e + 1)) - 1)/(p**k - 1) if k != 0
                             else e + 1 for p, e in factorint(n).items()])


def core(n, t=2):
    r"""Calculate core(n,t) = `core_t(n)` of a positive integer n.

    ``core_2(n)`` is equal to the squarefree part of n

    If n's prime factorization is:

    .. math ::
        n = \prod_{i=1}^\omega p_i^{m_i},

    then

    .. math ::
        core_t(n) = \prod_{i=1}^\omega p_i^{m_i \mod t}.

    Parameters
    ==========

    t : core(n,t) calculates the t-th power free part of n

        ``core(n, 2)`` is the squarefree part of ``n``
        ``core(n, 3)`` is the cubefree part of ``n``

        Default for t is 2.

    References
    ==========

    * https://en.wikipedia.org/wiki/Square-free_integer#Squarefree_core

    Examples
    ========

    >>> from diofant.ntheory.factor_ import core
    >>> core(24, 2)
    6
    >>> core(9424, 3)
    1178
    >>> core(379238)
    379238
    >>> core(15**11, 10)
    15

    See Also
    ========

    factorint, diofant.solvers.diophantine.square_factor

    """

    n = as_int(n)
    t = as_int(t)
    if n <= 0:
        raise ValueError("n must be a positive integer")
    elif t <= 1:
        raise ValueError("t must be >= 2")
    else:
        y = 1
        for p, e in factorint(n).items():
            y *= p**(e % t)
        return y
