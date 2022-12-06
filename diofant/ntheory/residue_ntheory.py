import collections
import math
import random

from ..core import Function, Integer
from ..core.compatibility import as_int
from ..core.numbers import igcdex, mod_inverse
from ..utilities.iterables import cantor_product
from .factor_ import factorint, multiplicity, totient, trailing
from .modular import crt, crt1, crt2
from .primetest import isprime


def n_order(a, n):
    """Returns the order of ``a`` modulo ``n``.

    The order of ``a`` modulo ``n`` is the smallest integer
    ``k`` such that ``a**k`` leaves a remainder of 1 with ``n``.

    Examples
    ========

    >>> n_order(3, 7)
    6
    >>> n_order(4, 7)
    3

    """
    a, n = as_int(a), as_int(n)
    if math.gcd(a, n) != 1:
        raise ValueError('The two numbers should be relatively prime')
    factors = collections.defaultdict(int)
    f = factorint(n)
    for px, kx in f.items():
        if kx > 1:
            factors[px] += kx - 1
        fpx = factorint(px - 1)
        for py, ky in fpx.items():
            factors[py] += ky
    group_order = 1
    for px, kx in factors.items():
        group_order *= px**kx
    order = 1
    if a > n:
        a = a % n
    for p, e in factors.items():
        exponent = group_order
        for f in range(e + 1):
            if pow(a, exponent, n) != 1:
                order *= p ** (e - f + 1)
                break
            exponent = exponent // p
    return order


def _primitive_root_prime_iter(p):
    """Generates the primitive roots for a prime ``p``

    References
    ==========

    * W. Stein "Elementary Number Theory" (2011), page 44

    Examples
    ========

    >>> list(_primitive_root_prime_iter(19))
    [2, 3, 10, 13, 14, 15]

    """
    p = as_int(p)
    v = [(p - 1) // i for i in factorint(p - 1)]
    a = 2
    while a < p:
        for pw in v:
            if pow(a, pw, p) == 1:
                break
        else:
            yield a
        a += 1


def primitive_root(p):
    """Returns the smallest primitive root or None.

    References
    ==========

    * W. Stein "Elementary Number Theory" (2011), page 44
    * P. Hackman "Elementary Number Theory" (2009),  Chapter C

    Parameters
    ==========

    p : positive integer

    Examples
    ========

    >>> primitive_root(19)
    2

    """
    p = as_int(p)
    if p < 1:
        raise ValueError('p is required to be positive')
    if p <= 2:
        return 1
    f = factorint(p)
    if len(f) > 2:
        return
    if len(f) == 2:
        if 2 not in f or f[2] > 1:
            return

        # case p = 2*p1**k, p1 prime
        for p1 in f:  # pragma: no branch
            if p1 != 2:
                break
        i = 1
        while i < p:  # pragma: no branch
            i += 2
            if i % p1 == 0:
                continue
            if is_primitive_root(i, p):
                return i

    else:
        if 2 in f:
            if p == 4:
                return 3
            return
        p1, n = list(f.items())[0]
        if n > 1:
            # see Ref [2], page 81
            g = primitive_root(p1)
            if is_primitive_root(g, p1**2):
                return g
            else:
                for i in range(2, g + p1 + 1):  # pragma: no branch
                    if math.gcd(i, p) == 1 and is_primitive_root(i, p):
                        return i

    return next(_primitive_root_prime_iter(p))


def is_primitive_root(a, p):
    """Returns True if ``a`` is a primitive root of ``p``

    ``a`` is said to be the primitive root of ``p`` if gcd(a, p) == 1 and
    totient(p) is the smallest positive number s.t.::

        a**totient(p) cong 1 mod(p)

    Examples
    ========

    >>> is_primitive_root(3, 10)
    True
    >>> is_primitive_root(9, 10)
    False
    >>> n_order(3, 10) == totient(10)
    True
    >>> n_order(9, 10) == totient(10)
    False

    """
    a, p = as_int(a), as_int(p)
    if math.gcd(a, p) != 1:
        raise ValueError('The two numbers should be relatively prime')
    if a > p:
        a = a % p
    return n_order(a, p) == totient(p)


def _sqrt_mod_tonelli_shanks(a, p):
    """Returns the square root in the case of ``p`` prime with ``p == 1 (mod 8)``

    References
    ==========

    * R. Crandall and C. Pomerance "Prime Numbers", 2nt Ed., page 101

    """
    s = trailing(p - 1)
    t = p >> s
    # find a non-quadratic residue
    while 1:
        d = random.randint(2, p - 1)
        r = legendre_symbol(d, p)
        if r == -1:
            break
    # assert legendre_symbol(d, p) == -1
    A = pow(a, t, p)
    D = pow(d, t, p)
    m = 0
    for i in range(s):
        adm = A*pow(D, m, p) % p
        adm = pow(adm, 2**(s - 1 - i), p)
        if adm % p == p - 1:
            m += 2**i
    # assert A*pow(D, m, p) % p == 1
    x = pow(a, (t + 1)//2, p)*pow(D, m//2, p) % p
    return x


def sqrt_mod(a, p, all_roots=False):
    """Find a root of ``x**2 = a mod p``.

    Parameters
    ==========

    a : integer
    p : positive integer
    all_roots : if True the list of roots is returned or None

    Notes
    =====

    If there is no root it is returned None; else the returned root
    is less or equal to ``p // 2``; in general is not the smallest one.
    It is returned ``p // 2`` only if it is the only root.

    Use ``all_roots`` only when it is expected that all the roots fit
    in memory; otherwise use ``sqrt_mod_iter``.

    Examples
    ========

    >>> sqrt_mod(11, 43)
    21
    >>> sqrt_mod(17, 32, True)
    [7, 9, 23, 25]

    """
    if all_roots:
        return sorted(sqrt_mod_iter(a, p))
    try:
        p = abs(as_int(p))
        it = sqrt_mod_iter(a, p)
        r = next(it)
        if r > p // 2:
            return p - r
        elif r < p // 2:
            return r
        else:
            try:
                r = next(it)
                if r > p // 2:
                    return p - r
            except StopIteration:
                pass
            return r
    except StopIteration:
        return


def sqrt_mod_iter(a, p, domain=int):
    """Iterate over solutions to ``x**2 = a mod p``.

    Parameters
    ==========

    a : integer
    p : positive integer
    domain : integer domain, ``int``, ``ZZ`` or ``Integer``

    Examples
    ========

    >>> list(sqrt_mod_iter(11, 43))
    [21, 22]

    """
    from ..domains import ZZ
    a, p = as_int(a), abs(as_int(p))
    if isprime(p):
        a = a % p
        if a == 0:
            res = _sqrt_mod1(a, p, 1)
        else:
            res = _sqrt_mod_prime_power(a, p, 1)
        if res:
            if domain is ZZ:
                for x in res:
                    yield x
            else:
                for x in res:
                    yield domain(x)
    else:
        f = factorint(p)
        v = []
        pv = []
        for px, ex in f.items():
            if a % px == 0:
                rx = _sqrt_mod1(a, px, ex)
                if not rx:
                    return
            else:
                rx = _sqrt_mod_prime_power(a, px, ex)
                if not rx:
                    return
            v.append(rx)
            pv.append(px**ex)
        mm, e, s = crt1(pv)
        if domain is ZZ:
            for vx in cantor_product(*v):
                r = crt2(pv, vx, mm, e, s)[0]
                yield r
        else:
            for vx in cantor_product(*v):
                r = crt2(pv, vx, mm, e, s)[0]
                yield domain(r)


def _sqrt_mod_prime_power(a, p, k):
    """Find the solutions to ``x**2 = a mod p**k`` when ``a % p != 0``.

    Parameters
    ==========

    a : integer
    p : prime number
    k : positive integer

    References
    ==========

    * P. Hackman "Elementary Number Theory" (2009),  page 160
    * http://www.numbertheory.org/php/squareroot.html
    * :cite:`Gathen1999modern`

    Examples
    ========

    >>> _sqrt_mod_prime_power(11, 43, 1)
    [21, 22]

    """
    from ..domains import ZZ

    assert k > 0

    pk = p**k
    a = a % pk

    if k == 1:
        if p == 2:
            return [ZZ(a)]
        if not is_quad_residue(a, p):
            return

        if p % 4 == 3:
            res = pow(a, (p + 1) // 4, p)
        elif p % 8 == 5:
            sign = pow(a, (p - 1) // 4, p)
            if sign == 1:
                res = pow(a, (p + 3) // 8, p)
            else:
                b = pow(4*a, (p - 5) // 8, p)
                x = (2*a*b) % p
                assert pow(x, 2, p) == a
                res = x
        else:
            res = _sqrt_mod_tonelli_shanks(a, p)

        # ``_sqrt_mod_tonelli_shanks(a, p)`` is not deterministic;
        # sort to get always the same result
        return sorted([ZZ(res), ZZ(p - res)])

    # see Ref.[2]
    if p == 2:
        if a % 8 != 1:
            return
        if k <= 3:
            s = set()
            for i in range(0, pk, 4):
                s.add(1 + i)
                s.add(-1 + i)
            return list(s)
        # according to Ref.[2] for k > 2 there are two solutions
        # (mod 2**k-1), that is four solutions (mod 2**k), which can be
        # obtained from the roots of x**2 = 0 (mod 8)
        rv = [ZZ(1), ZZ(3), ZZ(5), ZZ(7)]
        # hensel lift them to solutions of x**2 = 0 (mod 2**k)
        # if r**2 - a = 0 mod 2**nx but not mod 2**(nx+1)
        # then r + 2**(nx - 1) is a root mod 2**(nx+1)
        n = 3
        res = []
        for r in rv:
            nx = n
            while nx < k:
                r1 = (r**2 - a) >> nx
                if r1 % 2:
                    r = r + (1 << (nx - 1))
                assert (r**2 - a) % (1 << (nx + 1)) == 0
                nx += 1
            if r not in res:
                res.append(r)
            x = r + (1 << (k - 1))
            assert (x**2 - a) % pk == 0
            if x < (1 << nx) and x not in res:
                res.append(x)
        return res
    rv = _sqrt_mod_prime_power(a, p, 1)
    if not rv:
        return
    r = rv[0]
    fr = r**2 - a
    # hensel lifting with Newton iteration, see Ref.[3] chapter 9
    # with f(x) = x**2 - a; one has f'(a) != 0 (mod p) for p != 2
    n = 1
    px = p
    while 1:
        n1 = n
        n1 *= 2
        if n1 > k:
            break
        n = n1
        px = px**2
        frinv = igcdex(2*r, px)[0]
        r = (r - fr*frinv) % px
        fr = r**2 - a
    if n < k:
        px = p**k
        frinv = igcdex(2*r, px)[0]
        r = (r - fr*frinv) % px
    return [r, px - r]


def _sqrt_mod1(a, p, n):
    """Find solution to ``x**2 == a mod p**n`` when ``a % p == 0``

    References
    ==========

    * http://www.numbertheory.org/php/squareroot.html

    """
    pn = p**n
    a = a % pn
    if a == 0:
        # case gcd(a, p**k) = p**n
        m = n // 2
        if n % 2 == 1:
            pm1 = p**(m + 1)

            def _iter0a():
                i = 0
                while i < pn:
                    yield i
                    i += pm1
            return _iter0a()
        else:
            pm = p**m

            def _iter0b():
                i = 0
                while i < pn:
                    yield i
                    i += pm
            return _iter0b()

    # case gcd(a, p**k) = p**r, r < n
    f = factorint(a)
    r = f[p]
    if r % 2 == 1:
        return
    m = r // 2
    a1 = a >> r
    if p == 2:
        if n - r == 1:
            pnm1 = 1 << (n - m + 1)
            pm1 = 1 << (m + 1)

            def _iter1():
                k = 1 << (m + 2)
                i = 1 << m
                while i < pnm1:
                    j = i
                    while j < pn:
                        yield j
                        j += k
                    i += pm1
            return _iter1()
        if n - r == 2:
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return
            pnm = 1 << (n - m)

            def _iter2():
                s = set()
                for r in res:
                    i = 0
                    while i < pn:
                        x = (r << m) + i
                        assert x not in s
                        s.add(x)
                        yield x
                        i += pnm
            return _iter2()
        else:  # n - r > 2
            res = _sqrt_mod_prime_power(a1, p, n - r)
            if res is None:
                return
            pnm1 = 1 << (n - m - 1)

            def _iter3():
                s = set()
                for r in res:
                    i = 0
                    while i < pn:
                        x = ((r << m) + i) % pn
                        if x not in s:
                            s.add(x)
                            yield x
                        i += pnm1
            return _iter3()
    else:
        m = r // 2
        a1 = a // p**r
        res1 = _sqrt_mod_prime_power(a1, p, n - r)
        if res1 is None:
            return
        pm = p**m
        pnr = p**(n-r)
        pnm = p**(n-m)

        def _iter4():
            s = set()
            pm = p**m
            for rx in res1:
                i = 0
                while i < pnm:
                    x = ((rx + i) % pn)
                    assert x not in s
                    s.add(x)
                    yield x*pm
                    i += pnr

        return _iter4()


def is_quad_residue(a, p):
    """
    Returns True if ``a`` (mod ``p``) is in the set of squares mod ``p``,
    i.e a % p in {i**2 % p for i in range(p)}. If ``p`` is an odd
    prime, an iterative method is used to make the determination:

    >>> sorted({i**2 % 7 for i in range(7)})
    [0, 1, 2, 4]
    >>> [j for j in range(7) if is_quad_residue(j, 7)]
    [0, 1, 2, 4]

    See Also
    ========

    legendre_symbol, jacobi_symbol

    """
    a, p = as_int(a), as_int(p)
    if p < 1:
        raise ValueError('p must be > 0')
    if a >= p or a < 0:
        a = a % p
    if a < 2 or p < 3:
        return True
    if not isprime(p):
        if p % 2 and jacobi_symbol(a, p) == -1:
            return False
        r = sqrt_mod(a, p)
        if r is None:
            return False
        else:
            return True

    return pow(a, (p - 1) // 2, p) == 1


def is_nthpow_residue(a, n, m):
    """Returns True if ``x**n == a (mod m)`` has solutions.

    References
    ==========

    * P. Hackman "Elementary Number Theory" (2009),  page 76

    """
    a, n, m = map(as_int, (a, n, m))
    if m <= 0:
        raise ValueError('m must be > 0')
    if n < 0:
        raise ValueError('n must be >= 0')
    if a < 0:
        raise ValueError('a must be >= 0')
    if n == 0:
        if m == 1:
            return False
        return a == 1
    if n == 1:
        return True
    if n == 2:
        return is_quad_residue(a, m)
    return _is_nthpow_residue_bign(a, n, m)


def _is_nthpow_residue_bign(a, n, m):
    """Returns True if ``x**n == a (mod m)`` has solutions for n > 2."""
    assert n > 2
    assert a >= 0
    assert m > 0
    if primitive_root(m) is None:
        assert m >= 8
        for prime, power in factorint(m).items():
            if not _is_nthpow_residue_bign_prime_power(a, n, prime, power):
                return False
        return True
    f = totient(m)
    k = f // math.gcd(f, n)
    return pow(a, k, m) == 1


def _is_nthpow_residue_bign_prime_power(a, n, p, k):
    """Returns True/False if a solution for ``x**n == a (mod(p**k))``
    does/doesn't exist.

    """
    assert a >= 0
    assert n > 2
    assert isprime(p)
    assert k > 0
    if a % p:
        if p != 2:
            return _is_nthpow_residue_bign(a, n, pow(p, k))
        if n & 1:
            return True
        c = trailing(n)
        return a % pow(2, min(c + 2, k)) == 1
    else:
        a %= pow(p, k)
        if not a:
            return True
        mu = multiplicity(p, a)
        if mu % n:
            return False
        pm = pow(p, mu)
        return _is_nthpow_residue_bign_prime_power(a//pm, n, p, k - mu)


def _nthroot_mod2(s, q, p):
    f = factorint(q)
    v = []
    for b, e in f.items():
        v.extend([b]*e)
    for qx in v:
        s = _nthroot_mod1(s, qx, p, False)
    return s


def _nthroot_mod1(s, q, p, all_roots):
    """Root of ``x**q = s mod p``, ``p`` prime and ``q`` divides ``p - 1``.

    References
    ==========

    * A. M. Johnston "A Generalized qth Root Algorithm"

    """
    g = primitive_root(p)
    if not isprime(q):
        r = _nthroot_mod2(s, q, p)
    else:
        f = p - 1
        assert (p - 1) % q == 0
        # determine k
        k = 0
        while f % q == 0:
            k += 1
            f = f // q
        # find z, x, r1
        f1 = igcdex(-f, q)[0] % q
        z = f*f1
        x = (1 + z) // q
        r1 = pow(s, x, p)
        s1 = pow(s, f, p)
        h = pow(g, f*q, p)
        t = discrete_log(p, s1, h)
        g2 = pow(g, z*t, p)
        g3 = igcdex(g2, p)[0]
        r = r1*g3 % p
        assert pow(r, q, p) == s
    res = [r]
    h = pow(g, (p - 1) // q, p)
    assert pow(h, q, p) == 1
    hx = r
    for _ in range(q - 1):
        hx = (hx*h) % p
        res.append(hx)
    if all_roots:
        res.sort()
        return res
    return min(res)


def nthroot_mod(a, n, p, all_roots=False):
    """Find the solutions to ``x**n = a mod p``.

    Parameters
    ==========

    a : integer
    n : positive integer
    p : positive integer
    all_roots : if False returns the smallest root, else the list of roots

    Examples
    ========

    >>> nthroot_mod(11, 4, 19)
    8
    >>> nthroot_mod(11, 4, 19, True)
    [8, 11]
    >>> nthroot_mod(68, 3, 109)
    23

    """
    if n == 2:
        return sqrt_mod(a, p, all_roots)
    # see Hackman "Elementary Number Theory" (2009), page 76
    if not is_nthpow_residue(a, n, p):
        return
    if primitive_root(p) is None:
        raise NotImplementedError('Not Implemented for m without primitive root')

    if (p - 1) % n == 0:
        return _nthroot_mod1(a, n, p, all_roots)
    # The roots of ``x**n - a = 0 (mod p)`` are roots of
    # ``gcd(x**n - a, x**(p - 1) - 1) = 0 (mod p)``
    pa = n
    pb = p - 1
    b = 1
    if pa < pb:
        a, pa, b, pb = b, pb, a, pa
    while pb:
        # x**pa - a = 0; x**pb - b = 0
        # x**pa - a = x**(q*pb + r) - a = (x**pb)**q * x**r - a =
        #             b**q * x**r - a; x**r - c = 0; c = b**-q * a mod p
        q, r = divmod(pa, pb)
        c = pow(b, q, p)
        c = igcdex(c, p)[0]
        c = (c * a) % p
        pa, pb = pb, r
        a, b = b, c
    if pa == 1:
        if all_roots:
            res = [a]
        else:
            res = a
    elif pa == 2:
        return sqrt_mod(a, p, all_roots)
    else:
        res = _nthroot_mod1(a, pa, p, all_roots)
    return res


def quadratic_residues(p):
    """Returns the list of quadratic residues.

    Examples
    ========

    >>> quadratic_residues(7)
    [0, 1, 2, 4]

    """
    r = set()
    for i in range(p // 2 + 1):
        r.add(pow(i, 2, p))
    return sorted(r)


def legendre_symbol(a, p):
    r"""Returns the Legendre symbol `(a / p)`.

    For an integer ``a`` and an odd prime ``p``, the Legendre symbol is
    defined as

    .. math ::
        \genfrac(){}{}{a}{p} = \begin{cases}
             0 & \text{if } p \text{ divides } a\\
             1 & \text{if } a \text{ is a quadratic residue modulo } p\\
            -1 & \text{if } a \text{ is a quadratic nonresidue modulo } p
        \end{cases}

    Parameters
    ==========

    a : integer
    p : odd prime

    Examples
    ========

    >>> [legendre_symbol(i, 7) for i in range(7)]
    [0, 1, 1, -1, 1, -1, -1]
    >>> sorted({i**2 % 7 for i in range(7)})
    [0, 1, 2, 4]

    See Also
    ========

    is_quad_residue, jacobi_symbol

    References
    ==========

    * https://en.wikipedia.org/wiki/Legendre_symbol

    """
    a, p = as_int(a), as_int(p)
    if not isprime(p) or p == 2:
        raise ValueError('p should be an odd prime')
    a = a % p
    if not a:
        return 0
    if is_quad_residue(a, p):
        return 1
    return -1


def jacobi_symbol(m, n):
    r"""
    Returns the Jacobi symbol `(m / n)`.

    For any integer ``m`` and any positive odd integer ``n`` the Jacobi symbol
    is defined as the product of the Legendre symbols corresponding to the
    prime factors of ``n``:

    .. math ::
        \genfrac(){}{}{m}{n} =
            \genfrac(){}{}{m}{p^{1}}^{\alpha_1}
            \genfrac(){}{}{m}{p^{2}}^{\alpha_2}
            ...
            \genfrac(){}{}{m}{p^{k}}^{\alpha_k}
            \text{ where } n =
                p_1^{\alpha_1}
                p_2^{\alpha_2}
                ...
                p_k^{\alpha_k}

    Like the Legendre symbol, if the Jacobi symbol `\genfrac(){}{}{m}{n} = -1`
    then ``m`` is a quadratic nonresidue modulo ``n``.

    But, unlike the Legendre symbol, if the Jacobi symbol
    `\genfrac(){}{}{m}{n} = 1` then ``m`` may or may not be a quadratic residue
    modulo ``n``.

    Parameters
    ==========

    m : integer
    n : odd positive integer

    Examples
    ========

    >>> jacobi_symbol(45, 77)
    -1
    >>> jacobi_symbol(60, 121)
    1

    The relationship between the ``jacobi_symbol`` and ``legendre_symbol`` can
    be demonstrated as follows:

    >>> L = legendre_symbol
    >>> Integer(45).factors()
    {3: 2, 5: 1}
    >>> jacobi_symbol(7, 45) == L(7, 3)**2 * L(7, 5)**1
    True

    See Also
    ========

    is_quad_residue, legendre_symbol

    """
    m, n = as_int(m), as_int(n)
    if not n % 2:
        raise ValueError('n should be an odd integer')
    if m < 0 or m > n:
        m = m % n
    if not m:
        return int(n == 1)
    if n == 1 or m == 1:
        return 1
    if math.gcd(m, n) != 1:
        return 0

    j = 1
    s = trailing(m)
    m = m >> s
    if s % 2 and n % 8 in [3, 5]:
        j *= -1

    while m != 1:
        if m % 4 == 3 and n % 4 == 3:
            j *= -1
        m, n = n % m, m
        s = trailing(m)
        m = m >> s
        if s % 2 and n % 8 in [3, 5]:
            j *= -1
    return j


class mobius(Function):
    """Möbius function maps natural number to {-1, 0, 1}

    It is defined as follows:
        1) `1` if `n = 1`.
        2) `0` if `n` has a squared prime factor.
        3) `(-1)^k` if `n` is a square-free positive integer with `k`
           number of prime factors.

    It is an important multiplicative function in number theory
    and combinatorics.  It has applications in mathematical series,
    algebraic number theory and also physics (Fermion operator has very
    concrete realization with Möbius Function model).

    Parameters
    ==========

    n : positive integer

    Examples
    ========

    >>> mobius(13*7)
    1
    >>> mobius(1)
    1
    >>> mobius(13*7*5)
    -1
    >>> mobius(13**2)
    0

    References
    ==========

    * https://en.wikipedia.org/wiki/M%C3%B6bius_function
    * Thomas Koshy "Elementary Number Theory with Applications"

    """

    @classmethod
    def eval(cls, n):
        if n.is_integer:
            if n.is_positive is not True:
                raise ValueError('n should be a positive integer')
        else:
            raise TypeError('n should be an integer')
        if n.is_prime:
            return Integer(-1)
        elif n == 1:
            return Integer(1)
        elif n.is_Integer:
            a = factorint(n)
            if any(i > 1 for i in a.values()):
                return Integer(0)
            return Integer(-1)**len(a)


def _discrete_log_trial_mul(n, a, b, order=None):
    """
    Trial multiplication algorithm for computing the discrete logarithm of
    ``a`` to the base ``b`` modulo ``n``.

    The algorithm finds the discrete logarithm using exhaustive search. This
    naive method is used as fallback algorithm of ``discrete_log`` when the
    group order is very small.

    """
    a %= n
    b %= n
    if order is None:
        order = n
    x = 1
    for i in range(order):
        if x == a:
            return i
        x = x * b % n
    raise ValueError('Log does not exist')


def _discrete_log_shanks_steps(n, a, b, order=None):
    """
    Baby-step giant-step algorithm for computing the discrete logarithm of
    ``a`` to the base ``b`` modulo ``n``.

    The algorithm is a time-memory trade-off of the method of exhaustive
    search. It uses `O(sqrt(m))` memory, where `m` is the group order.

    """
    a %= n
    b %= n
    if order is None:
        order = n_order(b, n)
    m = math.isqrt(order) + 1
    T = {}
    x = 1
    for i in range(m):
        T[x] = i
        x = x * b % n
    z = mod_inverse(b, n)
    z = pow(z, m, n)
    x = a
    for i in range(m):
        if x in T:
            return i * m + T[x]
        x = x * z % n
    raise ValueError('Log does not exist')


def _discrete_log_pohlig_hellman(n, a, b, order=None):
    """
    Pohlig-Hellman algorithm for computing the discrete logarithm of ``a`` to
    the base ``b`` modulo ``n``.

    In order to compute the discrete logarithm, the algorithm takes advantage
    of the factorization of the group order. It is more efficient when the
    group order factors into many small primes.

    """
    a %= n
    b %= n

    if order is None:
        order = n_order(b, n)

    f = factorint(order)
    l = [0] * len(f)

    for i, (pi, ri) in enumerate(f.items()):
        for j in range(ri):
            gj = pow(b, l[i], n)
            aj = pow(a * mod_inverse(gj, n), order // pi**(j + 1), n)
            bj = pow(b, order // pi, n)
            cj = discrete_log(n, aj, bj, pi, True)
            l[i] += cj * pi**j

    d, _ = crt([pi**ri for pi, ri in f.items()], l)
    return d


def discrete_log(n, a, b, order=None, prime_order=None):
    """
    Compute the discrete logarithm of ``a`` to the base ``b`` modulo ``n``.

    This is a recursive function to reduce the discrete logarithm problem in
    cyclic groups of composite order to the problem in cyclic groups of
    prime order.

    Notes
    =====

    It employs different algorithms depending on the problem (subgroup order
    size, prime order or not):

    * Trial multiplication
    * Baby-step giant-step
    * Pohlig-Hellman

    References
    ==========

    * https://mathworld.wolfram.com/DiscreteLogarithm.html
    * :cite:`Menezes97`

    Examples
    ========

    >>> discrete_log(41, 15, 7)
    3

    """
    if order is None:
        order = n_order(b, n)

    if prime_order is None:
        prime_order = isprime(order)

    if order < 1000:
        return _discrete_log_trial_mul(n, a, b, order)
    elif prime_order:
        return _discrete_log_shanks_steps(n, a, b, order)

    return _discrete_log_pohlig_hellman(n, a, b, order)
