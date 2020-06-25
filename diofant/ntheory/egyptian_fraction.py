import math

from ..core import Integer


def egyptian_fraction(r, algorithm='Greedy'):
    """Compute an Egyptian fraction of the rational `r`.

    Returns
    =======

    list
        The list of denominators of an Egyptian
        fraction expansion.

    Parameters
    ==========

    r : Rational
        a positive rational number.
    algorithm : { "Greedy", "Graham Jewett", "Takenouchi", "Golomb" }, optional
        Denotes the algorithm to be used (the default is "Greedy").

    Examples
    ========

    >>> egyptian_fraction(Rational(3, 7))
    [3, 11, 231]
    >>> egyptian_fraction(Rational(3, 7), 'Graham Jewett')
    [7, 8, 9, 56, 57, 72, 3192]
    >>> egyptian_fraction(Rational(3, 7), 'Takenouchi')
    [4, 7, 28]
    >>> egyptian_fraction(Rational(3, 7), 'Golomb')
    [3, 15, 35]
    >>> egyptian_fraction(Rational(11, 5), 'Golomb')
    [1, 2, 3, 4, 9, 234, 1118, 2580]

    See Also
    ========

    diofant.core.numbers.Rational

    Notes
    =====

    Currently the following algorithms are supported:

    1) Greedy Algorithm

       Also called the Fibonacci-Sylvester algorithm.
       At each step, extract the largest unit fraction less
       than the target and replace the target with the remainder.

       It has some distinct properties:

       a) Given `p/q` in lowest terms, generates an expansion of maximum
          length `p`. Even as the numerators get large, the number of
          terms is seldom more than a handful.

       b) Uses minimal memory.

       c) The terms can blow up (standard examples of this are 5/121 and
          31/311).  The denominator is at most squared at each step
          (doubly-exponential growth) and typically exhibits
          singly-exponential growth.

    2) Graham Jewett Algorithm

       The algorithm suggested by the result of Graham and Jewett.
       Note that this has a tendency to blow up: the length of the
       resulting expansion is always ``2**(x/gcd(x, y)) - 1``.

    3) Takenouchi Algorithm

       The algorithm suggested by Takenouchi (1921).
       Differs from the Graham-Jewett algorithm only in the handling
       of duplicates.

    4) Golomb's Algorithm

       A method given by Golumb (1962), using modular arithmetic and
       inverses.  It yields the same results as a method using continued
       fractions proposed by Bleicher (1972).

    If the given rational is greater than or equal to 1, a greedy algorithm
    of summing the harmonic sequence 1/1 + 1/2 + 1/3 + ... is used, taking
    all the unit fractions of this sequence until adding one more would be
    greater than the given number.  This list of denominators is prefixed
    to the result from the requested algorithm used on the remainder.  For
    example, if r is 8/3, using the Greedy algorithm, we get [1, 2, 3, 4,
    5, 6, 7, 14, 420], where the beginning of the sequence, [1, 2, 3, 4, 5,
    6, 7] is part of the harmonic sequence summing to 363/140, leaving a
    remainder of 31/420, which yields [14, 420] by the Greedy algorithm.
    The result of egyptian_fraction(Rational(8, 3), "Golomb") is [1, 2, 3,
    4, 5, 6, 7, 14, 574, 2788, 6460, 11590, 33062, 113820], and so on.

    References
    ==========

    * https://en.wikipedia.org/wiki/Egyptian_fraction
    * https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions
    * https://www.ics.uci.edu/~eppstein/numth/egypt/conflict.html
    * https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_42_from129to134.pdf

    """
    if r <= 0:
        raise ValueError('Value must be positive')

    prefix, rem = egypt_harmonic(r)
    if rem == 0:
        return prefix
    x, y = rem.as_numer_denom()

    if algorithm == 'Greedy':
        return prefix + egypt_greedy(x, y)
    elif algorithm == 'Graham Jewett':
        return prefix + egypt_graham_jewett(x, y)
    elif algorithm == 'Takenouchi':
        return prefix + egypt_takenouchi(x, y)
    elif algorithm == 'Golomb':
        return prefix + egypt_golomb(x, y)
    else:
        raise ValueError('Entered invalid algorithm')


def egypt_greedy(x, y):
    if x == 1:
        return [y]
    else:
        a = (-y) % x
        b = y*(y//x + 1)
        c = math.gcd(a, b)
        if c > 1:
            num, denom = a//c, b//c
        else:
            num, denom = a, b
        return [y//x + 1] + egypt_greedy(num, denom)


def egypt_graham_jewett(x, y):
    l = [y] * x

    # l is now a list of integers whose reciprocals sum to x/y.
    # we shall now proceed to manipulate the elements of l without
    # changing the reciprocated sum until all elements are unique.

    while len(l) != len(set(l)):
        l.sort()  # so the list has duplicates. find a smallest pair
        for i in range(len(l) - 1):  # pragma: no branch
            if l[i] == l[i + 1]:
                break
        # we have now identified a pair of identical
        # elements: l[i] and l[i + 1].
        # now comes the application of the result of graham and jewett:
        l[i + 1] = l[i] + 1
        # and we just iterate that until the list has no duplicates.
        l.append(l[i]*(l[i] + 1))
    return sorted(l)


def egypt_takenouchi(x, y):
    l = [y] * x
    while len(l) != len(set(l)):
        l.sort()
        for i in range(len(l) - 1):  # pragma: no branch
            if l[i] == l[i + 1]:
                break
        k = l[i]
        if k % 2 == 0:
            l[i] = l[i] // 2
            del l[i + 1]
        else:
            l[i], l[i + 1] = (k + 1)//2, k*(k + 1)//2
    return sorted(l)


def egypt_golomb(x, y):
    from ..domains import ZZ
    if x == 1:
        return [y]
    xp = ZZ.invert(int(x), int(y))
    rv = [Integer(xp*y)]
    rv.extend(egypt_golomb((x*xp - 1)//y, xp))
    return sorted(rv)


def egypt_harmonic(r):
    rv = []
    d = Integer(1)
    acc = Integer(0)
    while acc + 1/d <= r:
        acc += 1/d
        rv.append(d)
        d += 1
    return rv, r - acc
