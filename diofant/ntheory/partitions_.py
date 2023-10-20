__all__ = 'npartitions',


def npartitions(n):
    """Give the number P(n) of unrestricted partitions of the integer n.

    The partition function P(n) computes the number of ways that n can be
    written as a sum of positive integers, where the order of addends is not
    considered significant.

    Examples
    ========

    >>> npartitions(25)
    1958

    References
    ==========

    * https://mathworld.wolfram.com/PartitionFunctionP.html

    """
    if n < 0:
        return 0
    p = [1] + [0] * n
    for i in range(1, n + 1):
        k = 1
        while True:
            k2 = 3*k**2
            s = 1 if k & 1 else -1
            j = (k2 - k) // 2
            if j > i:
                break
            p[i] += p[i - j]*s
            j = (k2 + k) // 2
            if j > i:
                break
            p[i] += p[i - j]*s
            k += 1
    return p[n]
