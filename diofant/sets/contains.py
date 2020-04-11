"""Implementation of the membership testing."""

from ..logic.boolalg import BooleanFunction


class Contains(BooleanFunction):
    """Asserts that x is an element of the set S

    Examples
    ========

    >>> Contains(Integer(2), S.Integers)
    true
    >>> Contains(Integer(-2), S.Naturals)
    false
    >>> Contains(n, S.Naturals)
    Contains(n, Naturals())

    References
    ==========

    * https://en.wikipedia.org/wiki/Element_%28mathematics%29

    """

    @classmethod
    def eval(cls, x, S):
        ret = S.contains(x)
        if not isinstance(ret, Contains):
            return ret
