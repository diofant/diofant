import math

from ...core import Function, Integer
from ...utilities import default_sort_key, has_dups


###############################################################################
# #################### Kronecker Delta, Levi-Civita etc. #################### #
###############################################################################


def Eijk(*args, **kwargs):
    """
    Represent the Levi-Civita symbol.

    This is just compatibility wrapper to ``LeviCivita()``.

    See Also
    ========

    diofant.functions.special.tensor_functions.LeviCivita

    """
    return LeviCivita(*args, **kwargs)


def eval_levicivita(*args):
    """Evaluate Levi-Civita symbol."""
    from .. import factorial
    n = len(args)
    return math.prod(
        math.prod(args[j] - args[i] for j in range(i + 1, n))
        / factorial(i) for i in range(n))
    # converting factorial(i) to int is slightly faster


class LeviCivita(Function):
    """Represent the Levi-Civita symbol.

    For even permutations of indices it returns 1, for odd permutations -1, and
    for everything else (a repeated index) it returns 0.

    Thus it represents an alternating pseudotensor.

    Examples
    ========

    >>> from diofant.abc import i, j
    >>> LeviCivita(1, 2, 3)
    1
    >>> LeviCivita(1, 3, 2)
    -1
    >>> LeviCivita(1, 2, 2)
    0
    >>> LeviCivita(i, j, k)
    LeviCivita(i, j, k)
    >>> LeviCivita(i, j, i)
    0

    See Also
    ========

    diofant.functions.special.tensor_functions.Eijk

    """

    is_integer = True

    @classmethod
    def eval(cls, *args):
        if all(isinstance(a, (int, Integer)) for a in args):
            return eval_levicivita(*args)
        if has_dups(args):
            return Integer(0)

    def doit(self, **hints):
        return eval_levicivita(*self.args)


class KroneckerDelta(Function):
    """The discrete, or Kronecker, delta function.

    A function that takes in two integers `i` and `j`. It returns `0` if `i` and `j` are
    not equal or it returns `1` if `i` and `j` are equal.

    Parameters
    ==========

    i : Number, Symbol
        The first index of the delta function.
    j : Number, Symbol
        The second index of the delta function.

    Examples
    ========

    A simple example with integer indices::

        >>> KroneckerDelta(1, 2)
        0
        >>> KroneckerDelta(3, 3)
        1

    Symbolic indices::

        >>> from diofant.abc import i, j
        >>> KroneckerDelta(i, j)
        KroneckerDelta(i, j)
        >>> KroneckerDelta(i, i)
        1
        >>> KroneckerDelta(i, i + 1)
        0
        >>> KroneckerDelta(i, i + 1 + k)
        KroneckerDelta(i, i + k + 1)

    See Also
    ========

    diofant.functions.special.tensor_functions.KroneckerDelta.eval
    diofant.functions.special.delta_functions.DiracDelta

    References
    ==========

    * https://en.wikipedia.org/wiki/Kronecker_delta

    """

    is_integer = True

    @classmethod
    def eval(cls, i, j):
        """
        Evaluates the discrete delta function.

        Examples
        ========

        >>> from diofant.abc import i, j

        >>> KroneckerDelta(i, j)
        KroneckerDelta(i, j)
        >>> KroneckerDelta(i, i)
        1
        >>> KroneckerDelta(i, i + 1)
        0
        >>> KroneckerDelta(i, i + 1 + k)
        KroneckerDelta(i, i + k + 1)

        """
        diff = i - j
        if diff.is_zero:
            return Integer(1)
        if diff.is_nonzero:
            return Integer(0)

        if i._assumptions.get('below_fermi') and \
                j._assumptions.get('above_fermi'):
            return Integer(0)
        if j._assumptions.get('below_fermi') and \
                i._assumptions.get('above_fermi'):
            return Integer(0)
        # to make KroneckerDelta canonical
        # following lines will check if inputs are in order
        # if not, will return KroneckerDelta with correct order
        if i is not min(i, j, key=default_sort_key):
            return cls(j, i)

    def _eval_power(self, other):
        if other.is_positive:
            return self
        if other.is_negative and -other != 1:
            return 1/self

    @property
    def is_above_fermi(self):
        """
        True if Delta can be non-zero above fermi

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> p = Symbol('p')
        >>> q = Symbol('q')
        >>> KroneckerDelta(p, a).is_above_fermi
        True
        >>> KroneckerDelta(p, i).is_above_fermi
        False
        >>> KroneckerDelta(p, q).is_above_fermi
        True

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.is_below_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_below_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_above_fermi

        """
        if self.args[0]._assumptions.get('below_fermi'):
            return False
        if self.args[1]._assumptions.get('below_fermi'):
            return False
        return True

    @property
    def is_below_fermi(self):
        """
        True if Delta can be non-zero below fermi

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> p = Symbol('p')
        >>> q = Symbol('q')
        >>> KroneckerDelta(p, a).is_below_fermi
        False
        >>> KroneckerDelta(p, i).is_below_fermi
        True
        >>> KroneckerDelta(p, q).is_below_fermi
        True

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.is_above_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_above_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_below_fermi

        """
        if self.args[0]._assumptions.get('above_fermi'):
            return False
        if self.args[1]._assumptions.get('above_fermi'):
            return False
        return True

    @property
    def is_only_above_fermi(self):
        """
        True if Delta is restricted to above fermi

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> p = Symbol('p')
        >>> q = Symbol('q')
        >>> KroneckerDelta(p, a).is_only_above_fermi
        True
        >>> KroneckerDelta(p, q).is_only_above_fermi
        False
        >>> KroneckerDelta(p, i).is_only_above_fermi
        False

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.is_above_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_below_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_below_fermi

        """
        return (self.args[0]._assumptions.get('above_fermi') or
                self.args[1]._assumptions.get('above_fermi') or False)

    @property
    def is_only_below_fermi(self):
        """
        True if Delta is restricted to below fermi

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> p = Symbol('p')
        >>> q = Symbol('q')
        >>> KroneckerDelta(p, i).is_only_below_fermi
        True
        >>> KroneckerDelta(p, q).is_only_below_fermi
        False
        >>> KroneckerDelta(p, a).is_only_below_fermi
        False

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.is_above_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_below_fermi
        diofant.functions.special.tensor_functions.KroneckerDelta.is_only_above_fermi

        """
        return (self.args[0]._assumptions.get('below_fermi') or
                self.args[1]._assumptions.get('below_fermi') or False)

    @property
    def indices_contain_equal_information(self):
        """
        Returns True if indices are either both above or below fermi.

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> p = Symbol('p')
        >>> q = Symbol('q')
        >>> KroneckerDelta(p, q).indices_contain_equal_information
        True
        >>> KroneckerDelta(p, q+1).indices_contain_equal_information
        True
        >>> KroneckerDelta(i, p).indices_contain_equal_information
        False

        """
        if (self.args[0]._assumptions.get('below_fermi') and
                self.args[1]._assumptions.get('below_fermi')):
            return True
        if (self.args[0]._assumptions.get('above_fermi')
                and self.args[1]._assumptions.get('above_fermi')):
            return True

        # if both indices are general we are True, else false
        return self.is_below_fermi and self.is_above_fermi

    @property
    def preferred_index(self):
        """
        Returns the index which is preferred to keep in the final expression.

        The preferred index is the index with more information regarding fermi
        level.  If indices contain same information, 'a' is preferred before
        'b'.

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> j = Symbol('j', below_fermi=True)
        >>> p = Symbol('p')
        >>> KroneckerDelta(p, i).preferred_index
        i
        >>> KroneckerDelta(p, a).preferred_index
        a
        >>> KroneckerDelta(i, j).preferred_index
        i

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.killable_index

        """
        if self._get_preferred_index():
            return self.args[1]
        return self.args[0]

    @property
    def killable_index(self):
        """
        Returns the index which is preferred to substitute in the final
        expression.

        The index to substitute is the index with less information regarding
        fermi level.  If indices contain same information, 'a' is preferred
        before 'b'.

        Examples
        ========

        >>> a = Symbol('a', above_fermi=True)
        >>> i = Symbol('i', below_fermi=True)
        >>> j = Symbol('j', below_fermi=True)
        >>> p = Symbol('p')
        >>> KroneckerDelta(p, i).killable_index
        p
        >>> KroneckerDelta(p, a).killable_index
        p
        >>> KroneckerDelta(i, j).killable_index
        j

        See Also
        ========

        diofant.functions.special.tensor_functions.KroneckerDelta.preferred_index

        """
        if self._get_preferred_index():
            return self.args[0]
        return self.args[1]

    def _get_preferred_index(self):
        """
        Returns the index which is preferred to keep in the final expression.

        The preferred index is the index with more information regarding fermi
        level.  If indices contain same information, index 0 is returned.

        """
        if not self.is_above_fermi:
            if self.args[0]._assumptions.get('below_fermi'):
                return 0
            return 1
        if not self.is_below_fermi:
            if self.args[0]._assumptions.get('above_fermi'):
                return 0
            return 1
        return 0

    @staticmethod
    def _latex_no_arg(printer):
        return r'\delta'
