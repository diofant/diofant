"""sympify -- convert objects Diofant internal format"""

from inspect import getmro

from .compatibility import iterable
from .evaluate import global_evaluate


class SympifyError(ValueError):
    """Generic sympification error."""

    def __init__(self, expr, base_exc=None):
        self.expr = expr
        self.base_exc = base_exc

    def __str__(self):
        try:
            s = str(self.expr)
        except TypeError:
            s = repr(self.expr)

        return ("Sympify of expression '%s' failed, because of exception being "
                "raised:\n%s: %s" % (s, self.base_exc.__class__.__name__,
                                     str(self.base_exc)))


converter = {}  # See sympify docstring.


class CantSympify:
    """
    Mix in this trait to a class to disallow sympification of its instances.

    Examples
    ========

    >>> class Something(dict):
    ...     pass
    ...
    >>> sympify(Something())
    {}

    >>> class Something(dict, CantSympify):
    ...     pass
    ...
    >>> sympify(Something())
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: {}

    """

    pass


def sympify(a, locals=None, convert_xor=True, strict=False, rational=False,
            evaluate=None):
    """Converts an arbitrary expression to a type that can be used inside Diofant.

    For example, it will convert Python ints into instance of diofant.Rational,
    floats into instances of diofant.Float, etc. It is also able to coerce symbolic
    expressions which inherit from Basic. This can be useful in cooperation
    with SAGE.

    It currently accepts as arguments:
       - any object defined in diofant
       - standard numeric python types: int, long, float, Decimal
       - booleans, including ``None`` (will leave ``None`` unchanged)
       - lists, sets or tuples containing any of the above

    If the argument is already a type that Diofant understands, it will do
    nothing but return that value. This can be used at the beginning of a
    function to ensure you are working with the correct type.

    >>> sympify(2).is_integer
    True
    >>> sympify(2).is_real
    True

    >>> sympify(2.0).is_real
    True

    *Locals*

    The sympification happens with access to everything that is loaded
    by ``from diofant import *``; anything used in a string that is not
    defined by that import will be converted to a symbol. In the following,
    the ``bitcount`` function is treated as a symbol and the ``O`` is
    interpreted as the Order object (used with series) and it raises
    an error when used improperly:

    >>> s = 'bitcount(42)'
    >>> sympify(s)
    bitcount(42)
    >>> sympify("O(x)")
    O(x)

    In order to have the ``O`` interpreted as a Symbol, identify it as such
    in the namespace dictionary. This can be done in a variety of ways; all
    three of the following are possibilities:

    >>> ns = {}
    >>> ns["O"] = Symbol("O")  # method 1
    >>> exec('from diofant.abc import O', ns)  # method 2
    >>> ns.update({O: Symbol("O")})  # method 3
    >>> sympify("O + 1", locals=ns)
    O + 1

    If you want *all* single-letter and Greek-letter variables to be symbols
    then you can use the clashing-symbols dictionaries that have been defined
    there as private variables: _clash1 (single-letter variables), _clash2
    (the multi-letter Greek names) or _clash (both single and multi-letter
    names that are defined in abc).

    >>> from diofant.abc import _clash1
    >>> _clash1
    {'E': E, 'I': I, 'N': N, 'O': O, 'S': S}
    >>> sympify('E & O', _clash1)
    E & O

    *Strict*

    If the option ``strict`` is set to ``True``, only the types for which an
    explicit conversion has been defined are converted. In the other
    cases, a SympifyError is raised.

    >>> print(sympify(None))
    None
    >>> sympify(None, strict=True)
    Traceback (most recent call last):
    ...
    SympifyError: SympifyError: None

    *Extending*

    To extend ``sympify`` to convert custom objects (not derived from ``Basic``),
    just define a ``_diofant_`` method to your class. You can do that even to
    classes that you do not own by subclassing or adding the method at runtime.

    >>> class MyList1:
    ...     def __iter__(self):
    ...         yield 1
    ...         yield 2
    ...         return
    ...
    ...     def __getitem__(self, i):
    ...         return list(self)[i]
    ...
    ...     def _diofant_(self):
    ...         return Matrix(self)
    >>> sympify(MyList1())
    Matrix([
    [1],
    [2]])

    If you do not have control over the class definition you could also use the
    ``converter`` global dictionary. The key is the class and the value is a
    function that takes a single argument and returns the desired Diofant
    object, e.g. ``converter[MyList] = lambda x: Matrix(x)``.

    >>> class MyList2:   # XXX Do not do this if you control the class!
    ...     def __iter__(self):  # Use _diofant_!
    ...         yield 1
    ...         yield 2
    ...         return
    ...
    ...     def __getitem__(self, i):
    ...         return list(self)[i]
    >>> converter[MyList2] = lambda x: Matrix(x)
    >>> sympify(MyList2())
    Matrix([
    [1],
    [2]])

    """
    from .basic import Basic
    from .symbol import Symbol

    if evaluate is None:
        evaluate = global_evaluate[0]
    try:
        if issubclass(a, Basic):
            return a
    except TypeError:  # Type of a is unhashable
        pass
    cls = a.__class__
    if issubclass(cls, Basic):
        return a
    if issubclass(cls, type(None)):
        if strict:
            raise SympifyError(a)
        else:
            return a

    try:
        return converter[cls](a)
    except KeyError:
        for superclass in getmro(cls):
            try:
                return converter[superclass](a)
            except KeyError:
                continue

    if isinstance(a, CantSympify):
        raise SympifyError(a)

    try:
        return a._diofant_()
    except AttributeError:
        pass

    if not isinstance(a, str):
        for coerce in (float, int):
            try:
                return sympify(coerce(a))
            except (TypeError, ValueError, AttributeError, SympifyError):
                continue

    if strict:
        raise SympifyError(a)

    if iterable(a):
        return type(a)([sympify(x, locals=locals, convert_xor=convert_xor,
                                rational=rational) for x in a])
    if isinstance(a, dict):
        return type(a)([sympify(x, locals=locals, convert_xor=convert_xor,
                                rational=rational) for x in a.items()])

    if isinstance(a, str):
        return Symbol(a)
    else:
        raise SympifyError(a)
