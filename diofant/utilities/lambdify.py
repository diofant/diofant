"""
This module provides convenient functions to transform diofant expressions to
lambda functions which can be used to calculate numerical values very fast.
"""

from __future__ import annotations

import inspect
import textwrap
import typing

from ..core.compatibility import is_sequence, iterable
from ..external import import_module  # noqa: F401
from .decorator import doctest_depends_on


# These are the namespaces the lambda functions will use.
MATH: dict[str, typing.Any] = {}
MPMATH: dict[str, typing.Any] = {}
NUMPY: dict[str, typing.Any] = {}
DIOFANT: dict[str, typing.Any] = {}

# Default namespaces, letting us define translations that can't be defined
# by simple variable maps, like I => 1j
# These are separate from the names above because the above names are modified
# throughout this file, whereas these should remain unmodified.
MATH_DEFAULT: dict[str, typing.Any] = {}
MPMATH_DEFAULT: dict[str, typing.Any] = {}
NUMPY_DEFAULT: dict[str, typing.Any] = {'I': 1j}
DIOFANT_DEFAULT: dict[str, typing.Any] = {}

# Mappings between diofant and other modules function names.
MATH_TRANSLATIONS = {
    'Abs': 'fabs',
    'ceiling': 'ceil',
    'E': 'e',
    'ln': 'log',
}

MPMATH_TRANSLATIONS = {
    'Abs': 'fabs',
    'elliptic_k': 'ellipk',
    'elliptic_f': 'ellipf',
    'elliptic_e': 'ellipe',
    'elliptic_pi': 'ellippi',
    'ceiling': 'ceil',
    'chebyshevt': 'chebyt',
    'chebyshevu': 'chebyu',
    'E': 'e',
    'I': 'j',
    'ln': 'log',
    # 'lowergamma': 'lower_gamma',
    'oo': 'inf',
    # 'uppergamma': 'upper_gamma',
    'LambertW': 'lambertw',
    'MutableDenseMatrix': 'matrix',
    'ImmutableMatrix': 'matrix',
    'conjugate': 'conj',
    'dirichlet_eta': 'altzeta',
    'Ei': 'ei',
    'Shi': 'shi',
    'Chi': 'chi',
    'Si': 'si',
    'Ci': 'ci',
    'Ynm': 'spherharm',
    'RisingFactorial': 'rf',
    'FallingFactorial': 'ff',
}

NUMPY_TRANSLATIONS = {
    'Abs': 'abs',
    'acos': 'arccos',
    'acosh': 'arccosh',
    'arg': 'angle',
    'asin': 'arcsin',
    'asinh': 'arcsinh',
    'atan': 'arctan',
    'atan2': 'arctan2',
    'atanh': 'arctanh',
    'ceiling': 'ceil',
    'E': 'e',
    'im': 'imag',
    'ln': 'log',
    'Mod': 'mod',
    'oo': 'inf',
    're': 'real',
    'SparseMatrix': 'array',
    'ImmutableSparseMatrix': 'array',
    'Matrix': 'array',
    'MutableDenseMatrix': 'array',
    'ImmutableMatrix': 'array',
    'ImmutableDenseMatrix': 'array',
}

# Available modules:
MODULES = {
    'math': (MATH, MATH_DEFAULT, MATH_TRANSLATIONS, ('from math import *',)),
    'mpmath': (MPMATH, MPMATH_DEFAULT, MPMATH_TRANSLATIONS, ('from mpmath import *',)),
    'numpy': (NUMPY, NUMPY_DEFAULT, NUMPY_TRANSLATIONS, ("import_module('numpy')",)),
    'diofant': (DIOFANT, DIOFANT_DEFAULT, {}, (
        'from diofant.functions import *',
        'from diofant.matrices import *',
        'from diofant import Sum, Integral, pi, oo, nan, zoo, E, I')),
}


def _import(module):
    """
    Creates a global translation dictionary for module.

    The argument module has to be one of the following strings: "math",
    "mpmath", "numpy", "diofant".
    These dictionaries map names of python functions to their equivalent in
    other modules.

    """
    try:
        namespace, namespace_default, translations, import_commands = MODULES[
            module]
    except KeyError:
        raise NameError(
            f"'{module}' module can't be used for lambdification")

    # Clear namespace or exit
    if namespace != namespace_default:
        namespace.clear()
        namespace.update(namespace_default)

    for import_command in import_commands:
        if import_command.startswith('import_module'):
            module = eval(import_command)

            if module is not None:
                namespace.update(module.__dict__)
                continue
        else:
            exec(import_command, {}, namespace)
            continue

        raise ImportError(
            f"can't import '{module}' with '{import_command}' command")

    # Add translated names to namespace
    for diofantname, translation in translations.items():
        namespace[diofantname] = namespace[translation]


@doctest_depends_on(modules=('numpy'))
def lambdify(args, expr, modules=None, printer=None, use_imps=True,
             dummify=True):
    """
    Returns a lambda function for fast calculation of numerical values.

    If not specified differently by the user, ``modules`` defaults to
    ``["numpy"]`` if NumPy is installed, and ``["math", "mpmath", "sympy"]``
    if it isn't, that is, Diofant functions are replaced as far as possible by
    either ``numpy`` functions if available, and Python's standard library
    ``math``, or ``mpmath`` functions otherwise. To change this behavior, the
    "modules" argument can be used. It accepts:

     - the strings "math", "mpmath", "numpy", "diofant"
     - any modules (e.g. math)
     - dictionaries that map names of diofant functions to arbitrary functions
     - lists that contain a mix of the arguments above, with higher priority
       given to entries appearing first.

    The default behavior is to substitute all arguments in the provided
    expression with dummy symbols. This allows for applied functions (e.g.
    f(t)) to be supplied as arguments. Call the function with dummify=False if
    dummy substitution is unwanted (and `args` is not a string). If you want
    to view the lambdified function or provide "diofant" as the module, you
    should probably set dummify=False.

    In previous releases ``lambdify`` replaced ``Matrix`` with ``numpy.matrix``
    by default. As of release 0.7.7 ``numpy.array`` is the default.
    To get the old default behavior you must pass in ``[{'ImmutableMatrix':
    numpy.matrix}, 'numpy']`` to the ``modules`` kwarg.

    (1) Use one of the provided modules:

        >>> f = lambdify(x, sin(x), 'math')

        Attention: Functions that are not in the math module will throw a name
                   error when the lambda function is evaluated! So this would
                   be better:

        >>> f = lambdify(x, sin(x)*gamma(x), ('math', 'mpmath', 'diofant'))

    (2) Use some other module:

        >>> import numpy
        >>> f = lambdify((x, y), tan(x*y), numpy)

        Attention: There are naming differences between numpy and diofant. So if
                   you simply take the numpy module, e.g. diofant.atan will not be
                   translated to numpy.arctan. Use the modified module instead
                   by passing the string "numpy":

        >>> f = lambdify((x, y), tan(x*y), 'numpy')
        >>> f(1, 2)
        -2.18503986326
        >>> from numpy import array
        >>> f(array([1, 2, 3]), array([2, 3, 5]))
        [-2.18503986 -0.29100619 -0.8559934 ]

    (3) Use a dictionary defining custom functions:

        >>> def my_cool_function(x):
        ...     return f'sin({x}) is cool'
        >>> myfuncs = {'sin': my_cool_function}
        >>> f = lambdify(x, sin(x), myfuncs)
        >>> f(1)
        'sin(1) is cool'

    Examples
    ========

    >>> from diofant.abc import w

    >>> f = lambdify(x, x**2)
    >>> f(2)
    4
    >>> f = lambdify((x, y, z), [z, y, x])
    >>> f(1, 2, 3)
    [3, 2, 1]
    >>> f = lambdify(x, sqrt(x))
    >>> f(4)
    2.0
    >>> f = lambdify((x, y), sin(x*y)**2)
    >>> f(0, 5)
    0.0
    >>> row = lambdify((x, y), Matrix((x, x + y)).T, modules='diofant')
    >>> row(1, 2)
    Matrix([[1, 3]])

    Tuple arguments are handled and the lambdified function should
    be called with the same type of arguments as were used to create
    the function.:

    >>> f = lambdify((x, (y, z)), x + y)
    >>> f(1, (2, 4))
    3

    A more robust way of handling this is to always work with flattened
    arguments:

    >>> args = w, (x, (y, z))
    >>> vals = 1, (2, (3, 4))
    >>> f = lambdify(flatten(args), w + x + y + z)
    >>> f(*flatten(vals))
    10

    Functions present in `expr` can also carry their own numerical
    implementations, in a callable attached to the ``_imp_``
    attribute.  Usually you attach this using the
    ``implemented_function`` factory:

    >>> f = implemented_function(Function('f'), lambda x: x+1)
    >>> func = lambdify(x, f(x))
    >>> func(4)
    5

    ``lambdify`` always prefers ``_imp_`` implementations to implementations
    in other namespaces, unless the ``use_imps`` input parameter is False.

    """
    from ..core import Symbol
    from ..printing.lambdarepr import MpmathPrinter, NumPyPrinter
    from .iterables import flatten

    module_provided = True

    # If the user hasn't specified any modules, use what is available.
    if modules is None:
        module_provided = False

        try:
            _import('numpy')
        except ImportError:
            # Use either numpy (if available) or python.math where possible.
            # XXX: This leads to different behaviour on different systems and
            #      might be the reason for irreproducible errors.
            modules = ['math', 'mpmath', 'diofant']
        else:
            modules = ['numpy']

    # Get the needed namespaces.
    namespaces = []
    # First find any function implementations
    if use_imps:
        namespaces.append(_imp_namespace(expr))
    # Check for dict before iterating
    if isinstance(modules, (dict, str)) or not hasattr(modules, '__iter__'):
        namespaces.append(modules)
    else:
        namespaces += list(modules)
    # fill namespace with first having highest priority
    namespace = {}
    for m in namespaces[::-1]:
        buf = _get_namespace(m)
        namespace.update(buf)

    if hasattr(expr, 'atoms'):
        # Try if you can extract symbols from the expression.
        # Move on if expr.atoms in not implemented.
        syms = expr.atoms(Symbol)
        for term in syms:
            namespace.update({str(term): term})

    if 'numpy' in namespaces and printer is None:
        printer = NumPyPrinter

    if 'mpmath' in namespaces and printer is None:
        printer = MpmathPrinter

    # Get the names of the args, for creating a docstring
    if not iterable(args):
        args = args,
    names = []
    for n, var in enumerate(args):
        if hasattr(var, 'name'):
            names.append(var.name)
        else:
            # Cannot infer name with certainty. arg_# will have to do.
            names.append('arg_' + str(n))

    # Create lambda function.
    lstr = lambdastr(args, expr, printer=printer, dummify=dummify)
    flat = '__flatten_args__'

    if flat in lstr:
        namespace.update({flat: flatten})
    func = eval(lstr, namespace)
    # For numpy lambdify, wrap all input arguments in arrays.
    if module_provided and 'numpy' in namespaces:
        def array_wrap(funcarg):
            def wrapper(*argsx, **kwargsx):
                return funcarg(*[namespace['asarray'](i) for i in argsx], **kwargsx)
            return wrapper
        func = array_wrap(func)
    # Apply the docstring
    sig = f"func({', '.join(str(i) for i in names)})"
    sig = textwrap.fill(sig, subsequent_indent=' '*8)
    expr_str = str(expr)
    if len(expr_str) > 78:
        expr_str = textwrap.wrap(expr_str, 75)[0] + '...'
    func.__doc__ = f"""Created with lambdify. Signature:

{sig}

Expression:

{expr_str}"""
    return func


def _get_namespace(m):
    """This is used by _lambdify to parse its arguments."""
    if isinstance(m, str):
        _import(m)
        return MODULES[m][0]
    elif isinstance(m, dict):
        return m
    elif hasattr(m, '__dict__'):
        return m.__dict__
    else:
        raise TypeError(f'Argument must be either a string, dict or module but it is: {m}')


def lambdastr(args, expr, printer=None, dummify=False):
    """
    Returns a string that can be evaluated to a lambda function.

    Examples
    ========

    >>> lambdastr(x, x**2)
    'lambda x: (x**2)'
    >>> lambdastr((x, y, z), [z, y, x])
    'lambda x,y,z: ([z, y, x])'

    Although tuples may not appear as arguments to lambda in Python 3,
    lambdastr will create a lambda function that will unpack the original
    arguments so that nested arguments can be handled:

    >>> lambdastr((x, (y, z)), x + y)
    'lambda _0,_1: (lambda x,y,z: (x + y))(*list(__flatten_args__([_0,_1])))'

    """
    # Transforming everything to strings.
    from ..core import Dummy, Function, Symbol, sympify
    from ..utilities import flatten

    if printer is not None:
        if inspect.isfunction(printer):
            lambdarepr = printer
        else:
            if inspect.isclass(printer):
                def lambdarepr(expr):
                    return printer().doprint(expr)
            else:
                def lambdarepr(expr):
                    return printer.doprint(expr)
    else:
        # XXX: This has to be done here because of circular imports
        from ..printing.lambdarepr import lambdarepr

    def sub_args(args, dummies_dict):
        if isinstance(args, str):
            return args
        elif iterable(args):
            dummies = flatten([sub_args(a, dummies_dict) for a in args])
            return ','.join(str(a) for a in dummies)
        else:
            # Sub in dummy variables for functions or symbols
            if isinstance(args, (Function, Symbol)):
                dummies = Dummy()
                dummies_dict.update({args: dummies})
                return lambdarepr(dummies)
            else:
                return lambdarepr(args)

    def sub_expr(expr, dummies_dict):
        try:
            expr = sympify(expr).xreplace(dummies_dict)
        except (TypeError, AttributeError):
            if isinstance(expr, dict):
                k = [sub_expr(sympify(a), dummies_dict) for a in expr]
                v = [sub_expr(sympify(a), dummies_dict) for a in expr.values()]
                expr = dict(zip(k, v))
            elif isinstance(expr, tuple):
                expr = tuple(sub_expr(sympify(a), dummies_dict) for a in expr)
            elif isinstance(expr, list):
                expr = [sub_expr(sympify(a), dummies_dict) for a in expr]
        return expr

    # Transform args
    def isiter(l):
        return iterable(l, exclude=(str,))

    if isiter(args) and any(isiter(i) for i in args):
        import re
        dum_args = [lambdarepr(Dummy(str(i))) for i in range(len(args))]
        iter_args = ','.join([i if isiter(a) else i
                              for i, a in zip(dum_args, args)])
        lstr = lambdastr(flatten(args), expr, printer=printer, dummify=dummify)
        flat = '__flatten_args__'
        rv = 'lambda %s: (%s)(*list(%s([%s])))' % (
            ','.join(dum_args), lstr, flat, iter_args)
        if len(re.findall(r'\b%s\b' % flat, rv)) > 1:
            raise ValueError(f'the name {flat} is reserved by lambdastr')
        return rv

    dummies_dict = {}
    if dummify:
        args = sub_args(args, dummies_dict)
    else:
        if iterable(args):
            args = ','.join(str(a) for a in args)

    # Transform expr
    if dummify:
        if isinstance(expr, str):
            pass
        else:
            expr = sub_expr(expr, dummies_dict)
    expr = lambdarepr(expr)

    return f'lambda {args}: ({expr})'


def _imp_namespace(expr, namespace=None):
    """Return namespace dict with function implementations

    We need to search for functions in anything that can be thrown at
    us - that is - anything that could be passed as `expr`.  Examples
    include diofant expressions, as well as tuples, lists and dicts that may
    contain diofant expressions.

    Parameters
    ==========

    expr : object
       Something passed to lambdify, that will generate valid code from
       ``str(expr)``.
    namespace : None or mapping
       Namespace to fill.  None results in new empty dict

    Returns
    =======

    namespace : dict
       dict with keys of implemented function names within `expr` and
       corresponding values being the numerical implementation of
       function

    Examples
    ========

    >>> f = implemented_function(Function('f'), lambda x: x+1)
    >>> g = implemented_function(Function('g'), lambda x: x*10)
    >>> namespace = _imp_namespace(f(g(x)))
    >>> sorted(namespace)
    ['f', 'g']

    """
    # Delayed import to avoid circular imports
    from ..core.function import FunctionClass
    if namespace is None:
        namespace = {}
    # tuples, lists, dicts are valid expressions
    if is_sequence(expr):
        for arg in expr:
            _imp_namespace(arg, namespace)
        return namespace
    elif isinstance(expr, dict):
        for key, val in expr.items():
            # functions can be in dictionary keys
            _imp_namespace(key, namespace)
            _imp_namespace(val, namespace)
        return namespace
    # diofant expressions may be Functions themselves
    func = getattr(expr, 'func', None)
    if isinstance(func, FunctionClass):
        imp = getattr(func, '_imp_', None)
        if imp is not None:
            name = expr.func.__name__
            if name in namespace and namespace[name] != imp:
                raise ValueError('We found more than one '
                                 'implementation with name '
                                 '"%s"' % name)
            namespace[name] = imp
    # and / or they may take Functions as arguments
    if hasattr(expr, 'args'):
        for arg in expr.args:
            _imp_namespace(arg, namespace)
    return namespace


def implemented_function(symfunc, implementation):
    """Add numerical ``implementation`` to function ``symfunc``.

    ``symfunc`` can be an ``UndefinedFunction`` instance, or a name string.
    In the latter case we create an ``UndefinedFunction`` instance with that
    name.

    Be aware that this is a quick workaround, not a general method to create
    special symbolic functions. If you want to create a symbolic function to be
    used by all the machinery of Diofant you should subclass the ``Function``
    class.

    Parameters
    ==========

    symfunc : ``str`` or ``UndefinedFunction`` instance
       If ``str``, then create new ``UndefinedFunction`` with this as
       name.  If `symfunc` is a diofant function, attach implementation to it.
    implementation : callable
       numerical implementation to be called by ``evalf()`` or ``lambdify``

    Returns
    =======

    afunc : diofant.FunctionClass instance
       function with attached implementation

    Examples
    ========

    >>> f = implemented_function(Function('f'), lambda x: x+1)
    >>> lam_f = lambdify(x, f(x))
    >>> lam_f(4)
    5

    """
    # Delayed import to avoid circular imports
    from ..core.function import UndefinedFunction

    # if name, create function to hold implementation
    if isinstance(symfunc, str):
        symfunc = UndefinedFunction(symfunc)
    elif not isinstance(symfunc, UndefinedFunction):
        raise ValueError('symfunc should be either a string or'
                         ' an UndefinedFunction instance.')
    # We need to attach as a method because symfunc will be a class
    symfunc._imp_ = staticmethod(implementation)
    return symfunc
