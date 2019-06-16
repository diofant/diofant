"""Printing subsystem driver

Diofant's printing system works the following way: Any expression can be
passed to a designated Printer who then is responsible to return an
adequate representation of that expression.

The basic concept is the following:
  1. Let the object print itself if it knows how.
  2. Take the best fitting method defined in the printer.
  3. As fall-back use the emptyPrinter method for the printer.

Some more information how the single concepts work and who should use which:

1. The object prints itself

    This was the original way of doing printing in diofant. Every class had
    its own latex, mathml, str and repr methods, but it turned out that it
    is hard to produce a high quality printer, if all the methods are spread
    out that far. Therefore all printing code was combined into the different
    printers, which works great for built-in diofant objects, but not that
    good for user defined classes where it is inconvenient to patch the
    printers.

    Nevertheless, to get a fitting representation, the printers look for a
    specific method in every object, that will be called if it's available
    and is then responsible for the representation. The name of that method
    depends on the specific printer and is defined under
    Printer.printmethod.

2. Take the best fitting method defined in the printer.

    The printer loops through expr classes (class + its bases), and tries
    to dispatch the work to _print_<EXPR_CLASS>

    e.g., suppose we have the following class hierarchy::

            Basic
            |
            Atom
            |
            Number
            |
        Rational

    then, for expr=Rational(...), in order to dispatch, we will try
    calling printer methods as shown in the figure below::

        p._print(expr)
        |
        |-- p._print_Rational(expr)
        |
        |-- p._print_Number(expr)
        |
        |-- p._print_Atom(expr)
        |
        `-- p._print_Basic(expr)

    if ._print_Rational method exists in the printer, then it is called,
    and the result is returned back.

    otherwise, we proceed with trying Rational bases in the inheritance
    order.

3. As fall-back use the emptyPrinter method for the printer.

    As fall-back self.emptyPrinter will be called with the expression. If
    not defined in the Printer subclass this will be the same as str(expr).
"""


class Printer:
    """Generic printer

    Its job is to provide infrastructure for implementing new printers easily.

    Basically, if you want to implement a printer, all you have to do is:

    1. Subclass Printer.

    2. Define Printer.printmethod in your subclass.
       If a object has a method with that name, this method will be used
       for printing.

    3. In your subclass, define ``_print_<CLASS>`` methods

       For each class you want to provide printing to, define an appropriate
       method how to do it. For example if you want a class FOO to be printed in
       its own way, define _print_FOO::

           def _print_FOO(self, e):
               ...

       this should return how FOO instance e is printed

       Also, if ``BAR`` is a subclass of ``FOO``, ``_print_FOO(bar)`` will
       be called for instance of ``BAR``, if no ``_print_BAR`` is provided.
       Thus, usually, we don't need to provide printing routines for every
       class we want to support -- only generic routine has to be provided
       for a set of classes.

       A good example for this are functions - for example ``PrettyPrinter``
       only defines ``_print_Function``, and there is no ``_print_sin``,
       ``_print_tan``, etc...

       On the other hand, a good printer will probably have to define
       separate routines for ``Symbol``, ``Atom``, ``Number``, ``Integral``,
       ``Limit``, etc...

    4. If convenient, override ``self.emptyPrinter``

       This callable will be called to obtain printing result as a last resort,
       that is when no appropriate print method was found for an expression.

    Examples
    ========

    Here we will overload ``StrPrinter``.

    >>> from diofant.printing.str import StrPrinter

    >>> class CustomStrPrinter(StrPrinter):
    ...     def _print_Derivative(self, expr):
    ...         return str(expr.args[0].func) + "'"*len(expr.args[1:])
    >>> def mystr(e):
    ...     return CustomStrPrinter().doprint(e)
    >>> t = Symbol('t')
    >>> x = Function('x')(t)
    >>> print(mystr(x.diff(t, 2)))
    x''

    """

    _global_settings = {}

    _default_settings = {}

    emptyPrinter = str
    printmethod = None

    def __init__(self, settings=None):
        import distutils
        import distutils.version
        from ..external import import_module

        self._str = str

        self._settings = self._default_settings.copy()

        for key, val in self._global_settings.items():
            if key in self._default_settings:
                self._settings[key] = val

        if settings is not None:
            self._settings.update(settings)

            if len(self._settings) > len(self._default_settings):
                for key in self._settings:  # pragma: no branch
                    if key not in self._default_settings:
                        raise TypeError("Unknown setting '%s'." % key)

        # _print_level is the number of times self._print() was recursively
        # called. See StrPrinter._print_Float() for an example of usage
        self._print_level = 0

        numpy = import_module("numpy")
        if numpy is not None:  # pragma: no cover
            kwargs = {'formatter': {'object': str}}
            if numpy.__version__ >= distutils.version.LooseVersion('1.14.0'):
                kwargs['legacy'] = "1.13"
            numpy.set_printoptions(**kwargs)

    @classmethod
    def set_global_settings(cls, **settings):
        """Set system-wide printing settings."""
        for key, val in settings.items():
            if val is not None:
                cls._global_settings[key] = val

    @property
    def order(self):
        if 'order' in self._settings:
            return self._settings['order']
        else:
            raise AttributeError("No order defined.")

    def doprint(self, expr):
        """Returns printer's representation for expr (as a string)."""
        return self._str(self._print(expr))

    def _print(self, expr, *args, **kwargs):
        """Internal dispatcher

        Tries the following concepts to print an expression:
            1. Let the object print itself if it knows how.
            2. Take the best fitting method defined in the printer.
            3. As fall-back use the emptyPrinter method for the printer.

        """
        self._print_level += 1
        try:
            # If the printer defines a name for a printing method
            # (Printer.printmethod) and the object knows for itself how it
            # should be printed, use that method.
            if (self.printmethod and hasattr(expr, self.printmethod)
                    and not isinstance(expr, type)):
                return getattr(expr, self.printmethod)(self, *args, **kwargs)

            # See if the class of expr is known, or if one of its super
            # classes is known, and use that print function
            for cls in type(expr).__mro__:
                printmethod = '_print_' + cls.__name__
                if hasattr(self, printmethod):
                    if cls.__class__.__name__ == 'UndefinedFunction':
                        # XXX a hack for sympy/sympy#6853
                        continue
                    return getattr(self, printmethod)(expr, *args, **kwargs)

            # Unknown object, fall back to the emptyPrinter.
            return self.emptyPrinter(expr)
        finally:
            self._print_level -= 1
