class DefaultPrinting:
    """
    The default implementation of printing for Diofant classes.

    This implements a hack that allows us to print elements of built-in
    Python containers in a readable way. Natively Python uses ``repr()``
    even if ``str()`` was explicitly requested. Mix in this trait into
    a class to get proper default printing.

    """

    # Note, we always use the default ordering (lex) in __str__ and __repr__,
    # regardless of the global setting. See issue sympy/sympy#5487.
    def __str__(self):
        from .str import sstr
        return sstr(self, order=None)

    def __repr__(self):
        from .repr import srepr
        return srepr(self, order=None)

    def _repr_pretty_(self, p, cycle):
        from .pretty import pretty
        p.text(pretty(self))

    def _repr_latex_(self):
        from .latex import latex
        return latex(self, mode='equation')
