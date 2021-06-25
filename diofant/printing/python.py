import keyword as kw

from .repr import ReprPrinter
from .str import StrPrinter


# A list of classes that should be printed using StrPrinter
STRPRINT = ('Add', 'Infinity', 'Integer', 'Mul', 'NegativeInfinity',
            'Pow', 'Zero')


class PythonPrinter(ReprPrinter, StrPrinter):
    """A printer which converts an expression into its Python interpretation."""

    def __init__(self, settings=None):
        ReprPrinter.__init__(self)
        StrPrinter.__init__(self, settings)
        self.symbols = []
        self.functions = []

        # Create print methods for classes that should use StrPrinter instead
        # of ReprPrinter.
        for name in STRPRINT:
            f_name = f'_print_{name}'
            f = getattr(StrPrinter, f_name)
            setattr(PythonPrinter, f_name, f)

    def _print_Function(self, expr):
        import diofant
        func = expr.func.__name__
        if not hasattr(diofant, func) and func not in self.functions:
            self.functions.append(func)
        return StrPrinter._print_Function(self, expr)

    # procedure (!) for defining symbols which have be defined in python()
    def _print_Symbol(self, expr):
        symbol = self._str(expr)
        if symbol not in self.symbols:
            self.symbols.append(symbol)
        return StrPrinter._print_Symbol(self, expr)
    _print_BaseSymbol = StrPrinter._print_BaseSymbol


def python(expr, **settings):
    """Return Python interpretation of passed expression
    (can be passed to the exec() function without any modifications)
    """
    from ..core import Function, Symbol

    printer = PythonPrinter(settings)
    exprp = printer.doprint(expr)

    result = ''
    # Returning found symbols and functions
    renamings = {}
    for symbolname in printer.symbols:
        newsymbolname = symbolname
        # Escape symbol names that are reserved python keywords
        if kw.iskeyword(newsymbolname):
            while True:
                newsymbolname += '_'
                if (newsymbolname not in printer.symbols and
                        newsymbolname not in printer.functions):
                    renamings[Symbol(symbolname)] = Symbol(newsymbolname)
                    break
        result += newsymbolname + " = Symbol('" + symbolname + "')\n"

    for functionname in printer.functions:
        newfunctionname = functionname
        # Escape function names that are reserved python keywords
        if kw.iskeyword(newfunctionname):
            while True:
                newfunctionname += '_'
                if (newfunctionname not in printer.symbols and
                        newfunctionname not in printer.functions):
                    renamings[Function(functionname)] = Function(newfunctionname)
                    break
        result += newfunctionname + " = Function('" + functionname + "')\n"

    if not len(renamings) == 0:
        exprp = expr.subs(renamings)
    result += 'e = ' + printer._str(exprp)
    return result
