"""Tools for setting up printing in interactive sessions. """

import builtins
import os
import sys

from diofant import latex as default_latex
from diofant.utilities.misc import debug


def _init_python_printing(stringify_func):
    """Setup printing in Python interactive session. """

    def _displayhook(arg):
        """Python's pretty-printer display hook.

        This function was adapted from PEP 217.
        """
        if arg is not None:
            builtins._ = None
            print(stringify_func(arg))
            builtins._ = arg

    sys.displayhook = _displayhook


def _init_ipython_printing(ip, stringify_func, use_latex,
                           print_builtin,
                           latex_printer):
    """Setup printing in IPython interactive session. """

    import IPython
    from diofant.core.basic import Basic
    from diofant.matrices.matrices import MatrixBase

    latex = latex_printer or default_latex

    def _print_plain(arg, p, cycle):
        """caller for pretty, for use in IPython"""
        if _can_print_latex(arg):
            p.text(stringify_func(arg))
        else:
            p.text(IPython.lib.pretty.pretty(arg))

    def _can_print_latex(o):
        """Return True if type o can be printed with LaTeX.

        If o is a container type, this is True if and only if every element of
        o can be printed with LaTeX.
        """
        if isinstance(o, (list, tuple, set, frozenset)):
            return all(_can_print_latex(i) for i in o)
        elif isinstance(o, dict):
            return all(_can_print_latex(i) and _can_print_latex(o[i]) for i in o)
        elif isinstance(o, bool):
            return False
        # TODO : Investigate if "elif hasattr(o, '_latex')" is more useful
        # to use here, than these explicit imports.
        elif isinstance(o, (Basic, MatrixBase)):
            return True
        elif isinstance(o, (float, int)) and print_builtin:
            return True
        return False

    def _print_latex_text(o):
        """
        A function to generate the latex representation of diofant expressions.
        """
        if _can_print_latex(o):
            s = latex(o, mode='plain')
            s = s.replace(r'\dag', r'\dagger')
            s = s.strip('$')
            return '$$%s$$' % s

    printable_types = [Basic, MatrixBase, float, tuple, list, set,
                       frozenset, dict, int]

    plaintext_formatter = ip.display_formatter.formatters['text/plain']

    for cls in printable_types:
        plaintext_formatter.for_type(cls, _print_plain)

    latex_formatter = ip.display_formatter.formatters['text/latex']
    if use_latex:
        debug("init_printing: using mathjax formatter")
        for cls in printable_types:
            latex_formatter.for_type(cls, _print_latex_text)
    else:
        debug("init_printing: not using text/latex formatter")
        for cls in printable_types:
            if cls in latex_formatter.type_printers:
                latex_formatter.type_printers.pop(cls)


def init_printing(pretty_print=True, order=None, use_unicode=None,
                  use_latex=None, wrap_line=None, num_columns=None,
                  no_global=False, ip=None,
                  print_builtin=True,
                  str_printer=None, pretty_printer=None,
                  latex_printer=None):
    """Initializes pretty-printer depending on the environment.

    Parameters
    ==========

    pretty_print: boolean
        If True, use pretty_print to stringify or the provided pretty
        printer; if False, use sstrrepr to stringify or the provided string
        printer.
    order: string or None
        There are a few different settings for this parameter:
        lex (default), which is lexographic order;
        grlex, which is graded lexographic order;
        grevlex, which is reversed graded lexographic order;
        None, which sets it to lex.
    use_unicode: boolean or None
        If True, use unicode characters;
        if False, do not use unicode characters.
    use_latex: string, boolean, or None
        If True, use default latex rendering in GUI interfaces;
        if False, do not use latex rendering;
    wrap_line: boolean
        If True, lines will wrap at the end; if False, they will not wrap
        but continue as one line. This is only relevant if `pretty_print` is
        True.
    num_columns: int or None
        If int, number of columns before wrapping is set to num_columns; if
        None, number of columns before wrapping is set to terminal width.
        This is only relevant if `pretty_print` is True.
    no_global: boolean
        If True, the settings become system wide;
        if False, use just for this console/session.
    ip: An interactive console
        This can either be an instance of IPython,
        or a class that derives from code.InteractiveConsole.
    print_builtin: boolean, optional, default=True
        If true then floats and integers will be printed. If false the
        printer will only print Diofant types.
    str_printer: function, optional, default=None
        A custom string printer function. This should mimic
        diofant.printing.sstrrepr().
    pretty_printer: function, optional, default=None
        A custom pretty printer. This should mimic diofant.printing.pretty().
    latex_printer: function, optional, default=None
        A custom LaTeX printer. This should mimic diofant.printing.latex()
        This should mimic diofant.printing.latex().

    Examples
    ========

    >>> from diofant.interactive import init_printing
    >>> from diofant import Symbol, sqrt
    >>> from diofant.abc import x, y
    >>> sqrt(5)
    sqrt(5)
    >>> init_printing(pretty_print=True)
    >>> sqrt(5)
      ___
    \/ 5
    >>> theta = Symbol('theta')
    >>> init_printing(use_unicode=True)
    >>> theta
    \u03b8
    >>> init_printing(use_unicode=False)
    >>> theta
    theta
    >>> init_printing(order='grevlex')
    >>> y + x + y**2 + x**2
     2    2
    x  + y  + x + y
    >>> init_printing(pretty_print=False, use_unicode=False, order='lex')
    >>> y + x + y**2 + x**2
    x**2 + x + y**2 + y
    """
    from diofant.printing.printer import Printer

    if pretty_print:
        if pretty_printer is not None:
            stringify_func = pretty_printer
        else:
            from diofant.printing import pretty as stringify_func
    else:
        if str_printer is not None:
            stringify_func = str_printer
        else:
            from diofant.printing import sstrrepr as stringify_func

    # Even if ip is not passed, double check that not in IPython shell
    in_ipython = False
    if ip is None:
        try:
            ip = get_ipython()
        except NameError:
            pass
        else:
            in_ipython = (ip is not None)
    else:
        in_ipython = True

    if in_ipython:
        if use_unicode is None:
            use_unicode = False if os.environ.get('TERM', '').endswith('linux') else True
        if use_latex is None:
            use_latex = True

    if not no_global:
        Printer.set_global_settings(order=order, use_unicode=use_unicode,
                                    wrap_line=wrap_line,
                                    num_columns=num_columns)
    else:
        _stringify_func = stringify_func

        if pretty_print:
            def stringify_func(expr):
                return _stringify_func(expr, order=order,
                                       use_unicode=use_unicode,
                                       wrap_line=wrap_line,
                                       num_columns=num_columns)
        else:
            def stringify_func(expr):
                return _stringify_func(expr, order=order)

    if in_ipython:
        _init_ipython_printing(ip, stringify_func, use_latex,
                               print_builtin, latex_printer)
    else:
        _init_python_printing(stringify_func)
