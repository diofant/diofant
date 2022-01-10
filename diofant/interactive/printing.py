import builtins
import fractions
import os
import sys

from ..printing import latex, pretty, sstrrepr
from ..printing.printer import Printer


def _init_python_printing(stringify_func):
    """Setup printing in Python interactive session."""

    def _displayhook(arg):
        """Python's pretty-printer display hook.

        This function was adapted from PEP 217.

        """
        if arg is not None:
            builtins._ = None
            if isinstance(arg, str):
                print(repr(arg))
            else:
                print(stringify_func(arg))
            builtins._ = arg

    sys.displayhook = _displayhook


def _init_ipython_printing(ip, stringify_func):
    """Setup printing in IPython interactive session."""

    def _print_plain_text(arg, p, cycle):
        p.text(stringify_func(arg))

    def _print_latex_text(o):
        return latex(o, mode='equation*')

    printable_types = [float, tuple, list, set, frozenset, dict,
                       int, fractions.Fraction]

    plaintext_formatter = ip.display_formatter.formatters['text/plain']
    latex_formatter = ip.display_formatter.formatters['text/latex']

    for cls in printable_types:
        plaintext_formatter.for_type(cls, _print_plain_text)

    for cls in printable_types:
        latex_formatter.for_type(cls, _print_latex_text)


def init_printing(no_global=False, pretty_print=None, **settings):
    r"""Initializes pretty-printer depending on the environment.

    Parameters
    ==========

    no_global : boolean
        If True, the settings become system wide;
        if False, use just for this console/session.
    pretty_print : boolean or None
        Enable pretty printer (turned on by default for IPython, but
        disabled for plain Python console).
    \*\*settings : dict
        A dictionary of default settings for printers.

    Notes
    =====

    This function runs automatically for wildcard imports (e.g.
    for ``from diofant import *``) in interactive sessions.

    Examples
    ========

    >>> from diofant.abc import theta
    >>> sqrt(5)
    sqrt(5)
    >>> init_printing(pretty_print=True, no_global=True)
    >>> sqrt(5)
      ___
    ╲╱ 5
    >>> theta
    θ
    >>> init_printing(pretty_print=True, use_unicode=False, no_global=True)
    >>> theta
    theta
    >>> init_printing(pretty_print=True, order='grevlex', no_global=True)
    >>> y + x + y**2 + x**2
     2    2
    x  + y  + x + y

    """
    try:
        ip = get_ipython()
    except NameError:
        ip = None

    # guess unicode support
    unicode_term = False
    TERM = os.environ.get('TERM', '')
    if TERM != '' and not TERM.endswith('linux'):
        unicode_term = True
    if settings.get('use_unicode') is None:
        settings['use_unicode'] = bool(unicode_term)

    if ip:
        stringify_func = pretty if pretty_print is not False else sstrrepr
    else:
        stringify_func = sstrrepr if not pretty_print else pretty

    if no_global:
        _stringify_func = stringify_func

        def stringify_func(expr):  # pylint: disable=function-redefined
            return _stringify_func(expr, **settings)
    else:
        Printer.set_global_settings(**settings)

    if ip:
        _init_ipython_printing(ip, stringify_func)
    else:
        _init_python_printing(stringify_func)
