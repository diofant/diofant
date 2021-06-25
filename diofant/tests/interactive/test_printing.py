import sys

from diofant import init_printing, sqrt
from diofant.abc import theta, x, y


def test_init_printing(capsys):
    init_printing()
    sys.displayhook(sqrt(5))
    assert capsys.readouterr().out == 'sqrt(5)\n'
    sys.displayhook('xyz')
    assert capsys.readouterr().out == "'xyz'\n"
    sys.displayhook(None)
    assert capsys.readouterr().out == ''

    init_printing(pretty_print=True, no_global=True)
    sys.displayhook(sqrt(5))
    assert capsys.readouterr().out == \
        """\
  ___\n\
╲╱ 5 \n"""
    sys.displayhook(theta)
    assert capsys.readouterr().out == 'θ\n'

    init_printing(pretty_print=True, use_unicode=False, no_global=True)
    sys.displayhook(theta)
    assert capsys.readouterr().out == 'theta\n'

    init_printing(pretty_print=True, order='grevlex', no_global=True)
    sys.displayhook(y + x + y**2 + x**2)
    assert capsys.readouterr().out == \
        """\
 2    2        \n\
x  + y  + x + y\n"""
