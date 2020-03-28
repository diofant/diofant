"""Implicit plotting tests."""

import tempfile
import warnings

import pytest

from diofant import (And, Eq, I, Or, cos, exp, pi, plot_implicit, re, sin,
                     symbols, tan)
from diofant.abc import x, y


__all__ = ()

matplotlib = pytest.importorskip('matplotlib', minversion='1.1.0')


def tmp_file(name=''):
    return tempfile.NamedTemporaryFile(suffix='.png').name


def plot_and_save(name):
    # implicit plot tests
    plot_implicit(Eq(y, cos(x)), (x, -5, 5), (y, -2, 2)).save(tmp_file(name))
    plot_implicit(Eq(y**2, x**3 - x), (x, -5, 5),
                  (y, -4, 4), show=False).save(tmp_file(name))
    plot_implicit(y > 1 / x, (x, -5, 5),
                  (y, -2, 2), show=False).save(tmp_file(name))
    plot_implicit(y < 1 / tan(x), (x, -5, 5),
                  (y, -2, 2), show=False).save(tmp_file(name))
    plot_implicit(y >= 2 * sin(x) * cos(x), (x, -5, 5),
                  (y, -2, 2), show=False).save(tmp_file(name))
    plot_implicit(y <= x**2, (x, -3, 3),
                  (y, -1, 5), show=False).save(tmp_file(name))
    plot_implicit(Or(And(x < y, x > y**2), sin(y) > x),
                  show=False).save(tmp_file(name))
    plot_implicit(Or(And(x < y, x > y**2), sin(y)),
                  show=False).save(tmp_file(name))
    plot_implicit(Or(x < y, sin(x) >= y), show=False).save(tmp_file(name))

    # Test all input args for plot_implicit
    plot_implicit(Eq(y**2, x**3 - x), show=False).save(tmp_file(name))
    plot_implicit(Eq(y**2, x**3 - x), adaptive=False, show=False).save(tmp_file(name))
    plot_implicit(Eq(y**2, x**3 - x), adaptive=False, points=500, show=False).save(tmp_file(name))
    plot_implicit(y > x, (x, -5, 5), show=False).save(tmp_file(name))
    plot_implicit(And(y > exp(x), y > x + 2), show=False).save(tmp_file(name))
    plot_implicit(Or(y > x, y > -x), show=False).save(tmp_file(name))
    plot_implicit(x**2 - 1, (x, -5, 5), show=False).save(tmp_file(name))
    plot_implicit(x**2 - 1, show=False).save(tmp_file(name))
    plot_implicit(y > x, depth=-5, show=False).save(tmp_file(name))
    plot_implicit(y > x, depth=5, show=False).save(tmp_file(name))
    plot_implicit(y > cos(x), adaptive=False, show=False).save(tmp_file(name))
    plot_implicit(y < cos(x), adaptive=False, show=False).save(tmp_file(name))
    plot_implicit(y - cos(pi / x), show=False).save(tmp_file(name))

    pytest.raises(ValueError, lambda: plot_implicit(y > x, (x, -1, 1, 2)))

    # issue sympy/sympy#17719
    plot_implicit(((x - 1)**2 + y**2 < 2) ^ ((x + 1)**2 + y**2 < 2),
                  show=False).save(tmp_file(name))


def test_line_color():
    x, y = symbols('x, y')
    p = plot_implicit(x**2 + y**2 - 1, line_color='green', show=False)
    assert p._series[0].line_color == 'green'
    p = plot_implicit(x**2 + y**2 - 1, line_color='r', show=False)
    assert p._series[0].line_color == 'r'


def test_matplotlib():
    plot_and_save('test')


@pytest.mark.xfail
def test_matplotlib2():
    name = 'test2'

    plot_implicit(And(y > cos(x), Or(y > x, Eq(y, x))), show=False).save(tmp_file(name))

    # Test plots which cannot be rendered using the adaptive algorithm
    # TODO: catch the warning.
    plot_implicit(Eq(y, re(cos(x) + I*sin(x))), show=False).save(tmp_file(name))

    with warnings.catch_warnings(record=True) as w:
        plot_implicit(x**2 - 1, legend='An implicit plot', show=False).save(tmp_file(name))
        assert len(w) == 1
        assert issubclass(w[-1].category, UserWarning)
        assert 'No labelled objects found' in str(w[0].message)
