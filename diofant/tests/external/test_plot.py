import errno
import functools
import os
import sys
import tempfile

import pytest

from diofant import (And, I, Integral, LambertW, Piecewise, cos, exp_polar,
                     log, meijerg, oo, pi, real_root, sin, sqrt, summation)
from diofant.abc import x, y, z
from diofant.plotting import (plot, plot3d, plot3d_parametric_line,
                              plot3d_parametric_surface, plot_parametric)
from diofant.plotting.plot import unset_show


__all__ = ()

matplotlib = pytest.importorskip('matplotlib', minversion='1.1.0')


class MockPrint:

    def write(self, s):
        pass

    def flush(self):
        pass


def disable_print(func, *args, **kwargs):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        sys.stdout = MockPrint()
        func(*args, **kwargs)
        sys.stdout = sys.__stdout__
    return wrapper


unset_show()


# XXX: We could implement this as a context manager instead
class TmpFileManager:

    tmp_files = []

    @classmethod
    def tmp_file(cls, name=''):
        cls.tmp_files.append(tempfile.NamedTemporaryFile(prefix=name, suffix='.png').name)
        return cls.tmp_files[-1]

    @classmethod
    def cleanup(cls):
        for f in cls.tmp_files:
            try:
                os.remove(f)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise


def test_matplotlib_intro():
    """Examples from the 'introduction' notebook."""
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        p = plot(x)
        p = plot(x*sin(x), x*cos(x))
        p.extend(p)
        p[0].line_color = lambda a: a
        p[1].line_color = 'b'
        p.title = 'Big title'
        p.xlabel = 'the x axis'
        p[1].label = 'straight line'
        p.legend = True
        p.aspect_ratio = (1, 1)
        p.xlim = (-15, 20)
        p.save(tmp_file('%s_basic_options_and_colors' % name))

        p.extend(plot(x + 1))
        p.append(plot(x + 3, x**2)[1])
        p.save(tmp_file('%s_plot_extend_append' % name))

        p[2] = plot(x**2, (x, -2, 3))
        p.save(tmp_file('%s_plot_setitem' % name))
        del p

        p = plot(sin(x), (x, -2*pi, 4*pi))
        p.save(tmp_file('%s_line_explicit' % name))
        del p

        p = plot(sin(x))
        p.save(tmp_file('%s_line_default_range' % name))
        del p

        p = plot((x**2, (x, -5, 5)), (x**3, (x, -3, 3)))
        p.save(tmp_file('%s_line_multiple_range' % name))
        del p

        pytest.raises(ValueError, lambda: plot(x, y))

        # parametric 2d plots.

        # Single plot with default range.
        plot_parametric(sin(x), cos(x)).save(tmp_file())

        # Single plot with range.
        p = plot_parametric(sin(x), cos(x), (x, -5, 5))
        p.save(tmp_file('%s_parametric_range' % name))
        del p

        # Multiple plots with same range.
        p = plot_parametric((sin(x), cos(x)), (x, sin(x)))
        p.save(tmp_file('%s_parametric_multiple' % name))
        del p

        # Multiple plots with different ranges.
        p = plot_parametric((sin(x), cos(x), (x, -3, 3)), (x, sin(x), (x, -5, 5)))
        p.save(tmp_file('%s_parametric_multiple_ranges' % name))
        del p

        # depth of recursion specified.
        p = plot_parametric(x, sin(x), depth=13)
        p.save(tmp_file('%s_recursion_depth' % name))
        del p

        # No adaptive sampling.
        p = plot_parametric(cos(x), sin(x), adaptive=False, nb_of_points=500)
        p.save(tmp_file('%s_adaptive' % name))
        del p

        # 3d parametric plots
        p = plot3d_parametric_line(sin(x), cos(x), x)
        p.save(tmp_file('%s_3d_line' % name))
        del p

        p = plot3d_parametric_line(
            (sin(x), cos(x), x, (x, -5, 5)), (cos(x), sin(x), x, (x, -3, 3)))
        p.save(tmp_file('%s_3d_line_multiple' % name))
        del p

        p = plot3d_parametric_line(sin(x), cos(x), x, nb_of_points=30)
        p.save(tmp_file('%s_3d_line_points' % name))
        del p

        # 3d surface single plot.
        p = plot3d(x * y)
        p.save(tmp_file('%s_surface' % name))
        del p

        # Multiple 3D plots with same range.
        p = plot3d(-x * y, x * y, (x, -5, 5))
        p.save(tmp_file('%s_surface_multiple' % name))
        del p

        # Multiple 3D plots with different ranges.
        p = plot3d(
            (x * y, (x, -3, 3), (y, -3, 3)), (-x * y, (x, -3, 3), (y, -3, 3)))
        p.save(tmp_file('%s_surface_multiple_ranges' % name))
        del p

        # Single Parametric 3D plot
        p = plot3d_parametric_surface(sin(x + y), cos(x - y), x - y)
        p.save(tmp_file('%s_parametric_surface' % name))
        del p

        # Multiple Parametric 3D plots.
        p = plot3d_parametric_surface(
            (x*sin(z), x*cos(z), z, (x, -5, 5), (z, -5, 5)),
            (sin(x + y), cos(x - y), x - y, (x, -5, 5), (y, -5, 5)))
        p.save(tmp_file('%s_parametric_surface' % name))
        del p
    finally:
        TmpFileManager.cleanup()


def test_matplotlib_colors():
    """Examples from the 'colors' notebook."""
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        p = plot(sin(x))
        p[0].line_color = lambda a: a
        p.save(tmp_file('%s_colors_line_arity1' % name))

        p[0].line_color = lambda a, b: b
        p.save(tmp_file('%s_colors_line_arity2' % name))

        p = plot(x*sin(x), x*cos(x), (x, 0, 10))
        p[0].line_color = lambda a: a
        p.save(tmp_file('%s_colors_param_line_arity1' % name))

        p[0].line_color = lambda a, b: a
        p.save(tmp_file('%s_colors_param_line_arity2a' % name))

        p[0].line_color = lambda a, b: b
        p.save(tmp_file('%s_colors_param_line_arity2b' % name))
        del p

        p = plot3d_parametric_line(sin(x) + 0.1*sin(x)*cos(7*x),
                                   cos(x) + 0.1*cos(x)*cos(7*x), 0.1*sin(7*x), (x, 0, 2*pi))
        p[0].line_color = lambda a: sin(4*a)
        p.save(tmp_file('%s_colors_3d_line_arity1' % name))
        p[0].line_color = lambda a, b: b
        p.save(tmp_file('%s_colors_3d_line_arity2' % name))
        p[0].line_color = lambda a, b, c: c
        p.save(tmp_file('%s_colors_3d_line_arity3' % name))
        del p

        p = plot3d(sin(x)*y, (x, 0, 6*pi), (y, -5, 5))
        p[0].surface_color = lambda a: a
        p.save(tmp_file('%s_colors_surface_arity1' % name))
        p[0].surface_color = lambda a, b: b
        p.save(tmp_file('%s_colors_surface_arity2' % name))
        p[0].surface_color = lambda a, b, c: c
        p.save(tmp_file('%s_colors_surface_arity3a' % name))
        p[0].surface_color = lambda a, b, c: sqrt((a - 3*pi)**2 + b**2)
        p.save(tmp_file('%s_colors_surface_arity3b' % name))
        del p

        p = plot3d_parametric_surface(x * cos(4 * y), x * sin(4 * y), y,
                                      (x, -1, 1), (y, -1, 1))
        p[0].surface_color = lambda a: a
        p.save(tmp_file('%s_colors_param_surf_arity1' % name))
        p[0].surface_color = lambda a, b: a*b
        p.save(tmp_file('%s_colors_param_surf_arity2' % name))
        p[0].surface_color = lambda a, b, c: sqrt(a**2 + b**2 + c**2)
        p.save(tmp_file('%s_colors_param_surf_arity3' % name))
        del p
    finally:
        TmpFileManager.cleanup()


@pytest.mark.xfail
def test_matplotlib_advanced():
    """Examples from the 'advanced' notebook."""
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        s = summation(1/x**y, (x, 1, oo))
        p = plot(s, (y, 2, 10))
        p.save(tmp_file('%s_advanced_inf_sum' % name))

        p = plot(summation(1/x, (x, 1, y)), (y, 2, 10), show=False)
        p[0].only_integers = True
        p[0].steps = True
        p.save(tmp_file('%s_advanced_fin_sum' % name))

        ###
        # Test expressions that can not be translated to np and
        # generate complex results.
        ###
        plot(sin(x) + I*cos(x)).save(tmp_file())
        plot(sqrt(sqrt(-x))).save(tmp_file())
        plot(LambertW(x)).save(tmp_file())
        plot(sqrt(LambertW(x))).save(tmp_file())

        # Characteristic function of a StudentT distribution with nu=10
        plot((meijerg(((1 / 2,), ()), ((5, 0, 1 / 2), ()),
                      5 * x**2 * exp_polar(-I*pi)/2)
              + meijerg(((1/2,), ()), ((5, 0, 1/2), ()),
                        5*x**2 * exp_polar(I*pi)/2)) / (48 * pi), (x, 1e-6, 1e-2)).save(tmp_file())
    finally:
        TmpFileManager.cleanup()


@pytest.mark.xfail
def test_matplotlib_advanced_2():
    """More examples from the 'advanced' notebook."""
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        i = Integral(log((sin(x)**2 + 1)*sqrt(x**2 + 1)), (x, 0, y))
        p = plot(i, (y, 1, 5))
        p.save(tmp_file('%s_advanced_integral' % name))
    finally:
        TmpFileManager.cleanup()


@disable_print
def test_append_sympyissue_7140():
    p1 = plot(x)
    p2 = plot(x**2)

    # append a series
    p2.append(p1[0])
    assert len(p2._series) == 2

    with pytest.raises(TypeError):
        p1.append(p2)

    with pytest.raises(TypeError):
        p1.append(p2._series)


@disable_print
def test_sympyissue_11461():
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        p = plot(real_root((log(x/(x-2))), 3), (x, 3, 4))
        p.save(tmp_file('%s_11461' % name))
        del p
    finally:
        TmpFileManager.cleanup()


@disable_print
def test_sympyissue_10925():
    try:
        name = 'test'
        tmp_file = TmpFileManager.tmp_file

        f = Piecewise((-1, x < -1), (x, And(-1 <= x, x < 0)),
                      (x**2, And(0 <= x, x < 1)), (x**3, True))
        p = plot(f, (x, -3, 3))
        p.save(tmp_file('%s_10925' % name))
        del p
    finally:
        TmpFileManager.cleanup()
