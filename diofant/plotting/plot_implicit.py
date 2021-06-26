"""Implicit plotting module for Diofant

The module implements a data series called ImplicitSeries which is used by
``Plot`` class to plot implicit plots for different backends.

See Also
========

diofant.plotting.plot

"""

from ..core import Dummy, Eq, Symbol, Tuple
from ..core.relational import Equality, GreaterThan, LessThan, Relational
from ..core.sympify import sympify
from ..external import import_module
from ..logic.boolalg import BooleanFunction
from ..polys.polyutils import _sort_gens
from ..utilities import flatten, lambdify
from ..utilities.decorator import doctest_depends_on
from .plot import BaseSeries, Plot


class ImplicitSeries(BaseSeries):
    """Representation for Implicit plot."""

    is_implicit = True

    def __init__(self, expr, var_start_end_x, var_start_end_y,
                 has_equality, use_interval_math, depth, nb_of_points,
                 line_color):
        super().__init__()
        self.expr = sympify(expr)
        self.var_x = sympify(var_start_end_x[0])
        self.start_x = float(var_start_end_x[1])
        self.end_x = float(var_start_end_x[2])
        self.var_y = sympify(var_start_end_y[0])
        self.start_y = float(var_start_end_y[1])
        self.end_y = float(var_start_end_y[2])
        self.get_points = self.get_raster

        # If the expression has equality, i.e. Eq, Greaterthan, LessThan.
        self.has_equality = has_equality

        self.nb_of_points = nb_of_points
        self.use_interval_math = use_interval_math
        self.depth = 4 + depth
        self.line_color = line_color

    def get_raster(self):
        return self._get_meshes_grid()

    def _get_meshes_grid(self):
        """Generates the mesh for generating a contour.

        In the case of equality, ``contour`` function of matplotlib can
        be used. In other cases, matplotlib's ``contourf`` is used.

        """
        equal = False
        if isinstance(self.expr, Equality):
            expr = self.expr.lhs - self.expr.rhs
            equal = True
        elif self.expr.has(Equality):  # pragma: no cover
            raise NotImplementedError('The expression is not supported for '
                                      'plotting in uniform meshed plot.')
        else:
            expr = self.expr
        np = import_module('numpy')
        xarray = np.linspace(self.start_x, self.end_x, self.nb_of_points)
        yarray = np.linspace(self.start_y, self.end_y, self.nb_of_points)
        x_grid, y_grid = np.meshgrid(xarray, yarray)

        func = lambdify((self.var_x, self.var_y), expr, 'numpy')
        z_grid = func(x_grid, y_grid)
        z_grid[np.ma.where(z_grid < 0)] = -1
        z_grid[np.ma.where(z_grid > 0)] = 1
        if equal:
            return xarray, yarray, z_grid, 'contour'
        else:
            return xarray, yarray, z_grid, 'contourf'


@doctest_depends_on(modules=('matplotlib',))
def plot_implicit(expr, x_var=None, y_var=None, **kwargs):
    """A plot function to plot implicit equations / inequalities.

    Parameters
    ==========

    ``expr`` : Expr
        The equation / inequality that is to be plotted.
    ``x_var`` : symbol or tuple, optional
        symbol to plot on x-axis or tuple giving symbol and range
        as ``(symbol, xmin, xmax)``
    ``y_var`` : symbol or tuple, optional
        symbol to plot on y-axis or tuple giving symbol and range
        as ``(symbol, ymin, ymax)``

    If neither ``x_var`` nor ``y_var`` are given then the free symbols in the
    expression will be assigned in the order they are sorted.

    The following keyword arguments can also be used:

    ``adaptive`` : Boolean, optional
        The default value is set to True. It has to be
        set to False if you want to use a mesh grid.

    ``depth`` : integer
        The depth of recursion for adaptive mesh grid.
        Default value is 0. Takes value in the range (0, 4).

    ``points`` : integer
        The number of points if adaptive mesh grid is not
        used. Default value is 200.

    ``title`` : str
        The title for the plot.

    ``xlabel`` : str
        The label for the x-axis

    ``ylabel`` : string
        The label for the y-axis

    Aesthetics options:

    ``line_color`` : float or str
        Specifies the color for the plot.  See ``Plot`` to see how to
        set color for the plots.

    plot_implicit, by default, uses interval arithmetic to plot functions. If
    the expression cannot be plotted using interval arithmetic, it defaults to
    a generating a contour using a mesh grid of fixed number of points. By
    setting adaptive to False, you can force plot_implicit to use the mesh
    grid. The mesh grid method can be effective when adaptive plotting using
    interval arithmetic, fails to plot with small line width.

    Examples
    ========

    Plot expressions:

    Without any ranges for the symbols in the expression

    >>> p1 = plot_implicit(Eq(x**2 + y**2, 5))

    With the range for the symbols

    >>> p2 = plot_implicit(Eq(x**2 + y**2, 3),
    ...                    (x, -3, 3), (y, -3, 3))

    With depth of recursion as argument.

    >>> p3 = plot_implicit(Eq(x**2 + y**2, 5),
    ...                    (x, -4, 4), (y, -4, 4), depth=2)

    Using mesh grid and not using adaptive meshing.

    >>> p4 = plot_implicit(Eq(x**2 + y**2, 5),
    ...                    (x, -5, 5), (y, -2, 2), adaptive=False)

    Using mesh grid with number of points as input.

    >>> p5 = plot_implicit(Eq(x**2 + y**2, 5),
    ...                    (x, -5, 5), (y, -2, 2),
    ...                    adaptive=False, points=400)

    Plotting regions.

    >>> p6 = plot_implicit(y > x**2)

    Plotting Using boolean conjunctions.

    >>> p7 = plot_implicit(And(y > x, y > -x))

    When plotting an expression with a single variable (y - 1, for example),
    specify the x or the y variable explicitly:

    >>> p8 = plot_implicit(y - 1, y_var=y)
    >>> p9 = plot_implicit(x - 1, x_var=x)

    """
    # Represents whether the expression contains an Equality,
    # GreaterThan or LessThan
    has_equality = False

    def arg_expand(bool_expr):
        """Recursively expands the arguments of an Boolean Function."""
        for arg in bool_expr.args:
            if isinstance(arg, BooleanFunction):
                arg_expand(arg)
            elif isinstance(arg, Relational):
                arg_list.append(arg)

    arg_list = []
    if isinstance(expr, BooleanFunction):
        arg_expand(expr)

    # Check whether there is an equality in the expression provided.
        if any(isinstance(e, (Equality, GreaterThan, LessThan))
               for e in arg_list):
            has_equality = True

    elif not isinstance(expr, Relational):
        expr = Eq(expr, 0)
        has_equality = True
    elif isinstance(expr, (Equality, GreaterThan, LessThan)):
        has_equality = True

    xyvar = [i for i in (x_var, y_var) if i is not None]
    free_symbols = expr.free_symbols
    range_symbols = Tuple(*flatten(xyvar)).free_symbols
    undeclared = free_symbols - range_symbols
    if len(free_symbols & range_symbols) > 2:  # pragma: no cover
        raise NotImplementedError('Implicit plotting is not implemented for '
                                  'more than 2 variables')

    # Create default ranges if the range is not provided.
    default_range = Tuple(-5, 5)

    def _range_tuple(s):
        if isinstance(s, (Dummy, Symbol)):
            return Tuple(s) + default_range
        if len(s) == 3:
            return Tuple(*s)
        raise ValueError(f'symbol or `(symbol, min, max)` expected but got {s!s}')

    if len(xyvar) == 0:
        xyvar = list(_sort_gens(free_symbols))
    var_start_end_x = _range_tuple(xyvar[0])
    x = var_start_end_x[0]
    if len(xyvar) != 2:
        if x in undeclared or not undeclared:
            xyvar.append(Dummy(f'f({x.name})'))
        else:
            xyvar.append(undeclared.pop())
    var_start_end_y = _range_tuple(xyvar[1])

    use_interval = kwargs.pop('adaptive', False)
    nb_of_points = kwargs.pop('points', 300)
    depth = kwargs.pop('depth', 0)
    line_color = kwargs.pop('line_color', 'blue')
    # Check whether the depth is greater than 4 or less than 0.
    if depth > 4:
        depth = 4
    elif depth < 0:
        depth = 0

    series_argument = ImplicitSeries(expr, var_start_end_x, var_start_end_y,
                                     has_equality, use_interval, depth,
                                     nb_of_points, line_color)
    show = kwargs.pop('show', True)

    # set the x and y limits
    kwargs['xlim'] = tuple(float(x) for x in var_start_end_x[1:])
    kwargs['ylim'] = tuple(float(y) for y in var_start_end_y[1:])
    # set the x and y labels
    kwargs.setdefault('xlabel', var_start_end_x[0].name)
    kwargs.setdefault('ylabel', var_start_end_y[0].name)
    p = Plot(series_argument, **kwargs)
    if show:
        p.show()
    return p
