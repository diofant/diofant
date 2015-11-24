Plotting
========

.. module:: sympy.plotting.plot

Introduction
------------

The plotting module allows you to make 2-dimensional and 3-dimensional plots.
Presently the plots are rendered using ``matplotlib`` as a backend.

The plotting module has the following functions:

* plot: Plots 2D line plots.
* plot_parametric: Plots 2D parametric plots.
* plot_implicit: Plots 2D implicit and region plots.
* plot3d: Plots 3D plots of functions in two variables.
* plot3d_parametric_line: Plots 3D line plots, defined by a parameter.
* plot3d_parametric_surface: Plots 3D parametric surface plots.

The above functions are only for convenience and ease of use. It is possible to
plot any plot by passing the corresponding ``Series`` class to ``Plot`` as
argument.

Plot Class
----------

.. autoclass:: sympy.plotting.plot.Plot
   :members:

Plotting Function Reference
---------------------------

.. autofunction:: plot

.. autofunction:: plot_parametric

.. autofunction:: plot3d

.. autofunction:: plot3d_parametric_line

.. autofunction:: plot3d_parametric_surface

.. autofunction:: sympy.plotting.plot_implicit.plot_implicit

Series Classes
--------------

.. autoclass:: sympy.plotting.plot.BaseSeries
   :members:

.. autoclass:: sympy.plotting.plot.Line2DBaseSeries
   :members:

.. autoclass:: sympy.plotting.plot.LineOver1DRangeSeries
   :members:

.. autoclass:: sympy.plotting.plot.Parametric2DLineSeries
   :members:

.. autoclass:: sympy.plotting.plot.Line3DBaseSeries
   :members:

.. autoclass:: sympy.plotting.plot.Parametric3DLineSeries
   :members:

.. autoclass:: sympy.plotting.plot.SurfaceBaseSeries
   :members:

.. autoclass:: sympy.plotting.plot.SurfaceOver2DRangeSeries
   :members:

.. autoclass:: sympy.plotting.plot.ParametricSurfaceSeries
   :members:

.. autoclass:: sympy.plotting.plot_implicit.ImplicitSeries
   :members:
