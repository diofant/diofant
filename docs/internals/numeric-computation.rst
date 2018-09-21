Numeric Computation
===================

Symbolic computer algebra systems like Diofant facilitate the construction and
manipulation of mathematical expressions.  Unfortunately when it comes time
to evaluate these expressions on numerical data, symbolic systems often have
poor performance.

Fortunately Diofant offers a number of easy-to-use hooks into other numeric
systems, allowing you to create mathematical expressions in Diofant and then
ship them off to the numeric system of your choice.  This page documents many
of the options available including the ``math`` library, the popular array
computing package ``numpy``, code generation in ``Fortran`` or ``C``.


Subs/evalf
----------

Subs is the slowest but simplest option.  It runs at Diofant speeds.
The ``.subs(...).evalf()`` method can substitute a numeric value
for a symbolic one and then evaluate the result within Diofant.


    >>> expr = sin(x)/x
    >>> expr.evalf(subs={x: 3.14}, strict=False)
    0.000507214304613640

This method is slow.  You should use this method production only if performance
is not an issue.  You can expect ``.subs`` to take tens of microseconds. It
can be useful while prototyping or if you just want to see a value once.


Lambdify
--------

The ``lambdify`` function translates Diofant expressions into Python functions,
leveraging a variety of numerical libraries.  It is used as follows:

    >>> expr = sin(x)/x
    >>> f = lambdify(x, expr)
    >>> f(3.14)
    0.000507214304614

Here lambdify makes a function that computes ``f(x) = sin(x)/x``.  By default
lambdify relies on implementations in the ``math`` standard library. This
numerical evaluation takes on the order of hundreds of nanoseconds, roughly two
orders of magnitude faster than the ``.subs`` method.  This is the speed
difference between Diofant and raw Python.

Lambdify can leverage a variety of numerical backends.  By default it uses the
``math`` library.  However it also supports ``mpmath`` and most notably,
``numpy``.  Using the ``numpy`` library gives the generated function access to
powerful vectorized ufuncs that are backed by compiled C code.

    >>> expr = sin(x)/x
    >>> f = lambdify(x, expr, "numpy")

    >>> import numpy
    >>> data = numpy.linspace(1, 10, 10000)
    >>> pprint(f(data))
    [ 0.84147098  0.84119981  0.84092844 ..., -0.05426074 -0.05433146
                               -0.05440211]

If you have array-based data this can confer a considerable speedup, on the
order of 10 nano-seconds per element. Unfortunately numpy incurs some start-up
time and introduces an overhead of a few microseconds.

uFuncify
--------

While NumPy operations are very efficient for vectorized data they sometimes
incur unnecessary costs when chained together. Consider the following operation

.. code:: python

    x = get_numpy_array(...)
    y = sin(x)/x

The operators ``sin`` and ``/`` call routines that execute tight for loops in
``C``. The resulting computation looks something like this

.. code:: c

    for(int i = 0; i < n; i++)
    {
        temp[i] = sin(x[i]);
    }
    for(int i = i; i < n; i++)
    {
        y[i] = temp[i] / x[i];
    }

This is slightly sub-optimal because

1.  We allocate an extra ``temp`` array
2.  We walk over ``x`` memory twice when once would have been sufficient

A better solution would fuse both element-wise operations into a single for loop

.. code:: c

    for(int i = i; i < n; i++)
    {
        y[i] = sin(x[i]) / x[i];
    }

Statically compiled projects like NumPy are unable to take advantage of such
optimizations. Fortunately, Diofant is able to generate efficient low-level C
or Fortran code. It can then depend on projects like ``Cython`` or ``f2py`` to
compile and reconnect that code back up to Python. Fortunately this process is
well automated and a Diofant user wishing to make use of this code generation
should call the ``ufuncify`` function

    >>> expr = sin(x)/x

    >>> from diofant.utilities.autowrap import ufuncify
    >>> f = ufuncify((x,), expr)

This function ``f`` consumes and returns a NumPy array. Generally ``ufuncify``
performs at least as well as ``lambdify``. If the expression is complicated
then ``ufuncify`` often significantly outperforms the NumPy backed solution.
Jensen has a good `blog post <https://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/>`_
on this topic.

So Which Should I Use?
----------------------

The options here were listed in order from slowest and least dependencies to
fastest and most dependencies.

+-----------------+-------+-----------------------------+---------------+
| Tool            | Speed | Qualities                   | Dependencies  |
+=================+=======+=============================+===============+
| subs/evalf      | 50us  | Simple                      | None          |
+-----------------+-------+-----------------------------+---------------+
| lambdify        | 1us   | Scalar functions            | math          |
+-----------------+-------+-----------------------------+---------------+
| lambdify-numpy  | 10ns  | Vector functions            | numpy         |
+-----------------+-------+-----------------------------+---------------+
| ufuncify        | 10ns  | Complex vector expressions  | f2py, Cython  |
+-----------------+-------+-----------------------------+---------------+
