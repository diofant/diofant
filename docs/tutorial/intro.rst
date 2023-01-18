==============
 Introduction
==============

Symbolic computation deals with the computation of mathematical
objects symbolically.  This means that the mathematical objects are
represented exactly, not approximately, and mathematical expressions
with unevaluated variables are left in symbolic form.

Let's take an example.  Start the Python interpreter::

   python

Say we wanted to use the built-in Python functions to compute square
roots.  We might do something like this

   >>> import math
   >>> math.sqrt(9)
   3.0

Here we got the exact answer --- 9 is a perfect square --- but usually
it will be an approximate result

   >>> math.sqrt(8)
   2.8284271247461903

This is where symbolic computation first comes in: with a symbolic
computation system like Diofant, square roots of numbers that are not
perfect squares are left unevaluated by default

   >>> import diofant
   >>> diofant.sqrt(3)
   sqrt(3)

Furthermore --- and this is where we start to see the real power of
symbolic computation --- results can be symbolically simplified.

   >>> diofant.sqrt(8)
   2*sqrt(2)

Yet we can also approximate this number with any precision

   >>> _.evalf(20)
   2.8284271247461900976

The above example starts to show how we can manipulate irrational
numbers exactly using Diofant.  Now we introduce symbols.

Let us define a symbolic expression, representing the mathematical
expression `x + 2y`.

   >>> x, y = diofant.symbols('x y')
   >>> expr = x + 2*y
   >>> expr
   x + 2*y

.. note::

   Unlike many symbolic manipulation systems you may have used, in
   Diofant symbols are not defined automatically.  To define symbols
   (instances of :class:`~diofant.core.symbol.Symbol`) you may
   use :func:`~diofant.core.symbol.symbols`.

Note that we wrote ``x + 2*y``, using Python's mathematical syntax,
just as we would if ``x`` and ``y`` were ordinary Python variables.
But in this case, instead of evaluating to something, the expression
remains as just ``x + 2*y``.  Now let us play around with it:

   >>> expr + 1
   x + 2*y + 1
   >>> expr - x
   2*y

Notice something in the above example.  When we typed ``expr - x``, we
did not get ``x + 2*y - x``, but rather just ``2*y``.  The ``x`` and
the ``-x`` automatically canceled one another.  This is similar to how
``sqrt(8)`` automatically turned into ``2*sqrt(2)`` above.

.. tip::

   Use :func:`~diofant.core.evaluate.evaluate` context or ``evaluate``
   flag to prevent automatic evaluation, for example:

       >>> diofant.sqrt(8, evaluate=False)
       sqrt(8)
       >>> _.doit()
       2*sqrt(2)

This isn't always the case in Diofant, however:

   >>> x*expr
   x*(x + 2*y)

Here, we might have expected `x(x + 2y)` to transform into `x^2 +
2xy`, but instead we see that the expression was left alone.  This is
a common theme in Diofant.  Aside from obvious simplifications like
`x - x = 0` and `\sqrt{8} = 2\sqrt{2}`, most simplifications are not
performed automatically.  This is because we might prefer the factored
form `x(x + 2y)`, or we might prefer the expanded form `x^2 + 2xy` ---
both forms are useful in different circumstances.  In Diofant, there
are functions to go from one form to the other

   >>> diofant.expand(x*expr)
   x**2 + 2*x*y
   >>> diofant.factor(_)
   x*(x + 2*y)

The real power of a symbolic computation system (which by the way, are
also often called computer algebra systems, or just CASs) such as
Diofant is the ability to do all sorts of computations symbolically:
simplify expressions, compute derivatives, integrals, and limits,
solve equations, work with matrices, and much more.  Diofant includes
modules for printing (like 2D pretty printed output of math
formulas, or `\LaTeX`), code generation, combinatorics,
number theory, logic, and more.  Here is a small sampling of the sort
of symbolic power Diofant is capable of, to whet your appetite.

.. note::

   From here on in this tutorial we assume that these statements were
   executed:

      >>> from diofant import *
      >>> a, b, c, d, t, x, y, z = symbols('a:d t x:z')
      >>> init_printing(pretty_print=True, use_unicode=True)

   Last one will make all further examples pretty print with unicode
   characters.

   ``import *`` has been used here to aid the readability of the
   tutorial, but is best to avoid such wildcard import statements in
   production code, as they make it unclear which names are present in
   the namespace.

Take the derivative of `\sin{(x)}e^x`.

   >>> diff(sin(x)*exp(x))
    x           x
   ℯ ⋅sin(x) + ℯ ⋅cos(x)

Compute `\int(e^x\sin{(x)} + e^x\cos{(x)})\,dx`.

   >>> integrate(exp(x)*sin(x) + exp(x)*cos(x))
    x
   ℯ ⋅sin(x)

Compute `\int_{-\infty}^\infty \sin{(x^2)}\,dx`.

   >>> integrate(sin(x**2), (x, -oo, oo))
     ___   ___
   ╲╱ 2 ⋅╲╱ π
   ───────────
        2

Find `\lim_{x\to 0^+}\frac{\sin{(x)}}{x}`.

   >>> limit(sin(x)/x, x, 0)
   1

Solve `x^2 - 2 = 0`.

   >>> solve(x**2 - 2, x)
   ⎡⎧      ___⎫  ⎧     ___⎫⎤
   ⎢⎨x: -╲╱ 2 ⎬, ⎨x: ╲╱ 2 ⎬⎥
   ⎣⎩         ⎭  ⎩        ⎭⎦

Solve the differential equation `f'' - f = e^x`.

   >>> f = symbols('f', cls=Function)
   >>> dsolve(Eq(f(x).diff((x, 2)) - f(x), exp(x)))
           x ⎛     x⎞    -x
   f(x) = ℯ ⋅⎜C₂ + ─⎟ + ℯ  ⋅C₁
             ⎝     2⎠

Find the eigenvalues of `\left[\begin{smallmatrix}1 & 2\\2 &
2\end{smallmatrix}\right]`.

   >>> Matrix([[1, 2], [2, 2]]).eigenvals()
   ⎧      ____         ____       ⎫
   ⎪3   ╲╱ 17        ╲╱ 17    3   ⎪
   ⎨─ + ──────: 1, - ────── + ─: 1⎬
   ⎪2     2            2      2   ⎪
   ⎩                              ⎭

Rewrite the Bessel function `J_y\left(z\right)` in terms of the
spherical Bessel function `j_y(z)`.

   >>> besselj(y, z).rewrite(jn)
     ___   ___
   ╲╱ 2 ⋅╲╱ z ⋅jn(y - 1/2, z)
   ──────────────────────────
               ___
             ╲╱ π

Print `\int_{0}^{\pi} \cos^{2}{\left (x \right )}\, dx` using `\LaTeX`.

   >>> latex(Integral(cos(x)**2, (x, 0, pi)))
   '\\int_{0}^{\\pi} \\cos^{2}{\\left (x \\right )}\\, dx'
