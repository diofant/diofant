.. _diophantine-docs:

Diophantine
===========

.. automodule:: diofant.solvers.diophantine

Diophantine equations
---------------------

The word "Diophantine" comes with the name Diophantus, a mathematician lived
in the great city of Alexandria sometime around 250 AD. Often referred to as
the "father of Algebra", Diophantus in his famous work "Arithmetica"
presented 150 problems that marked the early beginnings of number theory, the
field of study about integers and their properties. Diophantine equations play
a central and an important part in number theory.

We call a "Diophantine equation" to an equation of the form,
`f(x_1, x_2, \ldots x_n) = 0` where `n \geq 2` and `x_1, x_2, \ldots x_n` are
integer variables. If we can find `n` integers `a_1, a_2, \ldots a_n` such that
`x_1 = a_1, x_2 = a_2, \ldots x_n = a_n` satisfies the above equation, we say
that the equation is solvable.

Currently, following five types of Diophantine equations can be solved using
:py:meth:`~diofant.solvers.diophantine.diophantine` and other helper functions of
the Diophantine module.

- Linear Diophantine equations: `a_1x_1 + a_2x_2 + \ldots + a_nx_n = b`.
- General binary quadratic equation: `ax^2 + bxy + cy^2 + dx + ey + f = 0`
- Homogeneous ternary quadratic equation: `ax^2 + by^2 + cz^2 + dxy + eyz + fzx = 0`
- Extended Pythagorean equation: `a_{1}x_{1}^2 + a_{2}x_{2}^2 + \ldots + a_{n}x_{n}^2 = a_{n+1}x_{n+1}^2`
- General sum of squares: `x_{1}^2 + x_{2}^2 + \ldots + x_{n}^2 = k`

Module structure
----------------

This module contains :py:meth:`~diofant.solvers.diophantine.diophantine` and
helper functions that are needed to solve certain Diophantine equations. It's
structured in the following manner.

- :py:meth:`~diofant.solvers.diophantine.diophantine`

  - :py:meth:`~diofant.solvers.diophantine.diop_solve`

    - :py:meth:`~diofant.solvers.diophantine.classify_diop`
    - :py:meth:`~diofant.solvers.diophantine.diop_linear`
    - :py:meth:`~diofant.solvers.diophantine.diop_quadratic`
    - :py:meth:`~diofant.solvers.diophantine.diop_ternary_quadratic`
    - :py:meth:`~diofant.solvers.diophantine.diop_ternary_quadratic_normal`
    - :py:meth:`~diofant.solvers.diophantine.diop_general_pythagorean`
    - :py:meth:`~diofant.solvers.diophantine.diop_general_sum_of_squares`
    - :py:meth:`~diofant.solvers.diophantine.diop_general_sum_of_even_powers`

  - :py:meth:`~diofant.solvers.diophantine.merge_solution`

When an equation is given to :py:meth:`~diofant.solvers.diophantine.diophantine`,
it factors the equation(if possible) and solves the equation given by each
factor by calling :py:meth:`~diofant.solvers.diophantine.diop_solve` separately.
Then all the results are combined using :py:meth:`~diofant.solvers.diophantine.merge_solution`.

:py:meth:`~diofant.solvers.diophantine.diop_solve` internally uses
:py:meth:`~diofant.solvers.diophantine.classify_diop`
to find the type of the equation(and some other details) given to it and then
calls the appropriate solver function based on the type returned. For example,
if :py:meth:`~diofant.solvers.diophantine.classify_diop` returned "linear" as the
type of the equation, then :py:meth:`~diofant.solvers.diophantine.diop_solve`
calls :py:meth:`~diofant.solvers.diophantine.diop_linear` to solve the equation.

Each of the functions, :py:meth:`~diofant.solvers.diophantine.diop_linear`,
:py:meth:`~diofant.solvers.diophantine.diop_quadratic`,
:py:meth:`~diofant.solvers.diophantine.diop_ternary_quadratic`,
:py:meth:`~diofant.solvers.diophantine.diop_general_pythagorean`
and :py:meth:`~diofant.solvers.diophantine.diop_general_sum_of_squares` solves a
specific type of equations and the type can be easily guessed by it's name.

Apart from these functions, there are a considerable number of other functions
in the "Diophantine Module" and all of them are listed under User functions
and Internal functions.

Tutorial
--------

First, let's import the highest API of the Diophantine module.

>>> from diofant.solvers.diophantine import diophantine

Before we start solving the equations, we need to define the variables.

>>> x, y, z, t, p, q = symbols('x, y, z, t, p, q', integer=True)
>>> t1, t2, t3, t4, t5 = symbols('t1:6', integer=True)

Let's start by solving the easiest type of Diophantine equations, i.e. linear
Diophantine equations. Let's solve `2x + 3y = 5`. Note that although we
write the equation in the above form, when we input the equation to any of the
functions in Diophantine module, it needs to be in the form `eq = 0`.

>>> diophantine(2*x + 3*y - 5)
{(3*t_0 - 5, -2*t_0 + 5)}

Note that stepping one more level below the highest API, we can solve the very
same equation by calling :py:meth:`~diofant.solvers.diophantine.diop_solve`.

>>> from diofant.solvers.diophantine import diop_solve
>>> diop_solve(2*x + 3*y - 5)
(3*t_0 - 5, -2*t_0 + 5)

Note that it returns a tuple rather than a set.
:py:meth:`~diofant.solvers.diophantine.diophantine` always return a set of tuples.
But :py:meth:`~diofant.solvers.diophantine.diop_solve` may return a single tuple
or a set of tuples depending on the type of the equation given.

We can also solve this equation by calling :py:meth:`~diofant.solvers.diophantine.diop_linear`,
which is what :py:meth:`~diofant.solvers.diophantine.diop_solve` calls internally.

>>> from diofant.solvers.diophantine import diop_linear
>>> diop_linear(2*x + 3*y - 5)
(3*t_0 - 5, -2*t_0 + 5)

If the given equation has no solutions then the outputs will look like below.

>>> diophantine(2*x + 4*y - 3)
set()
>>> diop_solve(2*x + 4*y - 3)
(None, None)
>>> diop_linear(2*x + 4*y - 3)
(None, None)

Note that except for the highest level API, in case of no solutions, a tuple of
`None` are returned. Size of the tuple is the same as the number of variables.
Also, one can specifically set the parameter to be used in the solutions by
passing a customized parameter. Consider the following example:

>>> diop_solve(2*x + 3*y - 5, m)
(3*m_0 - 5, -2*m_0 + 5)

For linear Diophantine equations, the customized parameter is the prefix used
for each free variable in the solution. Consider the following example:

>>> diop_solve(2*x + 3*y - 5*z + 7, m)
(m_0, m_0 + 5*m_1 - 14, m_0 + 3*m_1 - 7)

In the solution above, m_0 and m_1 are independent free variables.

Please note that for the moment, users can set the parameter only for linear
Diophantine equations and binary quadratic equations.

Let's try solving a binary quadratic equation which is an equation with two
variables and has a degree of two. Before trying to solve these equations, an
idea about various cases associated with the equation would help a lot.
Let us define `\Delta = b^2 - 4ac` w.r.t. the binary quadratic
`ax^2 + bxy + cy^2 + dx + ey + f = 0`.

When `\Delta < 0`, there are either no solutions or only a finite number of solutions.

>>> diophantine(x**2 - 4*x*y + 8*y**2 - 3*x + 7*y - 5)
{(2, 1), (5, 1)}

In the above equation `\Delta = (-4)^2 - 4*1*8 = -16` and hence only a finite
number of solutions exist.

When `\Delta = 0` we might have either no solutions or parameterized solutions.

>>> diophantine(3*x**2 - 6*x*y + 3*y**2 - 3*x + 7*y - 5)
set()
>>> diophantine(x**2 - 4*x*y + 4*y**2 - 3*x + 7*y - 5)
{(-2*t**2 - 7*t + 10, -t**2 - 3*t + 5)}
>>> diophantine(x**2 + 2*x*y + y**2 - 3*x - 3*y)
{(t_0, -t_0), (t_0, -t_0 + 3)}

The most interesting case is when `\Delta > 0` and it is not a perfect square.
In this case, the equation has either no solutions or an infinite number of
solutions. Consider the below cases where `\Delta = 8`.

>>> diophantine(x**2 - 4*x*y + 2*y**2 - 3*x + 7*y - 5)
set()
>>> s = diophantine(x**2 - 2*y**2 - 2*x - 4*y, n)
>>> x_1, y_1 = s.pop()
>>> x_2, y_2 = s.pop()
>>> x_n = (-(-2*sqrt(2) + 3)**n/2 + sqrt(2)*(-2*sqrt(2) + 3)**n/2 -
...        sqrt(2)*(2*sqrt(2) + 3)**n/2 - (2*sqrt(2) + 3)**n/2 + 1)
>>> x_1 == x_n or x_2 == x_n
True
>>> y_n = (-sqrt(2)*(-2*sqrt(2) + 3)**n/4 + (-2*sqrt(2) + 3)**n/2 +
...        sqrt(2)*(2*sqrt(2) + 3)**n/4 + (2*sqrt(2) + 3)**n/2 - 1)
>>> y_1 == y_n or y_2 == y_n
True

Here `n` is an integer. Although x_n and y_n may not look like
integers, substituting in specific values for n (and simplifying) shows that they
are. For example consider the following example where we set n equal to 9.

>>> simplify(x_n.subs({n: 9}))
-9369318

Any binary quadratic of the form `ax^2 + bxy + cy^2 + dx + ey + f = 0` can be
transformed to an equivalent form `X^2 - DY^2 = N`.

>>> from diofant.solvers.diophantine import find_DN, diop_DN, transformation_to_DN
>>> find_DN(x**2 - 3*x*y + y**2 - 7*x + 5*y - 3)
(5, 920)

So, the above equation is equivalent to the equation `X^2 - 5Y^2 = 920` after
a linear transformation. If we want to find the linear transformation, we can
use :py:meth:`~diofant.solvers.diophantine.transformation_to_DN`

>>> A, B = transformation_to_DN(x**2 - 3*x*y + y**2 - 7*x + 5*y - 3)

Here A is a 2 X 2 matrix and B is a 2 X 1 matrix such that the transformation

.. math::

    \begin{bmatrix} X\\Y \end{bmatrix} = A \begin{bmatrix} x\\y \end{bmatrix} + B

gives the equation `X^2 -5Y^2 = 920`. Values of `A` and `B` are as belows.

>>> A
Matrix([
[1/10, 3/10],
[   0,  1/5]])
>>> B
Matrix([
[  1/5],
[-11/5]])

We can solve an equation of the form `X^2 - DY^2 = N` by passing `D` and `N` to
:py:meth:`~diofant.solvers.diophantine.diop_DN`

>>> diop_DN(5, 920)
[]

Unfortunately, our equation has no solution.

Now let's turn to homogeneous ternary quadratic equations. These equations are
of the form `ax^2 + by^2 + cz^2 + dxy + eyz + fzx = 0`. These type of equations
either have infinitely many solutions or no solutions (except the obvious
solution (0, 0, 0))

>>> diophantine(3*x**2 + 4*y**2 - 5*z**2 + 4*x*y + 6*y*z + 7*z*x)
{(0, 0, 0)}
>>> diophantine(3*x**2 + 4*y**2 - 5*z**2 + 4*x*y - 7*y*z + 7*z*x)
{(-16*p**2 + 28*p*q + 20*q**2, 3*p**2 + 38*p*q - 25*q**2, 4*p**2 - 24*p*q + 68*q**2)}

If you are only interested in a base solution rather than the parameterized
general solution (to be more precise, one of the general solutions), you can
use :py:meth:`~diofant.solvers.diophantine.diop_ternary_quadratic`.

>>> from diofant.solvers.diophantine import diop_ternary_quadratic
>>> diop_ternary_quadratic(3*x**2 + 4*y**2 - 5*z**2 + 4*x*y - 7*y*z + 7*z*x)
(-4, 5, 1)

:py:meth:`~diofant.solvers.diophantine.diop_ternary_quadratic` first converts the
given equation to an equivalent equation of the form `w^2 = AX^2 + BY^2` and
then it uses :py:meth:`~diofant.solvers.diophantine.descent` to solve the latter
equation. You can refer to the docs of
:py:meth:`~diofant.solvers.diophantine.transformation_to_normal` to find more on
this. The equation `w^2 = AX^2 + BY^2` can be solved more easily by using the
Aforementioned :py:meth:`~diofant.solvers.diophantine.descent`.

>>> from diofant.solvers.diophantine import descent
>>> descent(3, 1)  # solves the equation w**2 = 3*Y**2 + Z**2
(1, 0, 1)

Here the solution tuple is in the order (w, Y, Z)

The extended Pythagorean equation,
`a_{1}x_{1}^2 + a_{2}x_{2}^2 + \ldots + a_{n}x_{n}^2 = a_{n+1}x_{n+1}^2` and the
general sum of squares equation, `x_{1}^2 + x_{2}^2 + \ldots + x_{n}^2 = k` can
also be solved using the Diophantine module.

>>> from diofant.abc import e, f
>>> diophantine(9*a**2 + 16*b**2 + c**2 + 49*d**2 + 4*e**2 - 25*f**2)
{(70*t1**2 + 70*t2**2 + 70*t3**2 + 70*t4**2 - 70*t5**2, 105*t1*t5, 420*t2*t5, 60*t3*t5, 210*t4*t5, 42*t1**2 + 42*t2**2 + 42*t3**2 + 42*t4**2 + 42*t5**2)}

function :py:meth:`~diofant.solvers.diophantine.diop_general_pythagorean` can
also be called directly to solve the same equation. Either you can call
:py:meth:`~diofant.solvers.diophantine.diop_general_pythagorean` or use the high
level API. For the general sum of squares, this is also true, but one advantage
of calling :py:meth:`~diofant.solvers.diophantine.diop_general_sum_of_squares` is that
you can control how many solutions are returned.

>>> from diofant.solvers.diophantine import diop_general_sum_of_squares
>>> eq = a**2 + b**2 + c**2 + d**2 - 18
>>> diophantine(eq)
{(0, 0, 3, 3), (0, 1, 1, 4), (1, 2, 2, 3)}
>>> diop_general_sum_of_squares(eq, 2)
{(0, 0, 3, 3), (1, 2, 2, 3)}

The :py:meth:`~diofant.solvers.diophantine.sum_of_squares` routine will
providean iterator that returns solutions and one may control whether
the solutions contain zeros or not (and the solutions not containing
zeros are returned first):

>>> from diofant.solvers.diophantine import sum_of_squares
>>> sos = sum_of_squares(18, 4, zeros=True)
>>> next(sos)
(1, 2, 2, 3)
>>> next(sos)
(0, 0, 3, 3)

Simple Eqyptian fractions can be found with the Diophantine module, too.
For example, here are the ways that one might represent 1/2 as a sum of two
unit fractions:

>>> diophantine(Eq(1/x + 1/y, Rational(1, 2)))
{(-2, 1), (1, -2), (3, 6), (4, 4), (6, 3)}

To get a more thorough understanding of the Diophantine module, please
refer to the following blog.

https://thilinaatsympy.wordpress.com/


References
----------

* Andreescu, Titu. Andrica, Dorin. Cucurezeanu, Ion. An Introduction to
  Diophantine Equations
* Diophantine Equation, Wolfram Mathworld, [online]. Available:
  https://mathworld.wolfram.com/DiophantineEquation.html
* Methods to solve Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0,[online],
  Available: https://web.archive.org/web/20181231080858/https://www.alpertron.com.ar/METHODS.HTM
* Solving the equation ax^2+ bxy + cy^2 + dx + ey + f= 0, [online],
  Available: https://web.archive.org/web/20180831180321/http://www.jpr2718.org/ax2p.pdf

User Functions
--------------

This function is imported into the global namespace
with ``from diofant import *``:

diophantine
^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diophantine

And this function is imported with ``from diofant.solvers.diophantine import *``:

classify_diop
^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.classify_diop

Internal Functions
------------------
These functions are intended for internal use in the Diophantine module.

diop_solve
^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_solve

diop_linear
^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_linear

base_solution_linear
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.base_solution_linear

diop_quadratic
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_quadratic

diop_DN
^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_DN

cornacchia
^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.cornacchia

diop_bf_DN
^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_bf_DN

transformation_to_DN
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.transformation_to_DN

find_DN
^^^^^^^
.. autofunction:: diofant.solvers.diophantine.find_DN

diop_ternary_quadratic
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_ternary_quadratic

descent
^^^^^^^
.. autofunction:: diofant.solvers.diophantine.descent

diop_general_pythagorean
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_general_pythagorean

diop_general_sum_of_squares
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_general_sum_of_squares

diop_general_sum_of_even_powers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_general_sum_of_even_powers

partition
^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.partition

sum_of_three_squares
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.sum_of_three_squares

sum_of_four_squares
^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.sum_of_four_squares

power_representation
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.power_representation

.. function:: sum_of_powers

    alias of :func:`~diofant.solvers.diophantine.power_representation`

sum_of_squares
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.sum_of_squares

merge_solution
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.merge_solution

divisible
^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.divisible

PQa
^^^
.. autofunction:: diofant.solvers.diophantine.PQa

equivalent
^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.equivalent

parametrize_ternary_quadratic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.parametrize_ternary_quadratic

diop_ternary_quadratic_normal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.diop_ternary_quadratic_normal

ldescent
^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.ldescent

gaussian_reduce
^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.gaussian_reduce

holzer
^^^^^^
.. autofunction:: diofant.solvers.diophantine.holzer

prime_as_sum_of_two_squares
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.prime_as_sum_of_two_squares

square_factor
^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.square_factor

sqf_normal
^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.sqf_normal

reconstruct
^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.reconstruct

transformation_to_normal
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.diophantine.transformation_to_normal
