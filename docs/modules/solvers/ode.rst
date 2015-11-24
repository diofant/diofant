.. _ode-docs:

ODE
===

.. module::sympy.solvers.ode

User Functions
--------------
These are functions that are imported into the global namespace with ``from
sympy import *``.  These functions (unlike `Hint Functions`_, below) are
intended for use by ordinary users of SymPy.

dsolve
^^^^^^
.. autofunction:: sympy.solvers.ode.dsolve

classify_ode
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.classify_ode

checkodesol
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.checkodesol

homogeneous_order
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.homogeneous_order

infinitesimals
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.infinitesimals

checkinfsol
^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.checkinfsol

Hint Functions
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` and others.  Unlike `User Functions`_,
above, these are not intended for every-day use by ordinary SymPy users.
Instead, functions such as :py:meth:`~sympy.solvers.ode.dsolve` should be used.
Nonetheless, these functions contain useful information in their docstrings on
the various ODE solving methods. For this reason, they are documented here.

allhints
^^^^^^^^
.. autodata:: sympy.solvers.ode.allhints

odesimp
^^^^^^^
.. autofunction:: sympy.solvers.ode.odesimp

constant_renumber
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.constant_renumber

constantsimp
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.constantsimp

sol_simplicity
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_sol_simplicity

1st_exact
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_exact

1st_homogeneous_coeff_best
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_best

1st_homogeneous_coeff_subs_dep_div_indep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_subs_dep_div_indep

1st_homogeneous_coeff_subs_indep_div_dep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_homogeneous_coeff_subs_indep_div_dep

1st_linear
^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_linear

Bernoulli
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Bernoulli

Liouville
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Liouville

Riccati_special_minus2
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_Riccati_special_minus2

nth_linear_constant_coeff_homogeneous
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_homogeneous

nth_linear_constant_coeff_undetermined_coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_undetermined_coefficients

nth_linear_constant_coeff_variation_of_parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_nth_linear_constant_coeff_variation_of_parameters

separable
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_separable

almost_linear
^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_almost_linear

linear_coefficients
^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_linear_coefficients

separable_reduced
^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_separable_reduced

lie_group
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_lie_group

1st_power_series
^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_1st_power_series

2nd_power_series_ordinary
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_2nd_power_series_ordinary

2nd_power_series_regular
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.ode_2nd_power_series_regular

Lie heuristics
--------------
These functions are intended for internal use of the Lie Group Solver.
Nonetheless, they contain useful information in their docstrings on the algorithms
implemented for the various heuristics.

abaco1_simple
^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_abaco1_simple

abaco1_product
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_abaco1_product

bivariate
^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_bivariate

chi
^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_chi

abaco2_similar
^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_abaco2_similar

function_sum
^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_function_sum

abaco2_unique_unknown
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_abaco2_unique_unknown

abaco2_unique_general
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_abaco2_unique_general

linear
^^^^^^
.. autofunction:: sympy.solvers.ode.lie_heuristic_linear

System of ODEs
--------------
These functions are intended for internal use by
:py:meth:`~sympy.solvers.ode.dsolve` for system of differential equations.

system_of_odes_linear_2eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type1

system_of_odes_linear_2eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type2

system_of_odes_linear_2eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type3

system_of_odes_linear_2eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type4

system_of_odes_linear_2eq_order1_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type5

system_of_odes_linear_2eq_order1_type6
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type6

system_of_odes_linear_2eq_order1_type7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order1_type7

system_of_odes_linear_2eq_order2_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type1

system_of_odes_linear_2eq_order2_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type2

system_of_odes_linear_2eq_order2_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type3

system_of_odes_linear_2eq_order2_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type4

system_of_odes_linear_2eq_order2_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type5

system_of_odes_linear_2eq_order2_type6
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type6

system_of_odes_linear_2eq_order2_type7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type7

system_of_odes_linear_2eq_order2_type8
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type8

system_of_odes_linear_2eq_order2_type9
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type9

system_of_odes_linear_2eq_order2_type10
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type10

system_of_odes_linear_2eq_order2_type11
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_2eq_order2_type11

system_of_odes_linear_3eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_3eq_order1_type1

system_of_odes_linear_3eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_3eq_order1_type2

system_of_odes_linear_3eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_3eq_order1_type3

system_of_odes_linear_3eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_3eq_order1_type4

system_of_odes_linear_neq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._linear_neq_order1_type1

system_of_odes_nonlinear_2eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_2eq_order1_type1

system_of_odes_nonlinear_2eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_2eq_order1_type2

system_of_odes_nonlinear_2eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_2eq_order1_type3

system_of_odes_nonlinear_2eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_2eq_order1_type4

system_of_odes_nonlinear_2eq_order1_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_2eq_order1_type5

system_of_odes_nonlinear_3eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_3eq_order1_type1

system_of_odes_nonlinear_3eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_3eq_order1_type2

system_of_odes_nonlinear_3eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_3eq_order1_type3

system_of_odes_nonlinear_3eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_3eq_order1_type4

system_of_odes_nonlinear_3eq_order1_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: sympy.solvers.ode._nonlinear_3eq_order1_type5

Information on the ode module
-----------------------------

.. automodule:: sympy.solvers.ode

.. autofunction:: sympy.solvers.ode._undetermined_coefficients_match

.. autofunction:: sympy.solvers.ode._handle_Integral
