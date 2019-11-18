.. _ode-docs:

ODE
===

.. module::diofant.solvers.ode

User Functions
--------------
These are functions that are imported into the global namespace with ``from
diofant import *``.  These functions (unlike `Hint Functions`_, below) are
intended for use by ordinary users of Diofant.

dsolve
^^^^^^
.. autofunction:: diofant.solvers.ode.dsolve

classify_ode
^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.classify_ode

checkodesol
^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.checkodesol

homogeneous_order
^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.homogeneous_order

infinitesimals
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.infinitesimals

checkinfsol
^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.checkinfsol

Hint Functions
--------------
These functions are intended for internal use by
:py:meth:`~diofant.solvers.ode.dsolve` and others.  Unlike `User Functions`_,
above, these are not intended for every-day use by ordinary Diofant users.
Instead, functions such as :py:meth:`~diofant.solvers.ode.dsolve` should be used.
Nonetheless, these functions contain useful information in their docstrings on
the various ODE solving methods. For this reason, they are documented here.

allhints
^^^^^^^^
.. autodata:: diofant.solvers.ode.allhints

odesimp
^^^^^^^
.. autofunction:: diofant.solvers.ode.odesimp

constant_renumber
^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.constant_renumber

constantsimp
^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.constantsimp

sol_simplicity
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_sol_simplicity

1st_exact
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_exact

1st_homogeneous_coeff_best
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_homogeneous_coeff_best

1st_homogeneous_coeff_subs_dep_div_indep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_homogeneous_coeff_subs_dep_div_indep

1st_homogeneous_coeff_subs_indep_div_dep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_homogeneous_coeff_subs_indep_div_dep

1st_linear
^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_linear

Bernoulli
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_Bernoulli

Liouville
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_Liouville

Riccati_special_minus2
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_Riccati_special_minus2

nth_linear_constant_coeff_homogeneous
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_nth_linear_constant_coeff_homogeneous

nth_linear_constant_coeff_undetermined_coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_nth_linear_constant_coeff_undetermined_coefficients

nth_linear_constant_coeff_variation_of_parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_nth_linear_constant_coeff_variation_of_parameters

separable
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_separable

almost_linear
^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_almost_linear

linear_coefficients
^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_linear_coefficients

separable_reduced
^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_separable_reduced

lie_group
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_lie_group

1st_power_series
^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_1st_power_series

2nd_power_series_ordinary
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_2nd_power_series_ordinary

2nd_power_series_regular
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.ode_2nd_power_series_regular

Lie heuristics
--------------
These functions are intended for internal use of the Lie Group Solver.
Nonetheless, they contain useful information in their docstrings on the algorithms
implemented for the various heuristics.

abaco1_simple
^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_abaco1_simple

abaco1_product
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_abaco1_product

bivariate
^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_bivariate

chi
^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_chi

abaco2_similar
^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_abaco2_similar

function_sum
^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_function_sum

abaco2_unique_unknown
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_abaco2_unique_unknown

linear
^^^^^^
.. autofunction:: diofant.solvers.ode.lie_heuristic_linear

System of ODEs
--------------
These functions are intended for internal use by
:py:meth:`~diofant.solvers.ode.dsolve` for system of differential equations.

system_of_odes_linear_2eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order1_type3

system_of_odes_linear_2eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order1_type4

system_of_odes_linear_2eq_order1_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order1_type5

system_of_odes_linear_2eq_order1_type6
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order1_type6

system_of_odes_linear_2eq_order1_type7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order1_type7

system_of_odes_linear_2eq_order2_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type1

system_of_odes_linear_2eq_order2_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type2

system_of_odes_linear_2eq_order2_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type3

system_of_odes_linear_2eq_order2_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type5

system_of_odes_linear_2eq_order2_type6
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type6

system_of_odes_linear_2eq_order2_type7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type7

system_of_odes_linear_2eq_order2_type8
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type8

system_of_odes_linear_2eq_order2_type9
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type9

system_of_odes_linear_2eq_order2_type11
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_2eq_order2_type11

system_of_odes_linear_3eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._linear_3eq_order1_type4

system_of_odes_linear_neq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode.sysode_linear_neq_order1

system_of_odes_nonlinear_2eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_2eq_order1_type1

system_of_odes_nonlinear_2eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_2eq_order1_type2

system_of_odes_nonlinear_2eq_order1_type3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_2eq_order1_type3

system_of_odes_nonlinear_2eq_order1_type4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_2eq_order1_type4

system_of_odes_nonlinear_2eq_order1_type5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_2eq_order1_type5

system_of_odes_nonlinear_3eq_order1_type1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_3eq_order1_type1

system_of_odes_nonlinear_3eq_order1_type2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: diofant.solvers.ode._nonlinear_3eq_order1_type2

Information on the ode module
-----------------------------

.. automodule:: diofant.solvers.ode

.. autofunction:: diofant.solvers.ode._undetermined_coefficients_match

.. autofunction:: diofant.solvers.ode._handle_Integral
