.. _polys-reference:

=========================================
Polynomials Manipulation Module Reference
=========================================

Basic polynomial manipulation functions
=======================================

.. currentmodule:: diofant.polys.polytools

.. autofunction:: poly
.. autofunction:: poly_from_expr
.. autofunction:: parallel_poly_from_expr
.. autofunction:: degree
.. autofunction:: degree_list
.. autofunction:: LC
.. autofunction:: LM
.. autofunction:: LT
.. autofunction:: pdiv
.. autofunction:: prem
.. autofunction:: pquo
.. autofunction:: pexquo
.. autofunction:: div
.. autofunction:: rem
.. autofunction:: quo
.. autofunction:: exquo
.. autofunction:: half_gcdex
.. autofunction:: gcdex
.. autofunction:: invert
.. autofunction:: subresultants
.. autofunction:: resultant
.. autofunction:: discriminant
.. autofunction:: terms_gcd
.. autofunction:: cofactors
.. autofunction:: gcd
.. autofunction:: gcd_list
.. autofunction:: lcm
.. autofunction:: lcm_list
.. autofunction:: trunc
.. autofunction:: monic
.. autofunction:: content
.. autofunction:: primitive
.. autofunction:: compose
.. autofunction:: decompose
.. autofunction:: sturm
.. autofunction:: gff_list
.. autofunction:: gff
.. autofunction:: sqf_norm
.. autofunction:: sqf_part
.. autofunction:: sqf_list
.. autofunction:: sqf
.. autofunction:: factor_list
.. autofunction:: factor
.. autofunction:: intervals
.. autofunction:: refine_root
.. autofunction:: count_roots
.. autofunction:: real_roots
.. autofunction:: nroots
.. autofunction:: ground_roots
.. autofunction:: nth_power_roots_poly
.. autofunction:: cancel
.. autofunction:: reduced
.. autofunction:: groebner

.. autoclass:: Poly
   :members:

.. autoclass:: PurePoly
   :members:

.. autoclass:: GroebnerBasis
   :members:

Extra polynomial manipulation functions
=======================================

.. currentmodule:: diofant.polys.polyfuncs

.. autofunction:: symmetrize
.. autofunction:: horner
.. autofunction:: interpolate
.. autofunction:: viete

Domain constructors
===================

.. currentmodule:: diofant.polys.constructor

.. autofunction:: construct_domain

Algebraic number fields
=======================

.. currentmodule:: diofant.polys.numberfields

.. autofunction:: minimal_polynomial

.. function:: minpoly

    alias of :func:`minimal_polynomial`

.. autofunction:: minpoly_groebner

.. autofunction:: primitive_element
.. autofunction:: field_isomorphism

Monomials encoded as tuples
===========================

.. currentmodule:: diofant.polys.monomials

.. autoclass:: Monomial
.. autofunction:: itermonomials
.. autofunction:: monomial_count

Orderings of monomials
======================

.. currentmodule:: diofant.polys.orderings

.. autoclass:: LexOrder
.. autoclass:: GradedLexOrder
.. autoclass:: ReversedGradedLexOrder

Formal manipulation of roots of polynomials
===========================================

.. currentmodule:: diofant.polys.rootoftools

.. autoclass:: RootOf
    :members:

.. autoclass:: RootSum
    :members:

Symbolic root-finding algorithms
================================

.. currentmodule:: diofant.polys.polyroots

.. autofunction:: roots

Special polynomials
===================

.. currentmodule:: diofant.polys.specialpolys

.. autofunction:: swinnerton_dyer_poly
.. autofunction:: interpolating_poly
.. autofunction:: cyclotomic_poly
.. autofunction:: symmetric_poly
.. autofunction:: random_poly

Orthogonal polynomials
======================

.. currentmodule:: diofant.polys.orthopolys

.. autofunction:: chebyshevt_poly
.. autofunction:: chebyshevu_poly
.. autofunction:: gegenbauer_poly
.. autofunction:: hermite_poly
.. autofunction:: jacobi_poly
.. autofunction:: legendre_poly
.. autofunction:: laguerre_poly
.. autofunction:: spherical_bessel_fn

Manipulation of rational functions
==================================

.. currentmodule:: diofant.polys.rationaltools

.. autofunction:: together

Partial fraction decomposition
==============================

.. currentmodule:: diofant.polys.partfrac

.. autofunction:: apart
.. autofunction:: apart_list
.. autofunction:: assemble_partfrac_list

Dispersion of Polynomials
=========================

.. currentmodule:: diofant.polys.dispersion

.. autofunction:: dispersionset
.. autofunction:: dispersion
