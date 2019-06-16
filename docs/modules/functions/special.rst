.. _special-functions:

Special
=======

DiracDelta
----------
.. autoclass:: diofant.functions.special.delta_functions.DiracDelta
   :members:

Heaviside
---------
.. autoclass:: diofant.functions.special.delta_functions.Heaviside
   :members:

Gamma, Beta and related Functions
---------------------------------
.. module:: diofant.functions.special.gamma_functions

.. autoclass:: diofant.functions.special.gamma_functions.gamma
   :members:
.. autoclass:: diofant.functions.special.gamma_functions.loggamma
   :members:
.. autoclass:: diofant.functions.special.gamma_functions.polygamma
   :members:
.. autofunction:: diofant.functions.special.gamma_functions.digamma
.. autofunction:: diofant.functions.special.gamma_functions.trigamma
.. autoclass:: diofant.functions.special.gamma_functions.uppergamma
   :members:
.. autoclass:: diofant.functions.special.gamma_functions.lowergamma
   :members:
.. module:: diofant.functions.special.beta_functions
.. autoclass:: diofant.functions.special.beta_functions.beta
   :members:

Error Functions and Fresnel Integrals
-------------------------------------
.. module:: diofant.functions.special.error_functions

.. autoclass:: diofant.functions.special.error_functions.erf
.. autoclass:: diofant.functions.special.error_functions.erfc
.. autoclass:: diofant.functions.special.error_functions.erfi
.. autoclass:: diofant.functions.special.error_functions.erf2
.. autoclass:: diofant.functions.special.error_functions.erfinv
.. autoclass:: diofant.functions.special.error_functions.erfcinv
.. autoclass:: diofant.functions.special.error_functions.erf2inv

.. autoclass:: diofant.functions.special.error_functions.FresnelIntegral
   :members:

.. autoclass:: fresnels
.. autoclass:: fresnelc

Exponential, Logarithmic and Trigonometric Integrals
----------------------------------------------------

.. autoclass:: Ei
.. autoclass:: expint
.. autofunction:: E1
.. autoclass:: li
.. autoclass:: Li
.. autoclass:: Si
.. autoclass:: Ci
.. autoclass:: Shi
.. autoclass:: Chi

Bessel Type Functions
---------------------

.. module:: diofant.functions.special.bessel

.. autoclass:: diofant.functions.special.bessel.BesselBase
   :members:

.. autoclass:: diofant.functions.special.bessel.besselj
.. autoclass:: diofant.functions.special.bessel.bessely
.. _besseli:
.. autoclass:: diofant.functions.special.bessel.besseli
.. autoclass:: diofant.functions.special.bessel.besselk
.. autoclass:: diofant.functions.special.bessel.hankel1
.. autoclass:: diofant.functions.special.bessel.hankel2
.. autoclass:: diofant.functions.special.bessel.jn
.. autoclass:: diofant.functions.special.bessel.yn

.. autofunction:: diofant.functions.special.bessel.jn_zeros

Airy Functions
--------------

.. autoclass:: diofant.functions.special.bessel.AiryBase
   :members:

.. autoclass:: diofant.functions.special.bessel.airyai
.. autoclass:: diofant.functions.special.bessel.airybi
.. autoclass:: diofant.functions.special.bessel.airyaiprime
.. autoclass:: diofant.functions.special.bessel.airybiprime

B-Splines
---------

.. autofunction:: diofant.functions.special.bsplines.bspline_basis
.. autofunction:: diofant.functions.special.bsplines.bspline_basis_set

Riemann Zeta and Related Functions
----------------------------------
.. module:: diofant.functions.special.zeta_functions

.. autoclass:: zeta
.. autoclass:: dirichlet_eta
.. autoclass:: polylog
.. autoclass:: lerchphi

Hypergeometric Functions
------------------------
.. autoclass:: diofant.functions.special.hyper.hyper
   :members:

.. autoclass:: diofant.functions.special.hyper.meijerg
   :members:

Elliptic integrals
------------------
.. module:: diofant.functions.special.elliptic_integrals

.. autoclass:: elliptic_k
.. autoclass:: elliptic_f
.. autoclass:: elliptic_e
.. autoclass:: elliptic_pi

Orthogonal Polynomials
----------------------

.. automodule:: diofant.functions.special.polynomials

Jacobi Polynomials
++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.jacobi
   :members:

.. autofunction:: diofant.functions.special.polynomials.jacobi_normalized

Gegenbauer Polynomials
++++++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.gegenbauer
   :members:

Chebyshev Polynomials
+++++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.chebyshevt
   :members:

.. autoclass:: diofant.functions.special.polynomials.chebyshevu
   :members:

.. autoclass:: diofant.functions.special.polynomials.chebyshevt_root
   :members:

.. autoclass:: diofant.functions.special.polynomials.chebyshevu_root
   :members:

Legendre Polynomials
++++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.legendre
   :members:

.. autoclass:: diofant.functions.special.polynomials.assoc_legendre
   :members:

Hermite Polynomials
+++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.hermite
   :members:

Laguerre Polynomials
++++++++++++++++++++

.. autoclass:: diofant.functions.special.polynomials.laguerre
   :members:
.. autoclass:: diofant.functions.special.polynomials.assoc_laguerre
   :members:

Spherical Harmonics
-------------------

.. autoclass:: diofant.functions.special.spherical_harmonics.Ynm

.. autofunction:: diofant.functions.special.spherical_harmonics.Ynm_c

.. autoclass:: diofant.functions.special.spherical_harmonics.Znm

Tensor Functions
----------------

.. autofunction:: diofant.functions.special.tensor_functions.Eijk

.. autofunction:: diofant.functions.special.tensor_functions.eval_levicivita

.. autoclass:: diofant.functions.special.tensor_functions.LeviCivita
   :members:

.. autoclass:: diofant.functions.special.tensor_functions.KroneckerDelta
   :members:
