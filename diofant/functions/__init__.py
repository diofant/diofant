"""A functions module, includes all the standard functions.

Combinatorial - factorial, fibonacci, harmonic, bernoulli...
Elementary - hyperbolic, trigonometric, exponential, floor and ceiling, sqrt...
Special - gamma, zeta,spherical harmonics...
"""

from .combinatorial.factorials import (factorial, factorial2, rf, ff,  # noqa: F401
                                       binomial, RisingFactorial,
                                       FallingFactorial, subfactorial)
from .combinatorial.numbers import (fibonacci, lucas, harmonic, bernoulli,  # noqa: F401
                                    bell, euler, catalan, genocchi)
from .elementary.miscellaneous import (sqrt, root, Min, Max, Id,  # noqa: F401
                                       real_root, cbrt)
from .elementary.complexes import (re, im, sign, Abs, conjugate, arg,  # noqa: F401
                                   polar_lift, periodic_argument,
                                   unbranched_argument, principal_branch,
                                   transpose, adjoint, polarify, unpolarify)
from .elementary.trigonometric import (sin, cos, tan, sec, csc, cot, asin,  # noqa: F401
                                       acos, atan, asec, acsc, acot, atan2)
from .elementary.exponential import exp_polar, exp, log, LambertW  # noqa: F401
from .elementary.hyperbolic import (sinh, cosh, tanh, coth, sech, csch,  # noqa: F401
                                    asinh, acosh, atanh, acoth)
from .elementary.integers import floor, ceiling  # noqa: F401
from .elementary.piecewise import Piecewise, piecewise_fold  # noqa: F401
from .special.error_functions import (erf, erfc, erfi, erf2, erfinv,  # noqa: F401
                                      erfcinv, erf2inv, Ei, expint, E1, li,
                                      Li, Si, Ci, Shi, Chi, fresnels,
                                      fresnelc)
from .special.gamma_functions import (gamma, lowergamma, uppergamma,  # noqa: F401
                                      polygamma, loggamma, digamma,
                                      trigamma)
from .special.zeta_functions import dirichlet_eta, zeta, lerchphi, polylog  # noqa: F401
from .special.tensor_functions import Eijk, LeviCivita, KroneckerDelta  # noqa: F401
from .special.delta_functions import DiracDelta, Heaviside  # noqa: F401
from .special.bsplines import bspline_basis, bspline_basis_set  # noqa: F401
from .special.bessel import (besselj, bessely, besseli, besselk, hankel1,  # noqa: F401
                             hankel2, jn, yn, jn_zeros, airyai, airybi,
                             airyaiprime, airybiprime)
from .special.hyper import hyper, meijerg  # noqa: F401
from .special.polynomials import (legendre, assoc_legendre, hermite,  # noqa: F401
                                  chebyshevt, chebyshevu, chebyshevu_root,
                                  chebyshevt_root, laguerre, assoc_laguerre,
                                  gegenbauer, jacobi, jacobi_normalized)
from .special.spherical_harmonics import Ynm, Ynm_c, Znm  # noqa: F401
from .special.elliptic_integrals import (elliptic_k, elliptic_f, elliptic_e,  # noqa: F401
                                         elliptic_pi)
from .special.beta_functions import beta  # noqa: F401

ln = log
