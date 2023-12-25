============
Diofant 0.12
============

18 Jan 2021

New features
============

* Support modular exponentiation of :class:`~diofant.polys.rings.PolyElement`'s, see :pull:`1032`.
* :func:`~diofant.solvers.inequalities.reduce_inequalities` support solving linear inequalities with Fourier-Motzkin elimination algorithm, see :pull:`1063`.
* Added class ``FiniteRing`` for modular integers, see :pull:`876`.
* Implemented :meth:`~diofant.polys.fields.FracElement.compose` for functional composition in the fields of fractions, see :pull:`1100`.

Major changes
=============

* Module :mod:`~diofant.polys.sqfreetools` was ported to use sparse polynomial representation, see :pull:`1009`.
* Module :mod:`~diofant.polys.factortools` was ported to use sparse polynomial representation, see :pull:`1015`, :pull:`1018`, :pull:`1019`, :pull:`1020` and :pull:`1021`.
* Module :mod:`~diofant.polys.rootisolation` was ported to use sparse polynomial representation, finally the dense representation is used nowhere, see :pull:`1030`, :pull:`1031` and :pull:`1035`.
* :func:`~diofant.solvers.inequalities.reduce_inequalities` uses :class:`~diofant.sets.fancysets.ExtendedReals` subsets to solve inequalities, see :pull:`1067` and :pull:`1092`.
* Added new algorithm for factorization of multivariate polynomials over :class:`~diofant.domains.AlgebraicField`'s (uses Hensel lifting), see :pull:`876`.  Thanks to Katja Sophie Hotz.  Thanks to Kalevi Suominen for help with review.

Compatibility breaks
====================

* Removed ``vring()`` and ``vfield()`` functions, see :pull:`1016`.
* Drop support for ``from_list()`` initialization for multivariate polynomials, see :pull:`1035`.
* Drop ``to_dense()``, ``tail_degrees()``, ``almosteq()`` and ``degree_list()`` methods and ``is_monic``, ``is_primitive`` attributes of :class:`~diofant.polys.rings.PolyElement`, see :pull:`1035`, :pull:`1036` and :pull:`1051`.
* Drop ``is_monic``, ``is_primitive``, ``zero``, ``one`` and ``unit`` attributes and ``degree_list()`` method of :class:`~diofant.polys.polytools.Poly`, see :pull:`1036`, :pull:`1039` and :pull:`1051`.
* Drop ``sring()``, ``poly_from_expr()``, ``gcd_list()`` and ``lcm_list()`` functions, see :pull:`1037`, :pull:`1057` and :pull:`1086`.
* Functions and classes of the :mod:`~diofant.polys.polytools` module do not support anymore iterables as polynomial generator, see :pull:`1039`.
* Drop unused functions ``dispersion()``, ``dispersionset()`` and ``degree_list()``, see :pull:`1051` and :pull:`1053`.
* Drop rich comparison methods from the :class:`~diofant.polys.fields.FracElement`, see :pull:`1101`.
* :func:`~diofant.polys.polytools.Poly.from_list` support now ascending order of coefficients (i.e., the leading coefficient of univariate polynomial is coming last), see :pull:`1103`.
* Removed support for 3D geometry in the ``geometry`` module and ``Point.__getitem__()`` method, see :pull:`1105`.
* Drop ``coeff()``, ``coeffs()``, ``monoms()``, ``terms()`` and ``deflate()`` methods of :class:`~diofant.polys.rings.PolyElement`, use dictionary indexing, see :pull:`1108`.

Minor changes
=============

* Special case univariate polynomials with :class:`~diofant.polys.univar.UnivarPolynomialRing` and :class:`~diofant.polys.univar.UnivarPolyElement`, see :pull:`1024`.
* Implement :attr:`~diofant.domains.finitefield.ModularInteger.is_primitive`, see :pull:`1035`.
* Add :class:`~diofant.sets.fancysets.ExtendedReals` singleton, see :pull:`1067`.
* 100% test coverage for ``geometry`` module, see :pull:`1105`.  Overall test coverage is around 98%.

Developer changes
=================

* Depend on `flake8-sfs <https://github.com/peterjc/flake8-sfs>`_, see :pull:`983`.
* Depend on `mypy <http://mypy-lang.org/>`_, see :pull:`1046`.
* Drop dependency on strategies, see :pull:`1074`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/6?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`19630` ``rsolve`` gives None for linear homogeneous recurrence relation
* :sympyissue:`19076` modular exponentiation of poly
* :sympyissue:`19670` Poly(E**100000000) is slow to create
* :sympyissue:`19755` poly gives coercion error when integers and rationals are mixed
* :sympyissue:`19760` minimal_polynomial using Groebner basis can give wrong result
* :sympyissue:`19770` Limit involving cosine
* :sympyissue:`19766` Incorrect limit
* :sympyissue:`19774` evalf() doesn't evaluate terms in an exponential
* :sympyissue:`19988` Float loses precision after being pickled
* :sympyissue:`14874` Limit x --> oo for besselk
* :sympyissue:`19991` Wrong result from floor().evalf()
* :sympyissue:`10666` resultant misses the sign
* :sympyissue:`20163` Apart hangs with extension=[sqrt(3), I]
* :sympyissue:`9479` Cannot solve multivariate inequalities
* :sympyissue:`20365` Limit Bug
* :sympyissue:`20360` Incorrect definite integration of simple exponential involving pi
* :sympyissue:`20389` TypeError: Argument of Integer should be of numeric type, got -oo
* :sympyissue:`20391` Linear programming with simplex method
* :sympyissue:`19161` When applying simplify on a Poly it fails
* :sympyissue:`20397` bug in dividing polynomials by module
* :sympyissue:`19196` Slow f.factor_list
* :sympyissue:`20491` Inconsistencies in pretty printing in a notebook
* :sympyissue:`20490` LaTeX printing of negative constant PolyElement
* :sympyissue:`20484` Need more utility for polynomial substitution
* :sympyissue:`20485` Rational powers for non-monomial PolyElement
* :sympyissue:`20487` LaTeX printing errors for puiseux polynomial
* :sympyissue:`20610` Solve: GeneratorsNeeded with system involving constant equation
* :sympyissue:`20617` Complex exponentials are not recognized by domains
* :sympyissue:`20640` Multivariate polynomial division
* :sympyissue:`20704` Limit not terminating
