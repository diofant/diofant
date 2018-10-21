============
Diofant 0.10
============

Not Released Yet

New features
============

* New representation for elements of :class:`~diofant.domains.AlgebraicField`, see :pull:`619` and :pull:`631`.
* Ground domains for :class:`~diofant.domains.AlgebraicField` can be also instances of :class:`~diofant.domains.AlgebraicField`, see :pull:`653`.
* New subclasses of :class:`~diofant.domains.AlgebraicField`: :class:`~diofant.domains.RealAlgebraicField` and :class:`~diofant.domains.ComplexAlgebraicField`, see :pull:`669` and :pull:`630`.

Major changes
=============

* Stable enumeration of polynomial roots in :class:`~diofant.polys.rootoftools.RootOf`, see :pull:`633` and :pull:`658`.
* Support root isolation for polynomials with algebraic coefficients, see :pull:`673` and :pull:`630`.

Compatibility breaks
====================

* Removed ``DMF`` class, see :pull:`620`.
* Removed ``K[x, y, ...]`` sugar, use :meth:`~diofant.domains.domain.Domain.poly_ring` to create polynomial rings, see :pull:`622`.
* Allow only prime orders for :class:`~diofant.domains.FiniteField`, see :pull:`622`.
* Removed ``FracField`` class, see :pull:`622`.
* ``get_field()`` method for domains, derived from :class:`~diofant.domains.ring.Ring`, now is a property, e.g. :attr:`~diofant.domains.field.Field.field`, see :pull:`622`.
* Removed ``PolyRing`` class, see :pull:`621`.
* ``get_ring()`` method for domains, derived from :class:`~diofant.domains.ring.Ring`, now is a property, e.g. :attr:`~diofant.domains.ring.Ring.ring`, see :pull:`621`.
* Removed ``compose`` option for :func:`~diofant.polys.numberfields.minimal_polynomial`, use ``method`` instead, see :pull:`624`.
* :func:`~diofant.polys.numberfields.field_isomorphism` take fields as arguments, see :pull:`627`.
* Functions :func:`~diofant.polys.numberfields.minimal_polynomial` and :func:`~diofant.polys.numberfields.primitive_element` return :class:`~diofant.polys.polytools.PurePoly` instances, see :pull:`628`.
* Removed ``ANP`` class, see :pull:`619`.
* Removed ``to_number_field()``, use :meth:`~diofant.domains.domain.Domain.convert` instead, see :pull:`619`.
* Removed ``RealNumber`` alias, see :pull:`635`.
* Removed ``of_type()`` method of :class:`~diofant.domains.domain.Domain`, see :pull:`636`.
* Method ``characteristic()`` now is a property of :class:`~diofant.domains.characteristiczero.CharacteristicZero` and :class:`~diofant.domains.FiniteField`, see :pull:`636`.
* Removed ``abs()``, ``is_one()`` and ``unify_with_symbols()`` methods and ``has_CharacteristicZero`` attribute of :class:`~diofant.domains.domain.Domain`, see :pull:`637`.
* Removed ``is_unit()``, ``numer()`` and ``denom()`` methods of :class:`~diofant.domains.ring.Ring`, see :pull:`637`.
* ``from_<Foo>()`` methods of :class:`~diofant.domains.domain.Domain` now are private, see :pull:`637`.
* Method :meth:`~diofant.domains.domain.Domain.from_expr` was renamed from ``from_diofant()``, see :pull:`637`.
* Method :meth:`~diofant.domains.domain.Domain.to_expr` was renamed from ``to_diofant()``, see :pull:`637`.
* Removed ``AlgebraicNumber`` class, see :pull:`631`.
* Removed ``polys.distributedmodules`` module, see :pull:`648`.
* Removed ``p`` and ``q`` properties of :class:`~diofant.core.numbers.Rational`, see :pull:`654`.
* Removed ``@public`` decorator, see :pull:`666`.
* Removed ``dummy_eq()`` method from :class:`~diofant.core.basic.Basic`, see :pull:`666`.
* :class:`~diofant.core.function.Subs` now support only ``Subs(expr, (var1, val1), (var2, val2), ...)`` syntax, see :pull:`667`.
* :class:`~diofant.polys.rootoftools.RootOf` don't canonicalize anymore polynomials to have integer coefficients, use :func:`~diofant.core.function.expand_func` instead, see :pull:`679`.
* Removed `Theano <https://github.com/Theano/Theano/>`_ support, see :pull:`681`.
* Removed ``minpoly`` alias for :func:`~diofant.polys.numberfields.minimal_polynomial`, see :pull:`684`.
* Removed ``abs()`` and ``neg()`` methods of :class:`~diofant.polys.polytools.Poly`, see :pull:`688`.
* Method :meth:`~diofant.polys.polytools.GroebnerBasis.set_order` was renamed from ``fglm()``, see :pull:`688`.
* Removed ``row()``, ``col()``, ``row_del()`` and ``col_del()`` methods of :class:`~diofant.matrices.Matrix`, see :pull:`688`.
* Removed ``add()`` and ``mul()`` methods for :class:`~diofant.polys.rings.PolynomialRing`, see :pull:`697`.
* Removed ``itercoeffs()``, ``itermonoms()``, ``iterterms()``, ``listcoeffs()``, ``listmonoms()``, ``listterms()``, ``const()``, ``imul_num()`` and ``square()`` methods of :class:`~diofant.polys.rings.PolyElement`, see :pull:`697`.

Minor changes
=============

* Be sure that :func:`~diofant.polys.numberfields.minimal_polynomial` returns an irreducible polynomial over specified domain, see :pull:`622`.
* Support algebraic function fields in :func:`~diofant.polys.numberfields.minpoly_groebner`, see :pull:`623`.
* Added argument ``method`` for :func:`~diofant.polys.numberfields.minimal_polynomial` and ``MINPOLY_METHOD`` configuration option to select default algorithm, see :pull:`624`.
* Support derivatives of :class:`~diofant.polys.rootoftools.RootOf` instances, see :pull:`624`.
* :func:`~diofant.polys.numberfields.primitive_element` now return an algebraic integer and support algebraic fields, see :pull:`643`, :pull:`655` and :pull:`659`.
* Support :class:`~diofant.functions.elementary.complexes.conjugate`, :class:`~diofant.functions.elementary.complexes.Abs`, :class:`~diofant.functions.elementary.complexes.re` and :class:`~diofant.functions.elementary.complexes.im` in :func:`~diofant.polys.numberfields.minimal_polynomial`, see :pull:`661` and :pull:`668`.
* :meth:`~diofant.polys.rootoftools.RootOf.refine` method to refine interval for the root, see :pull:`670`.
* Support detection of imaginary roots in :class:`~diofant.polys.rootoftools.RootOf`, see :pull:`625`.
* Mutable matrices support indexed deletion with :meth:`~object.__delitem__`, see :pull:`688`.
* Integer powers of :class:`~diofant.polys.rootoftools.RootOf` instances are automatically reduced, according to their minimal polynomial, see :pull:`691`.

Developer changes
=================

* Removed cachetools dependence, see :pull:`647`.
* Depend on `pylint <https://pylint.readthedocs.io/en/latest/>`_, see :pull:`668`.

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/3?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`14384` An unspecified power of x is reported to be O(log(x)**6)
* :sympyissue:`14393` Incorrect limit
* :sympyissue:`14414` Should QQ[x, y, ...] syntax be removed?
* :sympyissue:`13886` Raise an exception for non-prime p in FiniteFIeld(p)
* :sympyissue:`14220` Should be there both PolyRing and PolynomialRing?
* :sympyissue:`7724` roots should find the roots of x**4*I + x**2 + I
* :sympyissue:`5850` minpoly() should use PurePoly
* :sympyissue:`14494` make better decisions for minpoly based on domain
* :sympyissue:`14389` AlgebraicNumber should be a domain element?
* :sympyissue:`14291` poly(((x - 1)**2 + 1)*((x - 1)**2 + 2)*(x - 1)).all_roots() hangs
* :sympyissue:`14590` limit((n**3*((n + 1)/n)**n)/((n + 1)*(n + 2)*(n + 3)), n, oo) is incorrect
* :sympyissue:`14645` Bug when solving multivariate polynomial systems with identical equations
* :sympyissue:`14294` to_number_field should be idempotent for single extension
* :sympyissue:`14721` solve can't find solution
* :sympyissue:`14293` Sorting of polynomial roots
* :sympyissue:`14380` AlgebraicField.numer() could return an algebraic integer
* :sympyissue:`14442` Should AlgebraicField be a Composite domain?
* :sympyissue:`14759` dup_isolate_real_roots_list() docstring is wrong
* :sympyissue:`14738` dup_count_complex_roots() can't handle degenerate cases
* :sympyissue:`14782` integrate(sqrt(-x**2 + 1)*(-x**2 + x), [x, -1, 1]) is incorrect
* :sympyissue:`14791` No solution is returned for solve(exp(log(5)*x) - exp(log(2)*x), x)
* :sympyissue:`14793` Limit involving log(factorial(x)) incorrect
* :sympyissue:`14811` Exception during evaluation of limit (only locally, not in the live version)
* :sympyissue:`14822` RisingFactorial cannot do numerical (floating point) evaluations
* :sympyissue:`14820` octave/matlab codegen wrong for two argument zeta
* :sympyissue:`14831` minpoly(-3*sqrt(12*sqrt(2) + 17) + 12*sqrt(2) + 17 -2*sqrt(2)*sqrt(12*sqrt(2) + 17), x) fails
* :sympyissue:`14476` QQ.algebraic_field(Rational) should be just QQ
* :sympyissue:`14885` Sympy series gives TypeError on x^(-3/2) * exp(x) at x = 0
* :sympyissue:`15055` Incorrect limit of n**3*((-n - 1)*sin(1/n) + (n + 2)*sin(1/(n + 1)))/(-n + 1)
* :sympyissue:`15056` dsolve: get_numbered_constants should consider Functions
* :sympyissue:`6938` Undefined Functions should not use the evalf name lookup scheme
* :sympyissue:`8945` integrate(sin(x)**3/x, (x, 0, 1)) can't do it
* :sympyissue:`15146` Incorrect limit (n/2) * (-2*n**3 - 2*(n**3 - 1) * n**2 * digamma(n**3 + 1) + 2*(n**3 - 1) * n**2 * digamma(n**3 +n + 1) + n + 3)
* :sympyissue:`5934` PolynomialError with minpoly()
* :sympyissue:`8210` Zero degree polynomial copy() error
* :sympyissue:`11775` TypeError: unorderable types: PolyElement() < mpz() from factor_list
* :sympyissue:`7047` Python and gmpy ground type specific stuff from "from sympy import \*"
* :sympyissue:`15323` limit of the derivative of (1-1/x)^x as x --> 1+ gives wrong answer
* :sympyissue:`15344` mathematica_code gives wrong output with Max
* :sympyissue:`12602` count_roots is extremely slow with Python ground types
* :sympyissue:`5595` Should mpmath use the polys ground types?
* :sympyissue:`5602` Poly should use free_symbols to check for variable dependence
* :sympyissue:`5555` Explain coefficient domain handling in groebner()'s docstring
* :sympyissue:`15407` BUG: dsolve fails for linear first order ODE with three equations
