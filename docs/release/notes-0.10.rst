============
Diofant 0.10
============

Not Released Yet

New features
============

Major changes
=============

* New representation for elements of :class:`~diofant.domains.AlgebraicField`, see :pull:`619`.

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
* Removed ``alias`` option for :class:`~diofant.core.numbers.AlgebraicNumber`, see :pull:`626`.
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

Minor changes
=============

* Be sure that :func:`~diofant.polys.numberfields.minimal_polynomial` returns an irreducible polynomial over specified domain, see :pull:`622`.
* Support algebraic function fields in :func:`~diofant.polys.numberfields.minpoly_groebner`, see :pull:`623`.
* Added argument ``method`` for :func:`~diofant.polys.numberfields.minimal_polynomial` and ``MINPOLY_METHOD`` configuration option to select default algorithm, see :pull:`624`.
* Support derivatives of :class:`~diofant.polys.rootoftools.RootOf` instances, see :pull:`624`.

Developer changes
=================

Issues closed
=============

See the `release milestone <https://github.com/diofant/diofant/milestone/3?closed=1>`_
for complete list of issues and pull requests involved in this release.

These Sympy issues also were addressed:

* :sympyissue:`14384` An unspecified power of x is reported to be `O(log(x)**6)`
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
