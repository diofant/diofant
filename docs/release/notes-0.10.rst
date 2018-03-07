============
Diofant 0.10
============

Not Released Yet

New features
============

Major changes
=============

Compatibility breaks
====================

* Removed ``DMF`` class, see :pull:`620`.
* Removed ``K[x, y, ...]`` sugar, use :meth:`~diofant.domains.domain.Domain.poly_ring` to create polynomial rings, see :pull:`622`.
* Allow only prime orders for :class:`~diofant.domains.FiniteField`, see :pull:`622`.
* Removed ``FracField`` class, see :pull:`622`.
* ``get_field()`` method for domains, derived from :class:`~diofant.domains.ring.Ring`, now is a property, e.g. :attr:`~diofant.domains.field.Field.field`, see :pull:`622`.
* Removed ``PolyRing`` class, see :pull:`621`.
* ``get_ring()`` method for domains, derived from :class:`~diofant.domains.ring.Ring`, now is a property, e.g. :attr:`~diofant.domains.ring.Ring.ring`, see :pull:`621`.

Minor changes
=============

* Be sure that :func:`~diofant.polys.numberfields.minimal_polynomial` returns an irreducible polynomial over specified domain, see :pull:`622`.

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
