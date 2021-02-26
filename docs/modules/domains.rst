=======
Domains
=======

.. module:: diofant.domains

Here we document the various implemented ground domains.  There are
three types: abstract domains, concrete domains, and "implementation
domains".  Abstract domains cannot be (usefully) instantiated at all,
and just collect together functionality shared by many other domains.
Concrete domains are those meant to be instantiated and used.  In some
cases, there are various possible ways to implement the data type the
domain provides.  For example, depending on what libraries are
available on the system, the integers are implemented either using the
python built-in integers, or using gmpy.  Note that various aliases
are created automatically depending on the libraries available.  As
such e.g. ``ZZ`` always refers to the most efficient implementation of
the integer ring available.

Abstract Domains
****************

.. autoclass:: diofant.domains.domain.Domain
   :members:

.. autoclass:: diofant.domains.field.Field
   :members:

.. autoclass:: diofant.domains.ring.CommutativeRing
   :members:

.. autoclass:: diofant.domains.simpledomain.SimpleDomain
   :members:

.. autoclass:: diofant.domains.compositedomain.CompositeDomain
   :members:

.. autoclass:: diofant.domains.characteristiczero.CharacteristicZero
   :members:

Concrete Domains
****************

.. autoclass:: IntegerModRing
   :members:

.. autoclass:: FiniteField
   :members:

.. autoclass:: IntegerRing
   :members:

.. autoclass:: RationalField
   :members:

.. autoclass:: AlgebraicField
   :members:

.. autoclass:: RealAlgebraicField
   :members:

.. autoclass:: ComplexAlgebraicField
   :members:

.. autoclass:: diofant.polys.rings.PolynomialRing
   :members:

.. autoclass:: diofant.polys.univar.UnivarPolynomialRing
   :members:

.. autoclass:: diofant.polys.fields.FractionField
   :members:

.. autoclass:: RealField
   :members:

.. autoclass:: ComplexField
   :members:

.. autoclass:: ExpressionDomain
   :members:

Implementation Domains
**********************

.. autoclass:: diofant.domains.finitefield.PythonFiniteField
.. autoclass:: diofant.domains.finitefield.GMPYFiniteField

.. autoclass:: diofant.domains.integerring.PythonIntegerRing
.. autoclass:: diofant.domains.integerring.GMPYIntegerRing

.. autoclass:: diofant.domains.rationalfield.PythonRationalField
.. autoclass:: diofant.domains.rationalfield.GMPYRationalField

Domain Elements
***************

.. autoclass:: diofant.domains.finitefield.ModularInteger
   :members:

.. autoclass:: diofant.domains.finitefield.GaloisFieldElement
   :members:
