from ..core import AtomicExpr, Integer, Pow
from ..matrices import ImmutableMatrix
from .basisdependent import (BasisDependent, BasisDependentAdd,
                             BasisDependentMul, BasisDependentZero)


class Dyadic(BasisDependent):
    """
    Super class for all Dyadic-classes.

    References
    ==========

    * https://en.wikipedia.org/wiki/Dyadic_tensor
    * Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    """

    _op_priority = 13.0

    @property
    def components(self):
        """
        Returns the components of this dyadic in the form of a
        Python dictionary mapping BaseDyadic instances to the
        corresponding measure numbers.

        """
        # The '_components' attribute is defined according to the
        # subclass of Dyadic the instance belongs to.
        return self._components

    def dot(self, other):
        """
        Returns the dot product(also called inner product) of this
        Dyadic, with another Dyadic or Vector.
        If 'other' is a Dyadic, this returns a Dyadic. Else, it returns
        a Vector (unless an error is encountered).

        Parameters
        ==========

        other : Dyadic/Vector
            The other Dyadic or Vector to take the inner product with

        Examples
        ========

        >>> from diofant.vector import CoordSysCartesian
        >>> N = CoordSysCartesian('N')
        >>> D1 = N.i.outer(N.j)
        >>> D2 = N.j.outer(N.j)
        >>> D1.dot(D2)
        (N.i|N.j)
        >>> D1.dot(N.j)
        N.i

        """

        from .vector import Vector

        if isinstance(other, BasisDependentZero):
            return Vector.zero
        elif isinstance(other, Vector):
            outvec = Vector.zero
            for k, v in self.components.items():
                vect_dot = k.args[1].dot(other)
                outvec += vect_dot * v * k.args[0]
            return outvec
        elif isinstance(other, Dyadic):
            outdyad = Dyadic.zero
            for k1, v1 in self.components.items():
                for k2, v2 in other.components.items():
                    vect_dot = k1.args[1].dot(k2.args[0])
                    outer_product = k1.args[0].outer(k2.args[1])
                    outdyad += vect_dot * v1 * v2 * outer_product
            return outdyad
        else:
            raise TypeError("Inner product is not defined for " +
                            str(type(other)) + " and Dyadics.")

    def __and__(self, other):
        return self.dot(other)
    __and__.__doc__ = dot.__doc__

    def cross(self, other):
        """
        Returns the cross product between this Dyadic, and a Vector, as a
        Vector instance.

        Parameters
        ==========

        other : Vector
            The Vector that we are crossing this Dyadic with

        Examples
        ========

        >>> from diofant.vector import CoordSysCartesian
        >>> N = CoordSysCartesian('N')
        >>> d = N.i.outer(N.i)
        >>> d.cross(N.j)
        (N.i|N.k)

        """

        from .vector import Vector

        if other == Vector.zero:
            return Dyadic.zero
        elif isinstance(other, Vector):
            outdyad = Dyadic.zero
            for k, v in self.components.items():
                cross_product = k.args[1].cross(other)
                outer = k.args[0].outer(cross_product)
                outdyad += v * outer
            return outdyad
        else:
            raise TypeError(str(type(other)) + " not supported for " +
                            "cross with dyadics")

    def __xor__(self, other):
        return self.cross(other)
    __xor__.__doc__ = cross.__doc__

    def to_matrix(self, system, second_system=None):
        """
        Returns the matrix form of the dyadic with respect to one or two
        coordinate systems.

        Parameters
        ==========

        system : CoordSysCartesian
            The coordinate system that the rows and columns of the matrix
            correspond to. If a second system is provided, this
            only corresponds to the rows of the matrix.
        second_system : CoordSysCartesian, optional, default=None
            The coordinate system that the columns of the matrix correspond
            to.

        Examples
        ========

        >>> from diofant.vector import CoordSysCartesian
        >>> N = CoordSysCartesian('N')
        >>> v = N.i + 2*N.j
        >>> d = v.outer(N.i)
        >>> d.to_matrix(N)
        Matrix([
        [1, 0, 0],
        [2, 0, 0],
        [0, 0, 0]])
        >>> q = Symbol('q')
        >>> P = N.orient_new_axis('P', q, N.k)
        >>> d.to_matrix(N, P)
        Matrix([
        [  cos(q),   -sin(q), 0],
        [2*cos(q), -2*sin(q), 0],
        [       0,         0, 0]])

        """

        if second_system is None:
            second_system = system

        return ImmutableMatrix([i.dot(self).dot(j) for i in system for j in
                                second_system]).reshape(3, 3)


class BaseDyadic(Dyadic, AtomicExpr):
    """Class to denote a base dyadic tensor component."""

    def __new__(cls, vector1, vector2):
        from .vector import Vector, BaseVector, VectorZero

        # Verify arguments
        if not isinstance(vector1, (BaseVector, VectorZero)) or \
                not isinstance(vector2, (BaseVector, VectorZero)):
            raise TypeError("BaseDyadic cannot be composed of non-base " +
                            "vectors")
        # Handle special case of zero vector
        elif vector1 == Vector.zero or vector2 == Vector.zero:
            return Dyadic.zero
        # Initialize instance
        obj = super().__new__(cls, vector1, vector2)
        obj._base_instance = obj
        obj._measure_number = 1
        obj._components = {obj: Integer(1)}
        obj._sys = vector1._sys
        obj._pretty_form = ('(' + vector1._pretty_form + '|' +
                            vector2._pretty_form + ')')
        obj._latex_form = ('(' + vector1._latex_form + "{|}" +
                           vector2._latex_form + ')')

        return obj

    def __str__(self, printer=None):
        return "(" + str(self.args[0]) + "|" + str(self.args[1]) + ")"

    _diofantstr = __str__
    _diofantrepr = _diofantstr


class DyadicMul(BasisDependentMul, Dyadic):
    """Products of scalars and BaseDyadics."""

    def __new__(cls, *args, **options):
        obj = BasisDependentMul.__new__(cls, *args, **options)
        return obj

    @property
    def base_dyadic(self):
        """The BaseDyadic involved in the product."""
        return self._base_instance

    @property
    def measure_number(self):
        """ The scalar expression involved in the definition of
        this DyadicMul.

        """
        return self._measure_number


class DyadicAdd(BasisDependentAdd, Dyadic):
    """Class to hold dyadic sums."""

    def __new__(cls, *args, **options):
        obj = BasisDependentAdd.__new__(cls, *args, **options)
        return obj

    def __str__(self, printer=None):
        ret_str = ''
        items = list(self.components.items())
        items.sort(key=lambda x: x[0].__str__())
        for k, v in items:
            temp_dyad = k * v
            ret_str += temp_dyad.__str__(printer) + " + "
        return ret_str[:-3]

    __repr__ = __str__
    _diofantstr = __str__


class DyadicZero(BasisDependentZero, Dyadic):
    """Class to denote a zero dyadic."""

    _op_priority = 13.1
    _pretty_form = '(0|0)'
    _latex_form = r'(\mathbf{\hat{0}}|\mathbf{\hat{0}})'

    def __new__(cls):
        obj = BasisDependentZero.__new__(cls)
        return obj


def _dyad_div(one, other):
    """Helper for division involving dyadics."""
    if isinstance(other, Dyadic):
        raise TypeError("Cannot divide two dyadics")
    else:
        return DyadicMul(one, Pow(other, -1))


Dyadic._expr_type = Dyadic
Dyadic._mul_func = DyadicMul
Dyadic._add_func = DyadicAdd
Dyadic._zero_func = DyadicZero
Dyadic._base_func = BaseDyadic
Dyadic._div_helper = _dyad_div
Dyadic.zero = DyadicZero()
