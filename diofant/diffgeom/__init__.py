"""
Package for differential geometry.
"""

from .diffgeom import (  # noqa: F401
    BaseCovarDerivativeOp, BaseScalarField, BaseVectorField, Commutator,
    contravariant_order, CoordSystem, CovarDerivativeOp, covariant_order,
    Differential, intcurve_diffequ, intcurve_series, LieDerivative,
    Manifold, metric_to_Christoffel_1st, metric_to_Christoffel_2nd,
    metric_to_Ricci_components, metric_to_Riemann_components, Patch,
    Point, TensorProduct, twoform_to_matrix, vectors_in_basis,
    WedgeProduct,
)
