"""Low-level linear systems solver."""

from ..matrices import Matrix, zeros


class RawMatrix(Matrix):
    """Dummy class with overriden sympify() helper."""

    _sympify = staticmethod(lambda x: x)


def eqs_to_matrix(eqs, ring):
    """Transform from equations to matrix form."""
    m = zeros(len(eqs), len(ring.gens) + 1, cls=RawMatrix)

    for j, e_j in enumerate(eqs):
        for i, x_i in enumerate(ring.gens):
            m[j, i] = e_j[x_i]
        m[j, -1] = -e_j[1]

    return m


def solve_lin_sys(eqs, ring):
    """Solve a system of linear equations."""
    # transform from equations to matrix form
    matrix = eqs_to_matrix(eqs, ring)

    # solve by row-reduction
    echelon, pivots = matrix.rref(iszerofunc=lambda x: not x)

    if not pivots:
        return {}
    elif pivots[-1] == len(ring.gens):
        return
    elif len(pivots) == len(ring.gens):
        sols = dict(zip(ring.gens, echelon[:, -1]))
    else:
        sols = {}
        for i, p in enumerate(pivots):
            vect = RawMatrix([[-x] for x in ring.gens[p+1:]] + [[ring.one]])
            sols[ring.gens[p]] = (echelon[i, p + 1:]*vect)[0]

    return {k: ring(v) for k, v in sols.items()}
