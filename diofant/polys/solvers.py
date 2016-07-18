"""Low-level linear systems solver. """

from diofant.matrices import Matrix, zeros


class RawMatrix(Matrix):
    _sympify = staticmethod(lambda x: x)


def eqs_to_matrix(eqs, ring):
    """Transform from equations to matrix form. """
    m = zeros(len(eqs), len(ring.gens) + 1, cls=RawMatrix)

    for j, e_j in enumerate(eqs):
        for i, x_i in enumerate(ring.gens):
            m[j, i] = e_j.coeff(x_i)
        m[j, -1] = -e_j.coeff(1)

    return m


def solve_lin_sys(eqs, ring):
    """Solve a system of linear equations. """

    assert ring.domain.has_Field

    # transform from equations to matrix form
    matrix = eqs_to_matrix(eqs, ring)

    # solve by row-reduction
    echelon, pivots = matrix.rref(iszerofunc=lambda x: not x)

    if not pivots:
        return {}
    elif pivots[-1] == len(ring.gens):
        return
    elif len(pivots) == len(ring.gens):
        sol = [ring.ground_new(s) for s in echelon[:, -1]]
        return dict(zip(ring.gens, sol))
    else:
        sols = {}
        for i, p in enumerate(pivots):
            vect = RawMatrix([[-x] for x in ring.gens[p+1:]] + [[ring.one]])
            sols[ring.gens[p]] = (echelon[i, p + 1:]*vect)[0]

        return sols
