import pytest

from diofant import Matrix, eye, sqrt, symbols
from diofant.tensor.tensor import (TensorIndexType, TensorType, tensor_indices,
                                   tensorhead, tensorsymmetry)


__all__ = ()

numpy = pytest.importorskip('numpy')


def _get_valued_base_test_variables():
    minkowski = Matrix((
        (1, 0, 0, 0),
        (0, -1, 0, 0),
        (0, 0, -1, 0),
        (0, 0, 0, -1),
    ))
    Lorentz = TensorIndexType('Lorentz', dim=4)
    Lorentz.data = minkowski

    i0, i1, i2, i3, i4 = tensor_indices('i0:5', Lorentz)

    E, px, py, pz = symbols('E px py pz')
    A = tensorhead('A', [Lorentz], [[1]])
    A.data = [E, px, py, pz]
    B = tensorhead('B', [Lorentz], [[1]], 'Gcomm')
    B.data = range(4)
    AB = tensorhead('AB', [Lorentz] * 2, [[1]]*2)
    AB.data = minkowski

    ba_matrix = Matrix((
        (1, 2, 3, 4),
        (5, 6, 7, 8),
        (9, 0, -1, -2),
        (-3, -4, -5, -6),
    ))

    BA = tensorhead('BA', [Lorentz] * 2, [[1]]*2)
    BA.data = ba_matrix

    # Let's test the diagonal metric, with inverted Minkowski metric:
    LorentzD = TensorIndexType('LorentzD')
    LorentzD.data = [-1, 1, 1, 1]
    mu0, mu1, mu2 = tensor_indices('mu0:3', LorentzD)
    C = tensorhead('C', [LorentzD], [[1]])
    C.data = [E, px, py, pz]

    # ### non-diagonal metric ###
    ndm_matrix = (
        (1, 1, 0,),
        (1, 0, 1),
        (0, 1, 0,),
    )
    ndm = TensorIndexType('ndm')
    ndm.data = ndm_matrix
    n0, n1, n2 = tensor_indices('n0:3', ndm)
    NA = tensorhead('NA', [ndm], [[1]])
    NA.data = range(10, 13)
    NB = tensorhead('NB', [ndm]*2, [[1]]*2)
    NB.data = [[i+j for j in range(10, 13)] for i in range(10, 13)]
    NC = tensorhead('NC', [ndm]*3, [[1]]*3)
    NC.data = [[[i+j+k for k in range(4, 7)] for j in range(1, 4)] for i in range(2, 5)]

    return (A, B, AB, BA, C, Lorentz, E, px, py, pz, LorentzD, mu0, mu1, mu2, ndm, n0, n1,
            n2, NA, NB, NC, minkowski, ba_matrix, ndm_matrix, i0, i1, i2, i3, i4)


def test_valued_tensor_iter():
    (A, _, _, BA, _, _, E, px, py, pz, _, _, _, _, _, _, _,
     _, _, _, _, _, ba_matrix, _, _, i1, i2, *_) = _get_valued_base_test_variables()

    # iteration on VTensorHead
    assert list(A) == [E, px, py, pz]
    assert list(ba_matrix) == list(BA)

    # iteration on VTensMul
    assert list(A(i1)) == [E, px, py, pz]
    assert list(BA(i1, i2)) == list(ba_matrix)
    assert list(3 * BA(i1, i2)) == [3 * i for i in list(ba_matrix)]
    assert list(-5 * BA(i1, i2)) == [-5 * i for i in list(ba_matrix)]

    # iteration on VTensAdd
    # A(i1) + A(i1)
    assert list(A(i1) + A(i1)) == [2*E, 2*px, 2*py, 2*pz]
    assert BA(i1, i2) - BA(i1, i2) == 0
    assert list(BA(i1, i2) - 2 * BA(i1, i2)) == [-i for i in list(ba_matrix)]


def test_valued_tensor_covariant_contravariant_elements():
    (A, _, AB, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, *_) = _get_valued_base_test_variables()

    assert A(-i0)[0] == A(i0)[0]
    assert A(-i0)[1] == -A(i0)[1]

    assert AB(i0, i1)[1, 1] == -1
    assert AB(i0, -i1)[1, 1] == 1
    assert AB(-i0, -i1)[1, 1] == -1
    assert AB(-i0, i1)[1, 1] == 1


def test_valued_tensor_get_matrix():
    (A, _, AB, _, _, _, E, px, py, pz, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, *_) = _get_valued_base_test_variables()

    matab = AB(i0, i1).get_matrix()
    assert matab == Matrix([
        [1,  0,  0,  0],
        [0, -1,  0,  0],
        [0,  0, -1,  0],
        [0,  0,  0, -1],
    ])
    # when alternating contravariant/covariant with [1, -1, -1, -1] metric
    # it becomes the identity matrix:
    assert AB(i0, -i1).get_matrix() == eye(4)

    # covariant and contravariant forms:
    assert A(i0).get_matrix() == Matrix([E, px, py, pz])
    assert A(-i0).get_matrix() == Matrix([E, -px, -py, -pz])


def test_valued_tensor_contraction():
    (A, B, AB, _, C, _, E, px, py, pz, _, mu0, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, *_) = _get_valued_base_test_variables()

    assert (A(i0) * A(-i0)).data == E ** 2 - px ** 2 - py ** 2 - pz ** 2
    assert (A(i0) * A(-i0)).data == A ** 2
    assert (A(i0) * A(-i0)).data == A(i0) ** 2
    assert (A(i0) * B(-i0)).data == -px - 2 * py - 3 * pz

    for i in range(4):
        for j in range(4):
            assert (A(i0) * B(-i1))[i, j] == [E, px, py, pz][i] * [0, -1, -2, -3][j]

    # test contraction on the alternative Minkowski metric: [-1, 1, 1, 1]
    assert (C(mu0) * C(-mu0)).data == -E ** 2 + px ** 2 + py ** 2 + pz ** 2

    contrexp = A(i0) * AB(i1, -i0)
    assert A(i0).rank == 1
    assert AB(i1, -i0).rank == 2
    assert contrexp.rank == 1
    for i in range(4):
        assert contrexp[i] == [E, px, py, pz][i]


def test_valued_tensor_self_contraction():
    (_, _, AB, BA, _, _, _, _, _, _, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, *_) = _get_valued_base_test_variables()

    assert AB(i0, -i0).data == 4
    assert BA(i0, -i0).data == 2


def test_valued_tensor_pow():
    (_, _, _, _, C, _, E, px, py, pz, _, mu0, _, _, _, _, _,
     *_) = _get_valued_base_test_variables()

    assert C**2 == -E**2 + px**2 + py**2 + pz**2
    assert C**1 == sqrt(-E**2 + px**2 + py**2 + pz**2)
    assert C(mu0)**2 == C**2
    assert C(mu0)**1 == C**1


def test_valued_tensor_expressions():
    (A, B, AB, BA, _, _, E, px, py, pz, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, i2, i3, i4) = _get_valued_base_test_variables()

    x1, x2, x3 = symbols('x1:4')

    # test coefficient in contraction:
    rank2coeff = x1 * A(i3) * B(i2)
    assert rank2coeff[1, 1] == x1 * px
    assert rank2coeff[3, 3] == 3 * pz * x1
    coeff_expr = ((x1 * A(i4)) * (B(-i4) / x2)).data

    assert coeff_expr.expand() == -px*x1/x2 - 2*py*x1/x2 - 3*pz*x1/x2

    add_expr = A(i0) + B(i0)

    assert add_expr[0] == E
    assert add_expr[1] == px + 1
    assert add_expr[2] == py + 2
    assert add_expr[3] == pz + 3

    sub_expr = A(i0) - B(i0)

    assert sub_expr[0] == E
    assert sub_expr[1] == px - 1
    assert sub_expr[2] == py - 2
    assert sub_expr[3] == pz - 3

    assert (add_expr * B(-i0)).data == -px - 2*py - 3*pz - 14

    expr1 = x1*A(i0) + x2*B(i0)
    expr2 = expr1 * B(i1) * (-4)
    expr3 = expr2 + 3*x3*AB(i0, i1)
    expr4 = expr3 / 2
    assert expr4 * 2 == expr3
    expr5 = expr4 * BA(-i1, -i0)

    assert expr5.data.expand() == 28*E*x1 + 12*px*x1 + 20*py*x1 + 28*pz*x1 + 136*x2 + 3*x3


def test_valued_tensor_add_scalar():
    (A, _, _, _, C, _, E, px, py, pz, _, mu0, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, *_) = _get_valued_base_test_variables()

    # one scalar summand after the contracted tensor
    expr1 = A(i0)*A(-i0) - (E**2 - px**2 - py**2 - pz**2)
    assert expr1.data == 0

    # multiple scalar summands in front of the contracted tensor
    expr2 = E**2 - px**2 - py**2 - pz**2 - A(i0)*A(-i0)
    assert expr2.data == 0

    # multiple scalar summands after the contracted tensor
    expr3 = A(i0)*A(-i0) - E**2 + px**2 + py**2 + pz**2
    assert expr3.data == 0

    # multiple scalar summands and multiple tensors
    expr4 = C(mu0)*C(-mu0) + 2*E**2 - 2*px**2 - 2*py**2 - 2*pz**2 - A(i0)*A(-i0)
    assert expr4.data == 0


def test_noncommuting_components():
    euclid = TensorIndexType('Euclidean')
    euclid.data = [1, 1]
    i1, i2, _ = tensor_indices('i1:4', euclid)

    a, b, c, d = symbols('a b c d', commutative=False)
    V1 = tensorhead('V1', [euclid] * 2, [[1]]*2)
    V1.data = [[a, b], (c, d)]
    V2 = tensorhead('V2', [euclid] * 2, [[1]]*2)
    V2.data = [[a, c], [b, d]]

    vtp = V1(i1, i2) * V2(-i2, -i1)

    assert vtp.data == a**2 + b**2 + c**2 + d**2
    assert vtp.data != a**2 + 2*b*c + d**2

    Vc = (b * V1(i1, -i1)).data
    assert Vc.expand() == b * a + b * d


def test_valued_non_diagonal_metric():
    (_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, n0, _,
     _, NA, _, _, _, _, ndm_matrix, *_) = _get_valued_base_test_variables()

    mmatrix = Matrix(ndm_matrix)
    assert (NA(n0)*NA(-n0)).data == (NA(n0).get_matrix().T * mmatrix * NA(n0).get_matrix())[0, 0]


def test_valued_assign_numpy_ndarray():
    (A, _, AB, _, _, _, E, px, py, pz, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, *_) = _get_valued_base_test_variables()

    # this is needed to make sure that a numpy.ndarray can be assigned to a
    # tensor.
    arr = [E+1, px-1, py, pz]
    A.data = numpy.array(arr)
    for i in range(4):
        assert A(i0).data[i] == arr[i]

    qx, qy, qz = symbols('qx qy qz')
    A(-i0).data = numpy.array([E, qx, qy, qz])
    for i in range(4):
        assert A(i0).data[i] == [E, -qx, -qy, -qz][i]
        assert A.data[i] == [E, -qx, -qy, -qz][i]

    # test on multi-indexed tensors.
    random_4x4_data = [[(i**3 - 3*i**2) % (j + 7) for i in range(4)] for j in range(4)]
    AB(-i0, -i1).data = random_4x4_data
    for i in range(4):
        for j in range(4):
            assert AB(i0, i1).data[i, j] == random_4x4_data[i][j]*(-1 if i else 1)*(-1 if j else 1)
            assert AB(-i0, i1).data[i, j] == random_4x4_data[i][j]*(-1 if j else 1)
            assert AB(i0, -i1).data[i, j] == random_4x4_data[i][j]*(-1 if i else 1)
            assert AB(-i0, -i1).data[i, j] == random_4x4_data[i][j]

    AB(-i0, i1).data = random_4x4_data
    for i in range(4):
        for j in range(4):
            assert AB(i0, i1).data[i, j] == random_4x4_data[i][j]*(-1 if i else 1)
            assert AB(-i0, i1).data[i, j] == random_4x4_data[i][j]
            assert AB(i0, -i1).data[i, j] == random_4x4_data[i][j]*(-1 if i else 1)*(-1 if j else 1)
            assert AB(-i0, -i1).data[i, j] == random_4x4_data[i][j]*(-1 if j else 1)


def test_valued_metric_inverse():
    (_, _, _, _, _, Lorentz, _, _, _, _, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, *_) = _get_valued_base_test_variables()

    # let's assign some fancy matrix, just to verify it:
    # (this has no physical sense, it's just testing diofant);
    # it is symmetrical:
    md = [[2, 2, 2, 1], [2, 3, 1, 0], [2, 1, 2, 3], [1, 0, 3, 2]]
    with pytest.raises(ValueError):
        Lorentz.data = [[[1, 2], [1, 2]], [[3, 4], [3, 4]]]
    with pytest.raises(ValueError):
        Lorentz.data = [[1, 2]]
    with pytest.raises(ValueError):
        Lorentz.data = [[1]]
    with pytest.raises(ValueError):
        Lorentz.data = [1]
    with pytest.raises(ValueError):
        Lorentz.data = [[1], [2]]
    Lorentz.data = md
    m = Matrix(md)
    metric = Lorentz.metric
    minv = m.inv()

    meye = eye(4)

    # the Kronecker Delta:
    KD = Lorentz.get_kronecker_delta()

    for i in range(4):
        for j in range(4):
            assert metric(i0, i1).data[i, j] == m[i, j]
            assert metric(-i0, -i1).data[i, j] == minv[i, j]
            assert metric(i0, -i1).data[i, j] == meye[i, j]
            assert metric(-i0, i1).data[i, j] == meye[i, j]
            assert metric(i0, i1)[i, j] == m[i, j]
            assert metric(-i0, -i1)[i, j] == minv[i, j]
            assert metric(i0, -i1)[i, j] == meye[i, j]
            assert metric(-i0, i1)[i, j] == meye[i, j]

            assert KD(i0, -i1)[i, j] == meye[i, j]


def test_valued_canon_bp_swapaxes():
    (A, B, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
     _, _, _, _, _, _, _, i0, i1, i2, *_) = _get_valued_base_test_variables()

    e1 = A(i1)*A(i0)
    e1.data[0, 1] = 44
    e2 = e1.canon_bp()
    assert e2 == A(i0)*A(i1)
    for i in range(4):
        for j in range(4):
            assert e1[i, j] == e2[j, i]
    o1 = B(i2)*A(i1)*B(i0)
    o2 = o1.canon_bp()
    for i in range(4):
        for j in range(4):
            for k in range(4):
                assert o1[i, j, k] == o2[j, i, k]


def test_contract_automatrix_and_data():
    L = TensorIndexType('L')
    S = TensorIndexType('S')
    G = tensorhead('G', [L, S, S], [[1]]*3, matrix_behavior=True)

    def G_data():
        G.data = [[[1]]]
    pytest.raises(ValueError, G_data)
    L.data = [1, -1]
    pytest.raises(ValueError, G_data)
    S.data = [[1, 0], [0, 2]]
    G.data = [
        [[1, 2],
         [3, 4]],
        [[5, 6],
         [7, 8]]
    ]
    m0, *_ = tensor_indices('m0:3', L)
    s0, s1, s2 = tensor_indices('s0:3', S)

    assert (G(-m0).data == numpy.array([
        [[1, 4],
         [3, 8]],
        [[-5, -12],
         [-7, -16]]])).all()

    c1 = G(m0, s0, -s1)*G(-m0, s1, -s2)
    c2 = G(m0) * G(-m0)

    assert (c1.data == c2.data).all()

    del L.data
    del S.data
    del G.data
    assert L.data is None
    assert S.data is None
    assert G.data is None


def test_valued_components_with_wrong_symmetry():
    IT = TensorIndexType('IT', dim=3)
    IT.data = [1, 1, 1]
    A_nosym = tensorhead('A', [IT]*2, [[1]]*2)
    A_sym = tensorhead('A', [IT]*2, [[1]*2])
    A_antisym = tensorhead('A', [IT]*2, [[2]])

    mat_nosym = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    mat_sym = mat_nosym + mat_nosym.T
    mat_antisym = mat_nosym - mat_nosym.T

    A_nosym.data = mat_nosym
    A_nosym.data = mat_sym
    A_nosym.data = mat_antisym

    def assign(A, dat):
        A.data = dat

    A_sym.data = mat_sym
    pytest.raises(ValueError, lambda: assign(A_sym, mat_nosym))
    pytest.raises(ValueError, lambda: assign(A_sym, mat_antisym))

    A_antisym.data = mat_antisym
    pytest.raises(ValueError, lambda: assign(A_antisym, mat_sym))
    pytest.raises(ValueError, lambda: assign(A_antisym, mat_nosym))

    A_sym.data = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    A_antisym.data = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]


def test_sympyissue_10972():
    Lorentz = TensorIndexType('Lorentz', metric=False, dummy_fmt='i', dim=2)
    Lorentz.data = [-1, 1]

    mu, nu, alpha, beta = tensor_indices('\\mu, \\nu, \\alpha, \\beta',
                                         Lorentz)

    Vec = TensorType([Lorentz], tensorsymmetry([1]))
    A2 = TensorType([Lorentz] * 2, tensorsymmetry([2]))

    u = Vec('u')
    u.data = [1, 0]

    F = A2('F')
    F.data = [[0, 1],
              [-1, 0]]

    mul_1 = F(mu, alpha) * u(-alpha) * F(nu, beta) * u(-beta)
    assert (mul_1.data == numpy.array([[0, 0], [0, 1]])).all()

    mul_2 = F(mu, alpha) * F(nu, beta) * u(-alpha) * u(-beta)
    assert (mul_2.data == mul_1.data).all()

    assert ((mul_1 + mul_1).data == 2 * mul_1.data).all()


def test_TensMul_data():
    Lorentz = TensorIndexType('Lorentz', metric=False, dummy_fmt='L', dim=4)
    Lorentz.data = [-1, 1, 1, 1]

    mu, nu, alpha, beta = tensor_indices('\\mu, \\nu, \\alpha, \\beta',
                                         Lorentz)

    Vec = TensorType([Lorentz], tensorsymmetry([1]))
    A2 = TensorType([Lorentz] * 2, tensorsymmetry([2]))

    u = Vec('u')
    u.data = [1, 0, 0, 0]

    F = A2('F')
    Ex, Ey, Ez, Bx, By, Bz = symbols('E_x E_y E_z B_x B_y B_z')
    F.data = [[0, Ex, Ey, Ez],
              [-Ex, 0, Bz, -By],
              [-Ey, -Bz, 0, Bx],
              [-Ez, By, -Bx, 0]]

    E = F(mu, nu) * u(-nu)

    assert ((E(mu) * E(nu)).data ==
            numpy.array([[0, 0, 0, 0],
                         [0, Ex ** 2, Ex * Ey, Ex * Ez],
                         [0, Ex * Ey, Ey ** 2, Ey * Ez],
                         [0, Ex * Ez, Ey * Ez, Ez ** 2]])).all()

    assert ((E(mu) * E(nu)).canon_bp().data == (E(mu) * E(nu)).data).all()

    assert ((F(mu, alpha) * F(beta, nu) * u(-alpha) * u(-beta)).data ==
            - (E(mu) * E(nu)).data).all()
    assert ((F(alpha, mu) * F(beta, nu) * u(-alpha) * u(-beta)).data ==
            (E(mu) * E(nu)).data).all()

    S2 = TensorType([Lorentz] * 2, tensorsymmetry([1] * 2))
    g = S2('g')
    g.data = Lorentz.data

    # tensor 'perp' is orthogonal to vector 'u'
    perp = u(mu) * u(nu) + g(mu, nu)

    mul_1 = u(-mu) * perp(mu, nu)
    assert (mul_1.data == numpy.array([0, 0, 0, 0])).all()

    mul_2 = u(-mu) * perp(mu, alpha) * perp(nu, beta)
    assert (mul_2.data == numpy.zeros(shape=(4, 4, 4))).all()

    Fperp = perp(mu, alpha) * perp(nu, beta) * F(-alpha, -beta)
    assert (Fperp.data[0, :] == numpy.array([0, 0, 0, 0])).all()
    assert (Fperp.data[:, 0] == numpy.array([0, 0, 0, 0])).all()

    mul_3 = u(-mu) * Fperp(mu, nu)
    assert (mul_3.data == numpy.array([0, 0, 0, 0])).all()


def test_sympyissue_11020():
    Lorentz = TensorIndexType('Lorentz', metric=False, dummy_fmt='i', dim=2)
    Lorentz.data = [-1, 1]

    a, b, c, d = tensor_indices('a, b, c, d', Lorentz)
    i0, _ = tensor_indices('i_0:2', Lorentz)

    Vec = TensorType([Lorentz], tensorsymmetry([1]))
    S2 = TensorType([Lorentz] * 2, tensorsymmetry([1] * 2))

    # metric tensor
    g = S2('g')
    g.data = Lorentz.data

    u = Vec('u')
    u.data = [1, 0]

    add_1 = g(b, c) * g(d, i0) * u(-i0) - g(b, c) * u(d)
    assert (add_1.data == numpy.zeros(shape=(2, 2, 2))).all()
    # Now let us replace index `d` with `a`:
    add_2 = g(b, c) * g(a, i0) * u(-i0) - g(b, c) * u(a)
    assert (add_2.data == numpy.zeros(shape=(2, 2, 2))).all()

    # some more tests
    # perp is tensor orthogonal to u^\mu
    perp = u(a) * u(b) + g(a, b)
    mul_1 = u(-a) * perp(a, b)
    assert (mul_1.data == numpy.array([0, 0])).all()

    mul_2 = u(-c) * perp(c, a) * perp(d, b)
    assert (mul_2.data == numpy.zeros(shape=(2, 2, 2))).all()
