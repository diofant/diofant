"""Gosper's algorithm for hypergeometric summation."""

from ..core import Dummy, Integer, symbols
from ..core.compatibility import is_sequence
from ..polys import parallel_poly_from_expr
from ..simplify import hypersimp
from ..solvers import solve


def gosper_normal(f, g, n):
    r"""Compute the Gosper's normal form of ``f`` and ``g``.

    Given relatively prime univariate polynomials ``f`` and ``g``,
    rewrite their quotient to a normal form defined as follows:

    .. math::
        \frac{f(n)}{g(n)} = Z \cdot \frac{A(n) C(n+1)}{B(n) C(n)}

    where ``Z`` is an arbitrary constant and ``A``, ``B``, ``C`` are
    monic polynomials in ``n`` with the following properties:

    1. `\gcd(A(n), B(n+h)) = 1 \forall h \in \mathbb{N}`
    2. `\gcd(B(n), C(n+1)) = 1`
    3. `\gcd(A(n), C(n)) = 1`

    This normal form, or rational factorization in other words, is a
    crucial step in Gosper's algorithm and in solving of difference
    equations. It can be also used to decide if two hypergeometric
    terms are similar or not.

    This procedure will return a tuple containing elements of this
    factorization in the form ``(Z*A, B, C)``.

    Examples
    ========

    >>> gosper_normal(4*n + 5, 2*(4*n + 1)*(2*n + 3), n)
    (Poly(1/4, n, domain='QQ'), Poly(n + 3/2, n, domain='QQ'),
     Poly(n + 1/4, n, domain='QQ'))

    """
    (C, p, q), _ = parallel_poly_from_expr((1, f, g), n, field=True)

    a, A = p.LC(), p.monic()
    b, B = q.LC(), q.monic()
    Z = a/b

    J = A.dispersionset(B)

    for i in sorted(J):
        d = A.gcd(B.shift(+i))

        A = A.quo(d)
        B = B.quo(d.shift(-i))

        for j in range(1, i + 1):
            C *= d.shift(-j)

    return A*Z, B, C


def gosper_term(f, n):
    r"""Compute Gosper's hypergeometric term for ``f``.

    Suppose ``f`` is a hypergeometric term such that:

    .. math::
        s_n = \sum_{k=0}^{n-1} f_k

    and `f_k` doesn't depend on `n`. Returns a hypergeometric
    term `g_n` such that `g_{n+1} - g_n = f_n`.

    Examples
    ========

    >>> gosper_term((4*n + 1)*factorial(n)/factorial(2*n + 1), n)
    (-n - 1/2)/(n + 1/4)

    """
    r = hypersimp(f, n)

    if r is None:
        return

    p, q = r.as_numer_denom()

    A, B, C = gosper_normal(p, q, n)
    B = B.shift(-1)

    N = Integer(A.degree())
    M = Integer(B.degree())
    K = Integer(C.degree())

    if (N != M) or (A.LC() != B.LC()):
        D = {K - max(N, M)}
    elif not N:
        D = {K - N + 1, Integer(0)}
    else:
        D = {K - N + 1, (B.coeff_monomial(n**(N - 1)) -
                         A.coeff_monomial(n**(N - 1)))/A.LC()}

    for d in set(D):
        if not d.is_Integer or d < 0:
            D.remove(d)

    if not D:
        return

    d = max(D)

    coeffs = symbols(f'c:{d + 1}', cls=Dummy)
    domain = A.domain.inject(*coeffs)

    x = sum(c*n**i for i, c in enumerate(reversed(coeffs)))
    x = x.as_poly(n, domain=domain)
    H = A*x.shift(1) - B*x - C

    solution = solve(H.coeffs(), coeffs)
    if solution:
        solution = solution[0]

    x = x.subs(solution).subs({_: 0 for _ in coeffs})

    if x != 0:
        return B.as_expr()*x.as_expr()/C.as_expr()


def gosper_sum(f, k):
    r"""Gosper's hypergeometric summation algorithm.

    Given a hypergeometric term ``f`` such that:

    .. math ::
        s_n = \sum_{k=0}^{n-1} f_k

    and `f(n)` doesn't depend on `n`, returns `g_{n} - g(0)` where
    `g_{n+1} - g_n = f_n`, or ``None`` if `s_n` can not be expressed
    in closed form as a sum of hypergeometric terms.

    Examples
    ========

    >>> gosper_sum((4*k + 1)*factorial(k)/factorial(2*k + 1), (k, 0, n))
    (-factorial(n) + 2*factorial(2*n + 1))/factorial(2*n + 1)

    References
    ==========

    * :cite:`Petkovsek1997AeqB`

    """
    indefinite = False

    if is_sequence(k):
        k, a, b = k
    else:
        indefinite = True

    g = gosper_term(f, k)

    if g is None:
        return

    if indefinite:
        result = f*g
    else:
        fg = (f*g).cancel()
        result = (f + fg).subs({k: b}) - fg.subs({k: a})

    return result.factor()
