from .matexpr import MatrixExpr


class Factorization(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: self.arg.shape)


class LofLU(Factorization):
    pass


class UofLU(Factorization):
    pass


class LofCholesky(LofLU):
    pass


class UofCholesky(UofLU):
    pass


class QofQR(Factorization):
    pass


class RofQR(Factorization):
    pass


class EigenVectors(Factorization):
    pass


class EigenValues(Factorization):
    pass


class UofSVD(Factorization):
    pass


class SofSVD(Factorization):
    pass


class VofSVD(Factorization):
    pass


def lu(expr):
    return LofLU(expr), UofLU(expr)


def qr(expr):
    return QofQR(expr), RofQR(expr)


def eig(expr):
    return EigenValues(expr), EigenVectors(expr)


def svd(expr):
    return UofSVD(expr), SofSVD(expr), VofSVD(expr)
