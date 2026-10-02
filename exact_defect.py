"""
Independent exact-arithmetic computation of delta_{3,5}(rho_*)^2 for the
radial counterexample d rho_* = (2/(3 pi)) (1 + |x|^2) dx on the unit disk.

All moments and both Gram matrices H (input gradients) and J (defect) are
exact rationals (sympy).  delta_{N,K}^2 is the largest eigenvalue of the
pencil (J, H), an algebraic number, evaluated with mpmath at 40 digits.
No floating-point quadrature or basis from diskpoly.py is used.

Requires sympy and mpmath (pip install sympy mpmath).
"""
from functools import lru_cache
import sympy as sp
from sympy import Rational as Q, rf, factorial

x, y = sp.symbols("x y")


def disk_mom(a, b):
    """Moments of the uniform probability measure on the unit disk."""
    if a % 2 or b % 2:
        return sp.Integer(0)
    i, j = a // 2, b // 2
    return rf(Q(1, 2), i) * rf(Q(1, 2), j) / factorial(i + j + 1)


class Measure:
    def __init__(self, mom):
        self.mom = lru_cache(maxsize=None)(mom)

    def integrate(self, p):
        p = sp.Poly(sp.expand(p), x, y)
        return sum(c * self.mom(a, b) for (a, b), c in p.terms())


def density_measure(base, w):
    terms = sp.Poly(sp.expand(w), x, y).terms()
    return Measure(lambda a, b: sp.Rational(sum(c * base.mom(a + i, b + j)
                                                for (i, j), c in terms)))


def monos(n, const=True):
    return [x**i * y**(k - i) for k in range(0 if const else 1, n + 1)
            for i in range(k, -1, -1)]


def grad(u):
    return [sp.diff(u, x), sp.diff(u, y)]


def vip(mu, F, G):
    return mu.integrate(sum(f * g for f, g in zip(F, G)))


def defect_matrices(mu, N, K):
    """Exact H (gradient Gram on Pi_K / R) and J (defect Gram)."""
    B = monos(K, const=False)
    n = len(B)
    H = sp.Matrix(n, n, lambda i, j: vip(mu, grad(B[i]), grad(B[j])))
    S = monos(N - 1)
    Vb = [[s, 0] for s in S] + [[0, s] for s in S]
    GV = sp.Matrix(len(Vb), len(Vb), lambda i, j: vip(mu, Vb[i], Vb[j]))
    RV = sp.Matrix(len(Vb), n, lambda i, j: vip(mu, Vb[i], grad(B[j])))
    CV = GV.LUsolve(RV)
    GB = [grad(b) for b in monos(N, const=False)]
    GG = sp.Matrix(len(GB), len(GB), lambda i, j: vip(mu, GB[i], GB[j]))
    RG = sp.Matrix(len(GB), n, lambda i, j: vip(mu, GB[i], grad(B[j])))
    CG = GG.LUsolve(RG)
    J = CV.T * GV * CV - CG.T * GG * CG          # ||P grad u||^2 - ||Q grad u||^2
    return H, J


def delta_sq_exact(N=3, K=5, digits=40):
    """delta_{N,K}(rho_*)^2 from exact rational matrices, as an mpmath number."""
    import mpmath as mp
    mp.mp.dps = digits
    sigma = Measure(disk_mom)
    rho = density_measure(sigma, Q(2, 3) * (1 + x**2 + y**2))
    H, J = defect_matrices(rho, N, K)
    L = mp.cholesky(mp.matrix(H.tolist()))
    Li = mp.inverse(L)
    A = Li * mp.matrix(J.tolist()) * Li.T
    return max(mp.eigsy(A)[0])


if __name__ == "__main__":
    import mpmath as mp
    print("delta_{3,5}(rho_*)^2 =", mp.nstr(delta_sq_exact(), 25))
