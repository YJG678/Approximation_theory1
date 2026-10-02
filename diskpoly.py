"""
Stable high-degree polynomial projections on the unit disk (d = 2).

Basis: generalized Zernike (disk) polynomials, orthogonal for
nu_alpha = c (1 - |x|^2)^alpha dx,
    Z_{n,m}(r, phi) = r^m P_k^{(alpha, m)}(2 r^2 - 1) * {cos m phi, sin m phi},
    k = (n - m)/2,   n = 0..K,  m = n, n-2, ..., >= 0,
ordered by total degree n, so the first dim(Pi_K) columns span Pi_K.
Gradients are analytic.  Because the basis is (close to) orthogonal for the
reference measure, Gram matrices of perturbations with c_- <= w <= c_+ have
condition number <= c_+/c_-; no monomials are ever formed.

Quadrature: Gauss-Jacobi in t = r^2 (weight (1-t)^alpha) times an angular
rule, either equispaced (for trigonometric integrands) or Gauss-Legendre on
user-given angular sectors (for weights with kinks along rays).
"""
import numpy as np
from scipy.special import roots_jacobi, eval_jacobi, roots_legendre


def disk_rule(alpha, n_rad, n_ang, sectors=None, t_breaks=None, radial_in_r=False):
    if radial_in_r:                         # alpha = 0: Gauss rule in r, weight 2r dr
        assert alpha == 0.0 and t_breaks is None
        s, ws = roots_jacobi(n_rad, 0.0, 1.0)
        t = ((1 + s) / 2) ** 2
        wt = ws / ws.sum()
    elif t_breaks is None:
        s, ws = roots_jacobi(n_rad, alpha, 0.0)
        t = (1 + s) / 2
        wt = ws / ws.sum()
    else:                                   # alpha = 0 only: Gauss-Legendre in t per piece
        assert alpha == 0.0
        g, gw = roots_legendre(n_rad)
        t, wt = [], []
        for a, b in zip(t_breaks[:-1], t_breaks[1:]):
            t.append((a + b) / 2 + (b - a) / 2 * g); wt.append((b - a) / 2 * gw)
        t, wt = np.concatenate(t), np.concatenate(wt)
    if sectors is None:
        phi = 2 * np.pi * (np.arange(n_ang) + 0.5) / n_ang
        wp = np.full(n_ang, 1.0 / n_ang)
    else:                                   # Gauss-Legendre on each sector
        g, gw = roots_legendre(n_ang)
        phi, wp = [], []
        for a, b in zip(sectors[:-1], sectors[1:]):
            phi.append((a + b) / 2 + (b - a) / 2 * g)
            wp.append((b - a) / 2 * gw / (2 * np.pi))
        phi, wp = np.concatenate(phi), np.concatenate(wp)
    r = np.sqrt(t)
    R, P = np.meshgrid(r, phi, indexing="ij")
    W = np.outer(wt, wp)
    return dict(r=R.ravel(), phi=P.ravel(), x=(R * np.cos(P)).ravel(),
                y=(R * np.sin(P)).ravel(), w=W.ravel(), alpha=alpha)


def _jac(k, a, b, s):
    return eval_jacobi(k, a, b, s) if k >= 0 else np.zeros_like(s)


class DiskBasis:
    """Values and gradients of generalized Zernike polynomials of degree <= K."""

    def __init__(self, Q, K, mmax=None, select=None):
        """select(m, kind) -> bool, kind in {"c", "s"}, keeps a subset of the
        basis (e.g. one symmetry class); the degree ordering is preserved."""
        r, phi, al = Q["r"], Q["phi"], Q["alpha"]
        s = 2 * r**2 - 1
        cols_v, cols_x, cols_y, deg, freq = [], [], [], [], []
        for n in range(K + 1):
            for m in range(n % 2, n + 1, 2):
                if mmax is not None and m > mmax:
                    continue
                k = (n - m) // 2
                P = _jac(k, al, m, s)
                dP = (k + al + m + 1) / 2 * _jac(k - 1, al + 1, m + 1, s)
                Rm = r**m * P                                   # R(r)
                dR = (m * r**(m - 1) if m else 0.0) * P + r**m * 4 * r * dP
                Rr = r**(m - 1) * P if m else None              # R(r)/r
                trig = [("c", np.cos(m * phi), -m * np.sin(m * phi))]
                if m:
                    trig.append(("s", np.sin(m * phi), m * np.cos(m * phi)))
                trig = [t_ for t_ in trig if select is None or select(m, t_[0])]
                for _, T, dT in trig:
                    v = Rm * T
                    if m:
                        gx = np.cos(phi) * dR * T - np.sin(phi) * Rr * dT
                        gy = np.sin(phi) * dR * T + np.cos(phi) * Rr * dT
                    else:
                        gx = np.cos(phi) * dR * T
                        gy = np.sin(phi) * dR * T
                    nrm = np.sqrt(np.sum(Q["w"] * v**2))
                    cols_v.append(v / nrm); cols_x.append(gx / nrm); cols_y.append(gy / nrm)
                    deg.append(n); freq.append(m)
        self.val = np.array(cols_v).T
        self.gx = np.array(cols_x).T
        self.gy = np.array(cols_y).T
        self.deg = np.array(deg); self.freq = np.array(freq)
        self.K = K

    def cols(self, n, const=True):
        return np.where((self.deg <= n) & ((self.deg >= 1) | const))[0]


# ------------------------------------------------------------------ algebra
def wproj(A, F, w):
    sw = np.sqrt(w)[:, None]
    Qm, _ = np.linalg.qr(sw * A)
    return (Qm @ (Qm.T @ (sw * F))) / sw


def as2d(f):
    return f[:, None] if f.ndim == 1 else f


def P_vec(B, m, F, w):
    A = B.val[:, B.cols(m)]
    F1, F2 = as2d(F[0]), as2d(F[1]); k = F1.shape[1]
    out = wproj(A, np.hstack([F1, F2]), w)
    return out[:, :k], out[:, k:]


def Q_grad(B, N, F, w, basis_for_G=None):
    Bg = basis_for_G or B
    c = Bg.cols(N, const=False)
    A = np.vstack([Bg.gx[:, c], Bg.gy[:, c]])
    F1, F2 = as2d(F[0]), as2d(F[1])
    out = wproj(A, np.vstack([F1, F2]), np.concatenate([w, w]))
    n = w.size
    return out[:n], out[n:]


def defect(B, N, F, w):
    PF = P_vec(B, N - 1, F, w)
    QF = Q_grad(B, N, F, w)
    return (PF[0] - QF[0], PF[1] - QF[1]), PF, QF


def vnorm(F, w):
    return np.sqrt(np.sum(w[:, None] * (as2d(F[0])**2 + as2d(F[1])**2), axis=0))


def delta_all_K(B, N, w, Ks, Binput=None):
    """delta_{N,K}(mu) for every K in Ks from one nested QR.

    Binput (optional) restricts the input space (e.g. to angular frequencies
    that can interact with the defect for radial measures).  B supplies
    Pi_{N-1} and G_N and must contain them."""
    Bi = Binput or B
    c = Bi.cols(max(Ks), const=False)            # graded by degree
    n = w.size
    A = np.vstack([Bi.gx[:, c], Bi.gy[:, c]])
    ww = np.concatenate([w, w]); sw = np.sqrt(ww)[:, None]
    Qm, Rm = np.linalg.qr(sw * A)
    G = Qm / sw
    Cf, _, _ = defect(B, N, (G[:n], G[n:]), w)
    M = sw * np.vstack([Cf[0], Cf[1]])
    degc = Bi.deg[c]
    Gm = M.T @ M                                  # small Gram matrix of the defects
    out = []
    for K in Ks:
        j = np.searchsorted(degc, K, side="right")
        out.append(np.sqrt(max(np.linalg.eigvalsh(Gm[:j, :j])[-1], 0.0)) if j else 0.0)
    return np.array(out)
