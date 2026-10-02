"""
Sobolev consequences (1.19): perturbed weights d mu = (1 + eps h) d nu_0 on the
unit disk with a discontinuous h (checkerboard sign(x1 x2)) and an input of
limited regularity, u = |x1|^{5/2}(1 + x2/2), which lies in H^m for m < 3.

For each theta and N we record
  E_N^mu(u)      = ||grad(u - R_N^mu u)||_mu,
  excess         = E_N^mu(u)^2 - E_{N,vec}^mu(u)^2 = ||C_N^mu u||_mu^2,
  L2 error       = ||u - R_N^mu u||_mu.
"""
import numpy as np, time, json, sys
from diskpoly import *

args = [a for a in sys.argv[1:] if not a.startswith("--")]
COARSE = "--coarse" in sys.argv          # robustness rerun quoted in Sec. 6.3
Nmax = int(args[0]) if args else 48
sect = [0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi]
if COARSE:   # Gauss rule in |x|^2 and fewer angular nodes
    Q = disk_rule(0.0, Nmax + 16, Nmax // 2 + 16, sectors=sect)
else:        # Gauss rule in |x|
    Q = disk_rule(0.0, Nmax + 24, Nmax + 24, sectors=sect, radial_in_r=True)
x, y, w0 = Q["x"], Q["y"], Q["w"]
B = DiskBasis(Q, Nmax)

beta = 2.5
ax = np.abs(x)
u = ax**beta * (1 + y / 2)
gu = (beta * np.sign(x) * ax**(beta - 1) * (1 + y / 2), ax**beta / 2)
h = np.sign(x * y)                                   # nu_0(h) = 0, ||h||_inf = 1

def energy_proj_vals(N, w):
    c = B.cols(N, const=False)
    A = np.vstack([B.gx[:, c], B.gy[:, c]])
    sw = np.sqrt(np.concatenate([w, w]))[:, None]
    Qm, Rm = np.linalg.qr(sw * A)
    F = np.concatenate(gu)[:, None]
    coef = np.linalg.solve(Rm, Qm.T @ (sw * F))
    vals = (B.val[:, c] @ coef)[:, 0]
    return vals + (np.sum(w * u) - np.sum(w * vals)) / np.sum(w)

thetas = [0.0, 0.25, 0.5, 0.75]
Ns = np.unique(np.round(np.geomspace(4, Nmax, 14)).astype(int))
res = {k: np.zeros((len(thetas), len(Ns))) for k in ("E", "Evec", "exc", "L2")}
t0 = time.time()
for i, th in enumerate(thetas):
    w = w0 * (1 + th * h)
    for j, N in enumerate(Ns):
        C, PF, QF = defect(B, N, gu, w)
        E = vnorm((gu[0][:, None] - QF[0], gu[1][:, None] - QF[1]), w)[0]
        Ev = vnorm((gu[0][:, None] - PF[0], gu[1][:, None] - PF[1]), w)[0]
        res["E"][i, j] = E; res["Evec"][i, j] = Ev
        res["exc"][i, j] = vnorm(C, w)[0] ** 2
        Ru = energy_proj_vals(N, w)
        res["L2"][i, j] = np.sqrt(np.sum(w * (u - Ru) ** 2))
    print("theta", th, "done %.1fs" % (time.time() - t0), flush=True)

np.savez("e4_sobolev_coarse.npz" if COARSE else "e4_sobolev.npz", thetas=thetas, Ns=Ns, beta=beta, **res)
for k in res:
    for i, th in enumerate(thetas):
        m = Ns >= 12
        sl = np.polyfit(np.log(Ns[m]), np.log(res[k][i, m]), 1)[0]
        print(k, th, "slope N>=12: %.3f" % sl, " values", " ".join("%.2e" % v for v in res[k][i, ::4]))

# ---- verify the numbers this script feeds into Section 6.3 (stored data only)
from verify_numerics import main as _verify
_verify(["--section", "6.3", "--quick", "--no-files"])
