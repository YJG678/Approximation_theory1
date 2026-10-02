"""
Run the forward-stability sweep and save the data to forward_stability_data.npz.

Panels:
 (a) worst-case finite defect delta_{N,K}(mu_eps), K = 24, versus theta,
     for a generic smooth non-radial perturbation h;
 (b) delta_{N,K}(mu_eps)/theta versus N at theta = 0.1 for the
     high-frequency perturbations h_l = Re (x + i y)^l;
 (c) excess squared error E_N^{mu}(u)^2 - E_{N,vec}^{mu}(u)^2, normalized by
     E_N^nu(u)^2, versus theta, for the analytic input u = 1/(2 - x - 0.6 y);
 (d) linearization remainder ||C_N^{mu_eps} u - eps L_N[h] u||_nu / E_N^nu(u).
"""
import numpy as np, time, json
from diskpoly import *

ALPHA = 0.0
K = 24
Qd = disk_rule(ALPHA, 40, 96)
X, Y, W = Qd["x"], Qd["y"], Qd["w"]
B = DiskBasis(Qd, K)
def delta_NK(B, N, K, w):
    return delta_all_K(B, N, w, [K])[0]

# input function u and its gradient
den = 2 - X - 0.6 * Y
gu = (1 / den**2, 0.6 / den**2)

# generic smooth non-radial perturbation, nu(h) = 0, ||h||_inf = 1
def h0f(x, y):
    return np.sin(3 * x + 2 * y + 0.5) + x * y**2 - 0.4 * y
rr = np.sqrt(np.linspace(0, 1, 1500)); pp = np.linspace(0, 2 * np.pi, 3000, endpoint=False)
Rg, Pg = np.meshgrid(rr, pp)
m0 = np.sum(W * h0f(X, Y))
hinf = np.abs(h0f(Rg * np.cos(Pg), Rg * np.sin(Pg)) - m0).max()
h = (h0f(X, Y) - m0) / hinf

checks = {}
checks["nu_mean_h"] = float(np.sum(W * h))
checks["delta_ref_max"] = float(max(delta_NK(B, N, K, W) for N in range(2, 17)))

thetas = np.logspace(-4, np.log10(0.9), 22)
Ns = np.arange(2, 11)

# ---- (a) worst-case finite defect
t0 = time.time()
dA = np.array([[delta_NK(B, N, K, W * (1 + th * h)) for th in thetas] for N in Ns])
print("panel a", time.time() - t0)
# K-convergence check at a few points
checks["K_conv"] = [[int(N), float(th), float(delta_NK(B, N, 20, W * (1 + th * h))),
                     float(delta_NK(B, N, 24, W * (1 + th * h)))] for N in (2, 6, 10) for th in (0.01, 0.5)]

# ---- (b) uniformity in N for high-frequency perturbations
z = X + 1j * Y
ells = [5, 9, 13, 17]
NsB = np.arange(2, 17)
thB = 0.1
dB = np.array([[delta_NK(B, N, K, W * (1 + thB * np.real(z**l))) / thB for N in NsB] for l in ells])
print("panel b", time.time() - t0)

# ---- (c), (d) fixed input u
exc = np.zeros((len(Ns), len(thetas)))
exc_direct = np.zeros_like(exc)
rem = np.zeros_like(exc)
EnuN = np.zeros(len(Ns))
for i, N in enumerate(Ns):
    Pnu = P_vec(B, N - 1, gu, W)
    e = (gu[0] - Pnu[0][:, 0], gu[1] - Pnu[1][:, 0])
    Enu = vnorm(e, W)[0]; EnuN[i] = Enu
    he = (h * e[0], h * e[1])
    Lh, _, _ = defect(B, N, he, W)                 # L_N[h]u = (P - Q)(h e)
    for j, th in enumerate(thetas):
        wm = W * (1 + th * h)
        Cm, Pm, Qm = defect(B, N, gu, wm)
        exc[i, j] = vnorm(Cm, wm)[0]**2 / Enu**2
        Emu2 = vnorm((gu[0] - Qm[0][:, 0], gu[1] - Qm[1][:, 0]), wm)[0]**2
        Evec2 = vnorm((gu[0] - Pm[0][:, 0], gu[1] - Pm[1][:, 0]), wm)[0]**2
        exc_direct[i, j] = (Emu2 - Evec2) / Enu**2
        rem[i, j] = vnorm((Cm[0] - th * Lh[0], Cm[1] - th * Lh[1]), W)[0] / Enu
mask = thetas >= 1e-2
checks["identity_rel_err_theta>=1e-2"] = float(np.max(np.abs(exc_direct - exc)[:, mask] / exc[:, mask]))
checks["E_nu"] = EnuN.tolist()
print("panels c,d", time.time() - t0)

np.savez("forward_stability_data.npz", thetas=thetas, Ns=Ns, dA=dA, ells=ells, NsB=NsB,
         thB=thB, dB=dB, exc=exc, exc_direct=exc_direct, rem=rem, EnuN=EnuN, K=K, alpha=ALPHA)
json.dump(checks, open("forward_stability_checks.json", "w"), indent=1)
print(json.dumps(checks, indent=1))

# ---- verify the numbers this script feeds into Section 6.1 (stored data only)
from verify_numerics import main as _verify
_verify(["--section", "6.1", "--quick", "--no-files"])
