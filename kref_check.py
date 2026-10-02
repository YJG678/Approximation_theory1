"""K_ref robustness for Figure 2: w_+ with K_ref = 80 (D4 symmetry split),
w_circ with K_ref = 300 (radial reduction)."""
import numpy as np, time, sys
from diskpoly import *
which = sys.argv[1] if len(sys.argv) > 1 else "both"
if which not in ("plus", "circle", "both"):
    sys.exit("usage: python kref_check.py [plus|circle|both]")
if which in ("plus", "both"):
    K = 80
    sect = [0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi]
    Q = disk_rule(0.0, K + 16, K + 16, sectors=sect, radial_in_r=True)
    x, y = Q["x"], Q["y"]
    w = Q["w"] * (1 + 0.5 * (np.abs(x) + np.abs(y) - 8 / (3 * np.pi)))
    B = DiskBasis(Q, 4)
    classes = {"A1": lambda m, k: m % 4 == 0 and k == "c", "A2": lambda m, k: m % 4 == 0 and k == "s",
               "B1": lambda m, k: m % 4 == 2 and k == "c", "B2": lambda m, k: m % 4 == 2 and k == "s",
               "E": lambda m, k: m % 2 == 1 and k == "c"}
    Ks = list(range(3, K + 1))
    out = {}
    for c, f in classes.items():
        t0 = time.time(); Bi = DiskBasis(Q, K, select=f)
        for N in (2, 3):
            out[(c, N)] = delta_all_K(B, N, w, Ks, Binput=Bi)
        del Bi
        print(c, "done %.0fs" % (time.time() - t0), flush=True)
    d2 = {N: np.max(np.array([out[(c, N)] for c in classes]), axis=0) ** 2 for N in (2, 3)}
    np.savez("kref_plus_K80.npz", Ks=Ks, d2_2=d2[2], d2_3=d2[3])
if which in ("circle", "both"):
    K = 300
    Q = disk_rule(0.0, 170, 32, t_breaks=[0, 0.5, 1.0])
    t = Q["x"]**2 + Q["y"]**2
    w = Q["w"] * (1 + 0.8 * (np.abs(t - 0.5) - 0.25))
    B = DiskBasis(Q, 4); Bi = DiskBasis(Q, K, mmax=3)
    Ks = list(range(3, K + 1))
    t0 = time.time(); d = delta_all_K(B, 3, w, Ks, Binput=Bi)
    print("radial done %.0fs" % (time.time() - t0))
    np.savez("kref_circle_K300.npz", Ks=Ks, d2_3=d**2)
