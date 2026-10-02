"""Figure: finite-degree truncation (Theorem 5.6)."""
import numpy as np
from plotstyle import *

R = np.load("kref_circle_K300.npz"); Kr, d2r = R["Ks"], R["d2_3"]          # K_ref = 300
P = np.load("kref_plus_K80.npz"); Kp = P["Ks"]                            # K_ref = 80
D = np.load("e2_nonradial_K56.npz"); Kn = D["Ks"]                         # analytic density

# model-free gaps: reference = largest computed K, plotted only for K <= Kref/2
def gap(Ks, d2, kmin, frac=0.5):
    Kref = Ks[-1]
    m = (Ks >= kmin) & (Ks <= frac * Kref)
    return Ks[m], d2[-1] - d2[m]

series = [
    (r"$w_\circ$, $N=3$", *gap(Kr, d2r, 5)),
    (r"$w_+$, $N=2$", *gap(Kp, P["d2_2"], 4)),
    (r"$w_+$, $N=3$", *gap(Kp, P["d2_3"], 4)),
]
slopes = {}
for lab, K, g in series:
    m = (K >= 10) & (K <= (K.max() * 2) / 3)   # fit range 10 <= K <= K_ref/3
    slopes[lab] = np.polyfit(np.log(K[m]), np.log(g[m]), 1)[0]
print("fitted slopes (10<=K<=K_ref/3):", slopes)

fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.55), constrained_layout=True,
                        gridspec_kw=dict(width_ratios=[1.35, 1]))
ax = axs[0]
for k, (lab, K, g) in enumerate(series):
    ax.loglog(K, g, color=CAT[k], lw=1.0, marker=MK[k], ms=2.6, mew=0,
              markevery=max(1, len(K) // 14), label=lab)
Kg = np.geomspace(6, 150, 50)
guide(ax, Kg, 10, 4e-4, -2, r"$K^{-2}$ (a priori)", where=0.55)
guide(ax, Kg, 10, 2.2e-5, -3, r"$K^{-3}$", where=0.55, above=False)
ax.set_xlabel(r"input degree $K$")
ax.set_ylabel(r"$\delta_N(\mu)^2-\delta_{N,K}(\mu)^2$")
ax.set_title(r"(a) Lipschitz densities", loc="left")
ax.legend(loc="lower left", frameon=False, handlelength=1.6, borderaxespad=0.2)
style(ax)

ax = axs[1]
K, g = Kn[Kn <= 26], (D["smooth_2"][-1] - D["smooth_2"])[Kn <= 26]
g = np.maximum(g, 1e-16)
ax.semilogy(K, g, color=CAT[0], lw=1.0, marker="o", ms=2.6, mew=0)
ax.axhline(1e-15, **GUIDE)
ax.text(K[-1], 1.6e-15, r"round-off", fontsize=7, color=INK2, ha="right", va="bottom")
ax.set_xlabel(r"input degree $K$")
ax.set_ylabel(r"$\delta_2(\mu)^2-\delta_{2,K}(\mu)^2$")
ax.set_title(r"(b) analytic density \eqref{eq:validated-example-density}" if False else r"(b) analytic density (6.1)", loc="left")
ax.set_ylim(1e-16, 1e-2)
style(ax)

for _a in [axs[0]]:
    _a.set_xticks([5, 10, 20, 50, 150]); _a.set_xticklabels([str(v) for v in [5, 10, 20, 50, 150]]); _a.minorticks_off() if False else None
    _a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
finalize(fig)
fig.savefig("fig_truncation.pdf")
fig.savefig("fig_truncation.png", dpi=220)
np.save("truncation_slopes.npy", np.array([slopes[s[0]] for s in series]))
# the exact series drawn in the figure, so the verifier can check them
np.savez("fig_truncation_data.npz",
         K_circle=series[0][1], g_circle=series[0][2], K_plus2=series[1][1], g_plus2=series[1][2],
         K_plus3=series[2][1], g_plus3=series[2][2], K_analytic=K, g_analytic=g)

# ---- verify the numbers this script feeds into Section 6.2 (stored data only)
from verify_numerics import main as _verify
_verify(["--section", "6.2", "--quick", "--no-files"])
