"""Figure: Sobolev consequences (1.19) under a discontinuous perturbation."""
import numpy as np
from plotstyle import *

D = np.load("e4_sobolev.npz")
th, Ns = D["thetas"], D["Ns"]
seq = [INK2] + CAT[:3]           # theta = 0 in neutral ink, perturbed in hues
lw = [1.3, 1.0, 1.0, 1.0]

fig, axs = plt.subplots(1, 3, figsize=(6.4, 2.35), constrained_layout=True)
Ng = np.geomspace(6, 48, 40)

ax = axs[0]
for i, t in enumerate(th):
    ax.loglog(Ns, D["E"][i], color=seq[i], lw=lw[i], marker=MK[i], ms=2.8, mew=0,
              label=rf"$\theta={t:g}$")
guide(ax, Ng, 10, 2.2e-2, -2, r"$N^{-2}$", where=0.45)
ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$E_N^{\mu}(u)$")
ax.set_title(r"(a) gradient error", loc="left")
ax.legend(loc="lower left", frameon=False, handlelength=1.5, borderaxespad=0.1)
style(ax)

ax = axs[1]
for i, t in enumerate(th):
    if t == 0:
        continue
    ax.loglog(Ns, D["exc"][i] / t**2, color=seq[i], lw=lw[i], marker=MK[i], ms=2.8, mew=0)
guide(ax, Ng, 10, 1.6e-3, -4, r"$N^{-4}$", where=0.45)
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\bigl(E_N^{\mu}(u)^2-E_{N,\mathrm{vec}}^{\mu}(u)^2\bigr)/\theta^2$", fontsize=8)
ax.set_title(r"(b) excess error$/\theta^2$", loc="left")
style(ax)

ax = axs[2]
for i, t in enumerate(th):
    ax.loglog(Ns, D["L2"][i], color=seq[i], lw=lw[i], marker=MK[i], ms=2.8, mew=0)
guide(ax, Ng, 10, 3.5e-4, -3, r"$N^{-3}$", where=0.5, above=False)
guide(ax, Ng, 10, 4.5e-3, -2, r"$N^{-2}$", where=0.5)
ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\|u-R_N^{\mu}u\|_{\mu}$")
ax.set_title(r"(c) $L^2$ error", loc="left")
style(ax)

for _a in axs:
    _a.set_xticks([4, 8, 16, 32]); _a.set_xticklabels([str(v) for v in [4, 8, 16, 32]]); _a.minorticks_off() if False else None
    _a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
finalize(fig)
fig.savefig("fig_sobolev.pdf")
fig.savefig("fig_sobolev.png", dpi=220)
m = Ns >= 12
for key in ("E", "exc", "L2"):
    print(key, [round(np.polyfit(np.log(Ns[m]), np.log(D[key][i, m]), 1)[0], 2) for i in range(len(th))])
print("excess ratios theta^2 check:", D["exc"][2, -1] / D["exc"][1, -1], D["exc"][3, -1] / D["exc"][1, -1])
print("E ratio mu/nu at N max:", D["E"][:, -1] / D["E"][0, -1])
