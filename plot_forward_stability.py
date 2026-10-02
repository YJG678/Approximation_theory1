"""Plot the forward-stability figure from forward_stability_data.npz."""
import numpy as np
from plotstyle import *

D = np.load("forward_stability_data.npz")
th, Ns = D["thetas"], D["Ns"]
K = int(D["K"])

seq = [mpl.colormaps["Blues"](v) for v in np.linspace(0.38, 0.95, len(Ns))]
cat, mk = CAT, MK


def line_label(ax, xs, ys, x_at, text, above=True):
    """Place text parallel to a curve (log-log), just above it."""
    ax.figure.canvas.draw()
    i = np.searchsorted(xs, x_at)
    p0 = ax.transData.transform((xs[i - 1], ys[i - 1]))
    p1 = ax.transData.transform((xs[i + 1], ys[i + 1]))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    y_at = np.interp(np.log(x_at), np.log(xs), np.log(ys))
    ax.annotate(text, (x_at, np.exp(y_at)), xytext=(0, 3.5 if above else -3.5),
                textcoords="offset points", rotation=ang, rotation_mode="anchor",
                ha="left", va="bottom" if above else "top", fontsize=7.5, color="#0b0b0b")


thd = np.logspace(-4, np.log10(0.9), 400)
fig, axs = plt.subplots(2, 2, figsize=(6.4, 5.0), constrained_layout=True)

# ---------------- (a) worst-case finite defect vs theta
ax = axs[0, 0]
for i, N in enumerate(Ns):
    ax.loglog(th, D["dA"][i], color=seq[i], lw=1.0)
bA = np.minimum(1, thd / (1 - thd))
ax.loglog(thd, bA, **BOUND)
ax.set_xlabel(r"$\theta=|\varepsilon|\,\|h\|_\infty$")
ax.set_ylabel(rf"$\delta_{{N,{K}}}(\mu_\varepsilon)$")
ax.set_title(r"(a) worst-case defect, smooth $h$", loc="left")
guide(ax, thd[:300], 1e-2, 4e-6, 1, r"slope $1$", where=0.35, above=False)
style(ax)

# ---------------- (b) uniformity in N for high-frequency h
ax = axs[0, 1]
NsB, thB = D["NsB"], float(D["thB"])
for k, (l, row) in enumerate(zip(D["ells"], D["dB"])):
    ax.plot(NsB, row, color=cat[k], lw=1.0, marker=mk[k], ms=3.2, mew=0,
            label=rf"$h=\operatorname{{Re}}(x+\mathrm iy)^{{{l}}}$")
ax.axhline(1 / (1 - thB), **BOUND)
ax.text(NsB[-1], 1 / (1 - thB) * 0.93, r"$1/(1-\theta)$ (Thm 4.1)", fontsize=7.5, ha="right", va="top")
ax.set_yscale("log")
ax.set_ylim(1e-2, 2.0)
ax.set_xlabel(r"output degree $N$")
ax.set_ylabel(rf"$\delta_{{N,{K}}}(\mu_\varepsilon)/\theta$")
ax.set_title(rf"(b) uniformity in $N$ at $\theta={thB:g}$", loc="left")
ax.legend(loc="lower right", frameon=False, handlelength=1.6, borderaxespad=0.2)
ax.set_xticks(NsB[::2])
style(ax)

# ---------------- (c) excess squared error vs theta
ax = axs[1, 0]
for i, N in enumerate(Ns):
    ax.loglog(th, D["exc"][i], color=seq[i], lw=1.0)
bC = thd**2 / (1 - thd)
ax.loglog(thd, bC, **BOUND)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel(r"$\bigl(E_N^{\mu_\varepsilon}(u)^2-E_{N,\mathrm{vec}}^{\mu_\varepsilon}(u)^2\bigr)\big/E_N^\nu(u)^2$", fontsize=8)
ax.set_title(r"(c) excess squared error, fixed $u$", loc="left")
guide(ax, thd[:300], 1e-2, 1e-11, 2, r"slope $2$", where=0.35, above=False)
style(ax)

# ---------------- (d) linearization remainder vs theta
ax = axs[1, 1]
for i, N in enumerate(Ns):
    ax.loglog(th, D["rem"][i], color=seq[i], lw=1.0)
bD = 2 * thd**2 / (1 - thd)
ax.loglog(thd, bD, **BOUND)
ax.set_xlabel(r"$\theta$")
ax.set_ylabel(r"$\|\mathcal C_N^{\mu_\varepsilon}u-\varepsilon\mathscr L_N[h]u\|_\nu\big/E_N^\nu(u)$", fontsize=8)
ax.set_title(r"(d) linearization remainder", loc="left")
guide(ax, thd[:300], 1e-2, 3e-9, 2, r"slope $2$", where=0.35, above=False)
style(ax)

# shared colour key for N (panels a, c, d)
sm = mpl.cm.ScalarMappable(cmap=mpl.colors.ListedColormap(seq),
                           norm=mpl.colors.BoundaryNorm(np.arange(Ns[0] - .5, Ns[-1] + 1), len(Ns)))
cb = fig.colorbar(sm, ax=axs[:, 0].tolist() + [axs[1, 1]], location="right",
                  shrink=0.45, aspect=18, pad=0.02, ticks=Ns)
cb.set_label(r"output degree $N$ in (a), (c), (d)", fontsize=8)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=7, width=0.4)

line_label(axs[0, 0], thd, bA, 1.4e-4, r"$\min\{1,\theta/(1-\theta)\}$ (Thm 4.1)")
line_label(axs[1, 0], thd, bC, 1.4e-4, r"$\theta^2/(1-\theta)$ (4.7)")
line_label(axs[1, 1], thd, bD, 1.4e-4, r"$2\theta^2/(1-\theta)$ (Prop 4.2)")
finalize(fig)
fig.savefig("fig_forward_stability.pdf")
fig.savefig("fig_forward_stability.png", dpi=220)
print("saved")
