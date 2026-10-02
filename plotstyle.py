"""Shared matplotlib style for the manuscript figures (Latin Modern via LaTeX)."""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath,amssymb,mathrsfs}",
    "font.family": "serif",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.labelcolor": "#0b0b0b",
})

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3de"
BOUND = dict(color=INK, lw=1.1, ls=(0, (4, 2)), zorder=5)
GUIDE = dict(color=INK2, lw=0.7, ls=(0, (1.2, 1.6)), zorder=4)
# validated categorical order (blue, orange, aqua, violet); markers give a
# second channel so identity never rests on colour alone
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]
MK = ["o", "s", "^", "D"]


def style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


_PENDING = []


def guide(ax, x, x0, y0, slope, label, where=0.5, dy=2.5, above=True):
    """Dotted reference line y = y0 (x/x0)^slope; its label is placed by finalize()."""
    x = np.asarray(x, float)
    y = y0 * (x / x0) ** slope
    ax.plot(x, y, **GUIDE)
    _PENDING.append((ax, x, y, label, where, dy, above))


def finalize(fig):
    """Place guide labels parallel to their lines, after the layout is fixed."""
    fig.canvas.draw()
    fig.set_layout_engine("none")
    for ax, x, y, label, where, dy, above in _PENDING:
        i = int(where * (len(x) - 1))
        p0 = ax.transData.transform((x[max(i - 1, 0)], y[max(i - 1, 0)]))
        p1 = ax.transData.transform((x[min(i + 1, len(x) - 1)], y[min(i + 1, len(x) - 1)]))
        ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
        sgn = 1 if above else -1
        off = (-np.sin(np.radians(ang)) * dy * sgn, np.cos(np.radians(ang)) * dy * sgn)
        ax.annotate(label, (x[i], y[i]), xytext=off, textcoords="offset points",
                    rotation=ang, rotation_mode="anchor", ha="center",
                    va="bottom" if above else "top", fontsize=7.5, color=INK2)
    _PENDING.clear()
