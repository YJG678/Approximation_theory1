"""
verify_numerics.py -- check every number stated in Section 6 of the manuscript
against the computed data, and report the result of each check.

Usage
  python verify_numerics.py                  all checks, incl. independent recomputations (~3-5 min)
  python verify_numerics.py --quick          skip recomputations that take more than a few seconds
  python verify_numerics.py --section 6.2    one part only: 6.0 (discretization), 6.1, 6.2, 6.3
  python verify_numerics.py --no-files       print only; do not write the report files

Inputs (written by the experiment scripts; see README.md)
  forward_stability_data.npz, forward_stability_checks.json   run_forward_stability.py
  kref_circle_K300.npz, kref_plus_K80.npz                      kref_check.py
  e2_nonradial_K56.npz                                         e2_nonradial.py 56
  fig_truncation_data.npz, truncation_slopes.npy               plot_truncation.py
  e4_sobolev.npz, e4_sobolev_coarse.npz                        e4_sobolev.py 48 [--coarse]

Outputs
  console table, verification_report.md, verification_report.json; exit code 1 if any FAIL.

Status
  PASS   the statement holds for the computed values exactly as written
  PASS*  holds only after rounding the computed values to the number of
         decimals the paper uses (e.g. 0.2525 reported as "0.25")
  FAIL   does not hold
  SKIP   input file or package missing, or recomputation skipped by --quick
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import Decimal, ROUND_HALF_UP

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# ===================================================================== records
RESULTS = []
ORDER = {"FAIL": 0, "PASS*": 1, "SKIP": 2, "PASS": 3}


def record(cid, section, claim, computed, criterion, status, note=""):
    RESULTS.append(dict(id=cid, section=section, claim=claim, computed=computed,
                        criterion=criterion, status=status, note=note))
    line = "[%-5s] %-7s %s" % (status, cid, claim)
    print(_ascii(line))
    print(_ascii("          computed: %s" % computed))
    if status != "PASS" and note:
        print(_ascii("          note: %s" % note))


def skip(cid, section, claim, why):
    record(cid, section, claim, "-", "-", "SKIP", why)


def _ascii(s):
    """Console output stays ASCII so it works in any Windows code page."""
    rep = {"theta": "theta", "θ": "theta", "δ": "delta", "≤": "<=",
           "≥": ">=", "−": "-", "–": "-", "²": "^2", "≈": "~",
           "ν": "nu", "μ": "mu", "∘": "o", "′": "'", "∞": "inf",
           "×": "x", "‖": "||", "ε": "eps", "ℓ": "l", "ρ": "rho",
           "α": "alpha"}
    for k, v in rep.items():
        s = s.replace(k, v)
    return s.encode("ascii", "replace").decode("ascii")


# ============================================================ claim predicates
def decimals_of(s):
    s = s.strip().lstrip("-+")
    return len(s.split(".")[1]) if "." in s else 0


def rnd(v, d):
    return float(Decimal(repr(float(v))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))


def rounds_to(v, stated):
    """A stated value like '-3.01' is correct if v rounds to it."""
    return rnd(v, decimals_of(stated)) == float(stated)


def range_status(values, lo_s, hi_s):
    """'between lo and hi' with rounded endpoints."""
    lo, hi = float(lo_s), float(hi_s)
    v = np.atleast_1d(np.asarray(values, float))
    if np.all((v >= lo) & (v <= hi)):
        return "PASS"
    # round at the precision of the endpoint being compared against
    dlo, dhi = decimals_of(lo_s), decimals_of(hi_s)
    if all(rnd(t, dlo) >= lo and rnd(t, dhi) <= hi for t in v):
        return "PASS*"
    return "FAIL"


def at_most_status(v, bound_s):
    """'at most b' with a rounded bound."""
    b = float(bound_s)
    if v <= b:
        return "PASS"
    if "e" not in bound_s.lower() and rnd(v, decimals_of(bound_s)) <= b:
        return "PASS*"
    return "FAIL"


def fmt(v, p=4):
    return np.format_float_positional(v, precision=p, unique=False, fractional=False, trim="k") \
        if abs(v) >= 1e-3 and abs(v) < 1e4 else "%.*e" % (p - 1, v)


def slope(xv, yv):
    return float(np.polyfit(np.log(xv), np.log(yv), 1)[0])


def load(name, producer):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        return None, "missing %s (run: %s)" % (name, producer)
    if name.endswith(".npz"):
        return np.load(p), ""
    if name.endswith(".npy"):
        return np.load(p), ""
    with open(p) as f:
        return json.load(f), ""


# ============================================= 6.0 discretization and solver
def section_60(quick):
    S = "6.0"
    from diskpoly import disk_rule, DiskBasis, delta_all_K, defect, vnorm

    # --- (4.22) reproduced
    t0 = time.time()
    Q = disk_rule(0.0, 40, 96); B = DiskBasis(Q, 10)
    t = Q["x"]**2 + Q["y"]**2
    h = 20 * t**3 - 30 * t**2 + 12 * t - 1
    gu = (t**2 + 4 * Q["x"]**2 * t, 4 * Q["x"] * Q["y"] * t)
    worst = 0.0
    for eps in (0.5, -0.3, 0.9):
        C, _, _ = defect(B, 3, gu, Q["w"] * (1 + eps * h))
        val = vnorm(C, Q["w"])[0]; exact = math.sqrt(6) / 70 * abs(eps)
        worst = max(worst, abs(val - exact) / exact)
    digits = -math.log10(max(worst, 1e-300))
    record("6.0-1", S, "The solver reproduces (4.22) [||C_3 u||_nu = sqrt(6)/70 |eps|] to at least twelve digits",
           "max relative error %.1e over eps in {0.5,-0.3,0.9} (%.1f digits)" % (worst, digits),
           "relative error <= 1e-12", "PASS" if worst <= 1e-12 else "FAIL")

    # --- delta_{3,5}(rho_*)^2 = 1.2152578150...e-3, exact rational arithmetic
    try:
        import mpmath as mp
        from exact_defect import delta_sq_exact
        ex = delta_sq_exact(3, 5)
        ex_s = mp.nstr(ex, 20)
        Qr = disk_rule(0.0, 40, 64)
        wr = Qr["w"] * (2 / 3) * (1 + Qr["x"]**2 + Qr["y"]**2)
        num = delta_all_K(DiskBasis(Qr, 12), 3, wr, [5])[0]**2
        rel = abs(num - float(ex)) / float(ex)
        mant = ("%.15e" % float(ex)).split("e")[0].replace(".", "")
        ok_digits = mant.startswith("12152578150")
        record("6.0-2", S, "Exact rational value delta_{3,5}(rho_*)^2 = 1.2152578150... x 10^-3",
               "exact (rational matrices, 40-digit eigenvalue) = %s" % ex_s,
               "leading digits 1.2152578150", "PASS" if ok_digits else "FAIL")
        record("6.0-3", S, "Solver agrees with the exact value to at least twelve digits",
               "solver %.15e, relative error %.1e (%.1f digits)" % (num, rel, -math.log10(rel)),
               "relative error <= 1e-12", "PASS" if rel <= 1e-12 else "FAIL")
    except ImportError as e:
        skip("6.0-2", S, "Exact rational value of delta_{3,5}(rho_*)^2", "needs sympy and mpmath: %s" % e)
        skip("6.0-3", S, "Solver agrees with exact value to twelve digits", "needs sympy and mpmath")

    # --- frequency restriction for radial measures; C_3 u = 0 for m >= 4
    Qr = disk_rule(0.0, 60, 64)
    wr = Qr["w"] * (2 / 3) * (1 + Qr["x"]**2 + Qr["y"]**2)
    Bf = DiskBasis(Qr, 24); B3 = DiskBasis(Qr, 24, mmax=3)
    df = delta_all_K(Bf, 3, wr, [8, 16, 24]); d3 = delta_all_K(Bf, 3, wr, [8, 16, 24], Binput=B3)
    rel = float(np.max(np.abs(df - d3) / df))
    hi = np.where(Bf.freq >= 4)[0]
    Cm, _, _ = defect(Bf, 3, (Bf.gx[:, hi], Bf.gy[:, hi]), wr)
    ratio = float(np.max(vnorm(Cm, wr) / vnorm((Bf.gx[:, hi], Bf.gy[:, hi]), wr)))
    record("6.0-4", S, "Radial mu: restricting inputs to angular frequency <= 3 does not change delta_{3,K}",
           "rho_*: max relative difference %.1e (K = 8, 16, 24); max ||C_3 u||/||grad u|| over %d inputs with m >= 4: %.1e"
           % (rel, len(hi), ratio),
           "relative difference <= 1e-12 and ratio <= 1e-12",
           "PASS" if rel <= 1e-12 and ratio <= 1e-12 else "FAIL")

    # --- K = 300 and K = 80 reached
    R, why = load("kref_circle_K300.npz", "python kref_check.py circle")
    if R is None:
        skip("6.0-5", S, "Radial reduction allows input degrees up to K = 300", why)
    else:
        k = int(R["Ks"][-1])
        record("6.0-5", S, "Radial reduction allows input degrees up to K = 300",
               "largest input degree in kref_circle_K300.npz: %d" % k, "= 300", "PASS" if k == 300 else "FAIL")
    P, why = load("kref_plus_K80.npz", "python kref_check.py plus")
    if P is None:
        skip("6.0-6", S, "Square symmetry: five isotypic classes, K = 80", why)
    else:
        k = int(P["Ks"][-1])
        classes = {"A1": lambda m, s: m % 4 == 0 and s == "c", "A2": lambda m, s: m % 4 == 0 and s == "s",
                   "B1": lambda m, s: m % 4 == 2 and s == "c", "B2": lambda m, s: m % 4 == 2 and s == "s",
                   "E": lambda m, s: m % 2 == 1 and s == "c"}
        # each basis function (m, cos/sin) must fall in exactly one class (E keeps one
        # partner of each odd-frequency pair, which suffices: the pair gives equal defects)
        bad = 0
        for m in range(0, 81):
            for s in (["c"] if m == 0 else ["c", "s"]):
                n = sum(f(m, s) for f in classes.values())
                bad += (n != 1) if not (m % 2 == 1 and s == "s") else (n != 0)
        record("6.0-6", S, "Square symmetry splits the inputs into five classes and allows K = 80",
               "largest input degree %d; classes %s; misassigned basis functions: %d" % (k, "/".join(classes), bad),
               "K = 80, five classes, no misassignment", "PASS" if k == 80 and bad == 0 else "FAIL")
    if quick:
        skip("6.0-7", S, "Five-class split reproduces the full defect", "skipped by --quick")
    else:
        sect = [0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi]
        Q = disk_rule(0.0, 48, 48, sectors=sect, radial_in_r=True)
        w = Q["w"] * (1 + 0.5 * (np.abs(Q["x"]) + np.abs(Q["y"]) - 8 / (3 * np.pi)))
        B4 = DiskBasis(Q, 4); Bf = DiskBasis(Q, 32)
        worst = 0.0
        for N in (2, 3):
            full = delta_all_K(B4, N, w, [10, 20, 32], Binput=Bf)
            per = [delta_all_K(B4, N, w, [10, 20, 32], Binput=DiskBasis(Q, 32, select=f)) for f in classes.values()]
            worst = max(worst, float(np.max(np.abs(np.max(per, axis=0) - full) / full)))
        record("6.0-7", S, "Five-class split reproduces the full defect (w_+, N = 2, 3, K <= 32)",
               "max relative difference %.1e" % worst, "<= 1e-12", "PASS" if worst <= 1e-12 else "FAIL")

    # --- Gegenbauer: delta_{N,60}(nu_alpha) <= 2e-13
    if quick:
        skip("6.0-8", S, "delta_{N,60}(nu_alpha) <= 2e-13 for alpha in {-1/2,0,2}, N in {2,3,6}", "skipped by --quick")
    else:
        vals = []
        for al in (-0.5, 0.0, 2.0):
            Q = disk_rule(al, 70, 150); B = DiskBasis(Q, 60)
            for N in (2, 3, 6):
                vals.append((al, N, delta_all_K(B, N, Q["w"], [60])[0]))
        mx = max(v[2] for v in vals)
        record("6.0-8", S, "delta_{N,60}(nu_alpha) <= 2 x 10^-13 for alpha in {-1/2,0,2}, N in {2,3,6}",
               "max %.2e; " % mx + ", ".join("a=%g,N=%d:%.1e" % v for v in vals),
               "all <= 2e-13", "PASS" if mx <= 2e-13 else "FAIL",
               "round-off level; the margin to 2e-13 depends on the BLAS/LAPACK build")


# ============================================================ 6.1 forward
def section_61(quick):
    S = "6.1"
    D, why = load("forward_stability_data.npz", "python run_forward_stability.py")
    if D is None:
        skip("6.1-*", S, "all forward-stability claims", why)
        return
    th, Ns = D["thetas"], D["Ns"]
    dA, exc, rem, exd = D["dA"], D["exc"], D["rem"], D["exc_direct"]
    K, al, thB = int(D["K"]), float(D["alpha"]), float(D["thB"])

    record("6.1-1", S, "Setting: nu = nu_0, K = 24, N = 2..10 (panels a,c,d), theta = 0.1 (panel b)",
           "alpha = %g, K = %d, N = %d..%d, theta_b = %g, l = %s" % (al, K, Ns[0], Ns[-1], thB,
                                                                     list(map(int, D["ells"]))),
           "as stated", "PASS" if (al == 0 and K == 24 and list(Ns) == list(range(2, 11)) and thB == 0.1) else "FAIL")

    # h: centred, ||h||_inf = 1
    from diskpoly import disk_rule
    Qd = disk_rule(0.0, 40, 96)
    h0 = lambda a, b: np.sin(3 * a + 2 * b + 0.5) + a * b**2 - 0.4 * b
    m0 = float(np.sum(Qd["w"] * h0(Qd["x"], Qd["y"])))
    rr = np.sqrt(np.linspace(0, 1, 1500)); pp = np.linspace(0, 2 * np.pi, 3000, endpoint=False)
    Rg, Pg = np.meshgrid(rr, pp)
    G = np.abs(h0(Rg * np.cos(Pg), Rg * np.sin(Pg)) - m0); hinf = float(G.max())
    mean_h = float(np.sum(Qd["w"] * (h0(Qd["x"], Qd["y"]) - m0) / hinf))
    from scipy.optimize import minimize
    best = hinf
    f = lambda v: -abs(h0(v[0] * np.cos(v[1]), v[0] * np.sin(v[1])) - m0)
    for i in np.argsort(G.ravel())[-40:]:
        r = minimize(f, [Rg.ravel()[i], Pg.ravel()[i]], bounds=[(0, 1), (None, None)], method="L-BFGS-B",
                     options=dict(ftol=1e-15, gtol=1e-12))
        best = max(best, -r.fun)
    dev = best / hinf - 1
    record("6.1-2", S, "h is centred and scaled: nu_0(h) = 0 and ||h||_inf = 1",
           "nu_0(h) = %.1e; ||h||_inf = 1 %+.1e (grid normalisation vs refined maximisation)" % (mean_h, dev),
           "|nu_0(h)| <= 1e-14 and | ||h||_inf - 1 | <= 1e-6",
           "PASS" if abs(mean_h) <= 1e-14 and abs(dev) <= 1e-6 else "FAIL")

    m = th <= 1e-2
    sA = [slope(th[m], dA[i, m]) for i in range(len(Ns))]
    record("6.1-3", S, "delta_{N,24}(mu_eps) has fitted slope 1.000 in theta for theta <= 1e-2, every N = 2..10",
           "slopes " + ", ".join("%.5f" % s for s in sA), "each rounds to 1.000",
           "PASS" if all(rounds_to(s, "1.000") for s in sA) else "FAIL")

    ratio = dA / (th / (1 - th))
    st = range_status([ratio.min(), ratio.max()], "0.011", "0.25")
    record("6.1-4", S, "delta_{N,24} lies between 0.011 and 0.25 times the bound theta/(1-theta)",
           "ratio min %.5f, max %.5f" % (ratio.min(), ratio.max()), "0.011 <= ratio <= 0.25", st,
           "max ratio %.4f exceeds 0.25 but rounds to it at two decimals" % ratio.max() if st == "PASS*" else "")

    dB, ells, NsB = D["dB"], [int(v) for v in D["ells"]], D["NsB"]
    peaks = dB.max(axis=1)
    record("6.1-5", S, "Panel (b): each ratio delta_{N,24}/theta rises to about 0.38",
           "maxima " + ", ".join("l=%d: %.4f" % (l, p) for l, p in zip(ells, peaks)) +
           "; first values " + ", ".join("%.3f" % r[0] for r in dB),
           "each maximum rounds to 0.38 and exceeds the first value",
           "PASS" if all(rounds_to(p, "0.38") and p > r[0] for p, r in zip(peaks, dB)) else "FAIL")
    ok, txt = True, []
    for l, row in zip(ells, dB):
        plateau = NsB[row >= 0.97 * row.max()]
        last = int(plateau.max())
        after = row[NsB > last]
        dec = bool(np.all(np.diff(after) < 0)) if len(after) > 1 else None
        good = abs(last - (l - 1)) <= 1 and (dec is not False)
        ok &= good
        txt.append("l=%d: plateau ends at N=%d%s" % (l, last, "" if dec is None else
                                                    ", then %s" % ("decreasing" if dec else "NOT decreasing")))
    record("6.1-6", S, "Panel (b): stays near the plateau until N ~ l - 1, then decreases",
           "; ".join(txt) + " (plateau = within 3% of the maximum)",
           "plateau end within 1 of l-1; strictly decreasing afterwards",
           "PASS" if ok else "FAIL",
           "for l = 17 the decrease lies beyond the tested range N <= 16")
    record("6.1-7", S, "Panel (b): all ratios remain below the degree-independent envelope 1/(1-theta)",
           "max %.4f vs 1/(1-0.1) = %.4f" % (dB.max(), 1 / (1 - thB)), "max <= 1/(1-theta)",
           "PASS" if dB.max() <= 1 / (1 - thB) else "FAIL")

    sC = [slope(th[m], exc[i, m]) for i in range(len(Ns))]
    sD = [slope(th[m], rem[i, m]) for i in range(len(Ns))]
    record("6.1-8", S, "Excess squared error and linearization remainder both have fitted slope 2.00",
           "excess %.4f..%.4f, remainder %.4f..%.4f (theta <= 1e-2)" % (min(sC), max(sC), min(sD), max(sD)),
           "each rounds to 2.00",
           "PASS" if all(rounds_to(s, "2.00") for s in sC + sD) else "FAIL")

    mk = th >= 1e-2
    relid = float(np.max(np.abs(exd - exc)[:, mk] / exc[:, mk]))
    J, _ = load("forward_stability_checks.json", "python run_forward_stability.py")
    record("6.1-9", S, "Direct difference of best-approximation errors agrees to a relative 7e-6 (theta >= 1e-2)",
           "%.2e (stored check: %s)" % (relid, "%.2e" % J["identity_rel_err_theta>=1e-2"] if J else "n/a"),
           "<= 7e-6", "PASS" if relid <= 7e-6 else "FAIL",
           "cancellation-limited; the value varies with the BLAS build (observed 1.8e-6 to 6.5e-6)")

    bA = np.minimum(1, th / (1 - th)); bC = th**2 / (1 - th); bD = 2 * th**2 / (1 - th)
    rA, rC, rD = (dA / bA).max(), (exc / bC).max(), (rem / bD).max()
    below = rA <= 1 and dB.max() <= 1 / (1 - thB) and rC <= 1 and rD <= 1
    record("6.1-10", S, "All finite-input defects and fixed-input errors remain below their envelopes",
           "max data/bound: (a) %.3f, (b) %.3f, (c) %.2e, (d) %.2e" % (rA, dB.max() * (1 - thB), rC, rD),
           "all <= 1", "PASS" if below else "FAIL")
    mins = {"a": (dA / bA).min(), "c": (exc / bC).min(), "d": (rem / bD).min()}
    orders = {k: -math.log10(v) for k, v in mins.items()}
    record("6.1-11", S, "... in some panels by several orders of magnitude",
           "min data/bound: " + ", ".join("(%s) %.1e = %.1f orders" % (k, mins[k], orders[k]) for k in mins),
           "some panel >= 2 orders below", "PASS" if max(orders.values()) >= 2 else "FAIL")


# ============================================================ 6.2 truncation
def section_62(quick):
    S = "6.2"
    from diskpoly import disk_rule, DiskBasis, delta_all_K
    sect = [0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi]

    # --- the three measures
    Qc = disk_rule(0.0, 80, 64, t_breaks=[0, 0.5, 1.0])
    Qp = disk_rule(0.0, 80, 80, sectors=sect, radial_in_r=True)
    Qa = disk_rule(0.0, 40, 96)
    wc = lambda Q: 1 + 0.8 * (np.abs(Q["x"]**2 + Q["y"]**2 - 0.5) - 0.25)
    wp = lambda Q: 1 + 0.5 * (np.abs(Q["x"]) + np.abs(Q["y"]) - 8 / (3 * np.pi))
    wa = lambda Q: 1 + 0.1 * np.real((Q["x"] + 1j * Q["y"])**5)      # (5.20) after scaling to the unit disk
    worst = 0.0; txt = []
    for name, Q, w in (("w_o", Qc, wc), ("w_+", Qp, wp), ("(5.20)", Qa, wa)):
        W = Q["w"] * w(Q); x, y = Q["x"], Q["y"]
        e = [abs(W.sum() - 1), abs(np.sum(W * x)), abs(np.sum(W * y)), abs(np.sum(W * x * y)),
             abs(np.sum(W * x * x) - np.sum(W * y * y))]
        worst = max(worst, max(e)); txt.append("%s %.0e" % (name, max(e)))
    record("6.2-1", S, "The three measures are probability measures, centred, and isotropic after scalar rescaling",
           "max |mass-1|, |mean|, |m_xy|, |m_xx-m_yy|: " + ", ".join(txt), "all <= 1e-14",
           "PASS" if worst <= 1e-14 else "FAIL")
    lo = min(0.8, 1 - 0.5 * 8 / (3 * np.pi), 0.9)
    hi = max(1.2, 1 + 0.5 * (math.sqrt(2) - 8 / (3 * np.pi)), 1.1)
    r = np.sqrt(np.linspace(0, 1, 801)); p = np.linspace(0, 2 * np.pi, 1601)
    Rg, Pg = np.meshgrid(r, p); grid = dict(x=(Rg * np.cos(Pg)).ravel(), y=(Rg * np.sin(Pg)).ravel())
    vals = np.concatenate([wc(grid), wp(grid), wa(grid)])
    record("6.2-2", S, "Their densities lie between 0.57 and 1.29",
           "exact range [%.5f, %.5f]; sampled range [%.5f, %.5f]" % (lo, hi, vals.min(), vals.max()),
           "0.57 <= w <= 1.29", range_status([lo, hi, vals.min(), vals.max()], "0.57", "1.29"))

    if quick:
        skip("6.2-3", S, "w_o nu_0 is rotationally invariant, so delta_2 = 0", "skipped by --quick")
    else:
        Q = disk_rule(0.0, 60, 96, t_breaks=[0, 0.5, 1.0]); B = DiskBasis(Q, 24)
        d = delta_all_K(B, 2, Q["w"] * wc(Q), [8, 16, 24])
        record("6.2-3", S, "w_o nu_0 is rotationally invariant, so delta_2 = 0",
               "delta_{2,K} = " + ", ".join("%.1e" % v for v in d) + " (K = 8, 16, 24)", "<= 1e-12",
               "PASS" if d.max() <= 1e-12 else "FAIL")

    R, w1 = load("kref_circle_K300.npz", "python kref_check.py circle")
    P, w2 = load("kref_plus_K80.npz", "python kref_check.py plus")
    A, w3 = load("e2_nonradial_K56.npz", "python e2_nonradial.py 56")
    if R is None or P is None or A is None:
        skip("6.2-4..", S, "K_ref values, fitted exponents, analytic density", "; ".join(w for w in (w1, w2, w3) if w))
        return
    kr = (int(R["Ks"][-1]), int(P["Ks"][-1]), int(A["Ks"][-1]))
    record("6.2-4", S, "K_ref = 300 for w_o, 80 for w_+, 56 for the analytic density",
           "K_ref in data: %d, %d, %d" % kr, "= 300, 80, 56", "PASS" if kr == (300, 80, 56) else "FAIL")

    F, wf = load("fig_truncation_data.npz", "python plot_truncation.py")
    if F is None:
        skip("6.2-5", S, "Remainder plotted only for K <= K_ref/2", wf)
    else:
        mx = [int(F["K_circle"].max()), int(F["K_plus2"].max()), int(F["K_plus3"].max())]
        ok = mx[0] <= 150 and mx[1] <= 40 and mx[2] <= 40
        record("6.2-5", S, "Remainder plotted only for K <= K_ref/2 (panel a)",
               "largest plotted K: w_o %d, w_+ N=2 %d, w_+ N=3 %d" % tuple(mx), "<= 150, 40, 40",
               "PASS" if ok else "FAIL")

    def expo(Ks, d2, Kref, hi):
        ref = d2[int(np.where(Ks == Kref)[0][0])]
        m = (Ks >= 10) & (Ks <= hi)
        return slope(Ks[m], ref - d2[m])

    cases = [("w_o, N=3", R["Ks"], R["d2_3"], 300, 200, "-3.01"),
             ("w_+, N=2", P["Ks"], P["d2_2"], 80, 56, "-2.92"),
             ("w_+, N=3", P["Ks"], P["d2_3"], 80, 56, "-3.02")]
    prim, allv, maxchg, rows = [], [], 0.0, []
    for lab, Ks, d2, Kr, Ks_, stated in cases:
        p0 = expo(Ks, d2, Kr, Kr / 3)
        var = [expo(Ks, d2, Kr, Kr / 2), expo(Ks, d2, Ks_, Ks_ / 3), expo(Ks, d2, Ks_, Ks_ / 2)]
        prim.append((lab, p0, stated)); allv += [p0] + var
        maxchg = max(maxchg, max(abs(v - p0) for v in var))
        rows.append("%s: %.4f | Kref/2 %.4f | Kref=%d %.4f | Kref=%d,/2 %.4f" % (lab, p0, var[0], Ks_, var[1], Ks_, var[2]))
    record("6.2-6", S, "Exponents fitted over 10 <= K <= K_ref/3: -3.01 (w_o), -2.92 (w_+, N=2), -3.02 (w_+, N=3)",
           ", ".join("%s %.4f" % (l, v) for l, v, _ in prim), "each rounds to the stated value",
           "PASS" if all(rounds_to(v, s) for _, v, s in prim) else "FAIL")
    T, _ = load("truncation_slopes.npy", "python plot_truncation.py")
    if T is not None:
        d = float(np.max(np.abs(np.array(T) - np.array([v for _, v, _ in prim]))))
        record("6.2-7", S, "The exponents used for the figure (plot_truncation.py) are the same numbers",
               "max |difference| %.1e" % d, "<= 1e-10", "PASS" if d <= 1e-10 else "FAIL")
    st = at_most_status(maxchg, "0.09")
    record("6.2-8", S, "Decreasing K_ref (200 for w_o, 56 for w_+) or fitting to K_ref/2 changes the exponents by at most 0.09",
           "max change %.4f; " % maxchg + " ; ".join(rows), "<= 0.09", st,
           "max change %.4f is above 0.09 but rounds to it" % maxchg if st == "PASS*" else "")
    st = range_status(allv, "-3.10", "-2.87")
    record("6.2-9", S, "Every fitted exponent lies between -3.10 and -2.87",
           "range [%.4f, %.4f] over 12 fits" % (min(allv), max(allv)), "-3.10 <= exponent <= -2.87", st,
           "extremes %.4f and %.4f lie outside but round to the endpoints" % (min(allv), max(allv)) if st == "PASS*" else "")
    record("6.2-10", S, "The remainder decays like K^-3 rather than K^-2",
           "all 12 exponents closer to -3 than to -2", "|e+3| < |e+2|",
           "PASS" if all(abs(v + 3) < abs(v + 2) for v in allv) else "FAIL")

    # --- analytic density
    Ks = A["Ks"]; g2 = A["smooth_2"][-1] - A["smooth_2"]; g3 = A["smooth_3"][-1] - A["smooth_3"]
    sel = (Ks >= 5) & (Ks <= 21)
    c = np.corrcoef(Ks[sel], np.log10(g2[sel]))[0, 1]
    mono = bool(np.all(np.diff(g2[sel]) <= 1e-18))
    rate = np.polyfit(Ks[sel], np.log10(g2[sel]), 1)[0]
    record("6.2-11", S, "Analytic density: the remainder decays geometrically",
           "log10-linear fit K = 5..21: %.2f decades per degree, correlation %.3f, nonincreasing: %s" % (rate, c, mono),
           "correlation <= -0.95 and nonincreasing", "PASS" if c <= -0.95 and mono else "FAIL")
    first = int(Ks[np.argmax(g2 < 1e-12)])
    record("6.2-12", S, "... falls below 1e-12 at K = 20",
           "first K with remainder < 1e-12: %d (K=19: %.2e, K=20: %.2e)" % (first, g2[Ks == 19][0], g2[Ks == 20][0]),
           "= 20", "PASS" if first == 20 else "FAIL")

    # noise floor of the squared remainder: stored N=2 vs N=3 for K >= 7 (equal in exact
    # arithmetic there), and, unless --quick, a rerun with a different quadrature
    floor = float(np.max(np.abs(g2 - g3)[Ks >= 7]))
    how = "N=2 vs N=3 spread for K >= 7"
    if not quick:
        Q = disk_rule(0.0, 64, 64, sectors=sect, radial_in_r=True)
        Bq = DiskBasis(Q, 40); kk = list(range(3, 41))
        d2q = delta_all_K(Bq, 2, Q["w"] * wa(Q), kk)**2
        floor = max(floor, float(np.max(np.abs(d2q - A["smooth_2"][:len(kk)]))))
        how += " and a rerun with another quadrature"
    g22 = float(g2[Ks == 22][0])
    reach = int(Ks[np.argmax(g2 <= 10 * floor)])
    record("6.2-13", S, "... and reaches round-off at K = 22",
           "remainder at K=22: %.2e; numerical floor %.1e (%s); floor first reached at K = %d" % (g22, floor, how, reach),
           "remainder(22) <= 10 x floor", "PASS" if g22 <= 10 * floor else "FAIL",
           "the value at K = 22 is %.0f times the floor and reproducible; the geometric decay "
           "continues to about K = %d" % (g22 / floor, reach))
    d = float(np.max(np.abs(g2 - g3)[(Ks >= 20) & (Ks <= 22)]))
    record("6.2-14", S, "... the two levels (N = 2, 3) give identical values there (K ~ 20-22)",
           "max |difference| %.1e" % d, "<= numerical floor %.1e" % floor, "PASS" if d <= floor else "FAIL")
    pl = (Ks >= 3) & (Ks <= 26)
    diff = np.abs(g2 - g3) / np.maximum(g2, 1e-300)
    badK = [int(k) for k, r_, a_ in zip(Ks[pl], diff[pl], np.abs(g2 - g3)[pl]) if r_ > 1e-6 and a_ > floor]
    record("6.2-15", S, "Caption of Fig. 2(b): the case N = 3 coincides (plotted range K = 3..26)",
           "N=2 and N=3 differ at K = %s (e.g. K=4: %.3e vs %.3e; K=6: %.3e vs %.3e); identical for K >= 7"
           % (badK, g2[Ks == 4][0], g3[Ks == 4][0], g2[Ks == 6][0], g3[Ks == 6][0]) if badK else "identical",
           "relative difference <= 1e-6 or below the floor at every plotted K",
           "FAIL" if badK else "PASS")

    if quick:
        skip("6.2-16", S, "Refining quadrature 40 -> 72 nodes changes delta_{2,24}(w_+ nu_0)^2 by < 1e-15", "skipped by --quick")
    else:
        v = {}
        for n in (40, 72):
            Q = disk_rule(0.0, n, n, sectors=sect, radial_in_r=True)
            v[n] = delta_all_K(DiskBasis(Q, 24), 2, Q["w"] * wp(Q), [24])[0]**2
        dv = abs(v[40] - v[72])
        record("6.2-16", S, "Refining the quadrature from 40 to 72 nodes per direction changes delta_{2,24}(w_+ nu_0)^2 by less than 1e-15",
               "%.15e vs %.15e: change %.2e" % (v[40], v[72], dv), "< 1e-15",
               "PASS" if dv < 1e-15 else "FAIL",
               "the change is %.1e, i.e. about 1e-15 (relative %.0e)" % (dv, dv / v[72]))
    worst = 0.0
    for N in (2, 3):
        a = A["kink_%d" % N]; b = P["d2_%d" % N][:len(a)]
        sel = (A["Ks"] >= 5) & (A["Ks"] <= 56)
        worst = max(worst, float(np.max(np.abs(a - b)[sel] / np.abs(a[sel]))))
    record("6.2-17", S, "w_+ runs with K_ref = 56 and 80 (different quadratures) agree to a relative 1e-12 for 5 <= K <= 56",
           "max relative difference %.1e" % worst, "<= 1e-12", "PASS" if worst <= 1e-12 else "FAIL")


# ============================================================ 6.3 Sobolev
def section_63(quick):
    S = "6.3"
    D, why = load("e4_sobolev.npz", "python e4_sobolev.py 48")
    if D is None:
        skip("6.3-*", S, "all Sobolev claims", why)
        return
    th, Ns = [float(t) for t in D["thetas"]], D["Ns"]
    E, exc, L2 = D["E"], D["exc"], D["L2"]
    from diskpoly import disk_rule
    sect = [0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi]
    Q = disk_rule(0.0, 72, 72, sectors=sect, radial_in_r=True)
    mh = float(np.sum(Q["w"] * np.sign(Q["x"] * Q["y"])))
    record("6.3-1", S, "Setting: h = sign(x1 x2) with nu_0(h) = 0, ||h||_inf = 1; theta in {0,1/4,1/2,3/4}; N up to 48",
           "nu_0(h) = %.1e; theta = %s; N = %s" % (mh, th, list(map(int, Ns))), "as stated",
           "PASS" if abs(mh) <= 1e-14 and th == [0, 0.25, 0.5, 0.75] and int(Ns[-1]) == 48 else "FAIL")
    m = Ns >= 12
    sl = lambda Y: [slope(Ns[m], Y[i, m]) for i in range(len(th))]
    sE, sX, sL = sl(E), sl(exc), sl(L2)

    record("6.3-2", S, "nu_0 best-approximation errors decay at the H^{3-} rates (N^-2 gradient, N^-3 L2)",
           "theta=0: gradient %.3f, L2 %.3f" % (sE[0], sL[0]), "within 0.05 of -2 and -3",
           "PASS" if abs(sE[0] + 2) <= 0.05 and abs(sL[0] + 3) <= 0.05 else "FAIL")
    st = range_status(sE, "-2.05", "-2.01")
    record("6.3-3", S, "Gradient error slopes between -2.01 and -2.05 for all four theta (N >= 12)",
           ", ".join("%.4f" % s for s in sE), "-2.05 <= slope <= -2.01", st,
           "extreme %.4f lies outside but rounds to the endpoint" % min(sE) if st == "PASS*" else "")
    r48 = E[:, -1] / E[0, -1]
    st = range_status(r48, "0.78", "1")
    record("6.3-4", S, "At N = 48 the gradient error lies between 0.78 and 1 times its theta = 0 value",
           ", ".join("%.4f" % v for v in r48), "0.78 <= ratio <= 1", st,
           "smallest ratio %.4f is below 0.78 but rounds to it" % r48.min() if st == "PASS*" else "")
    q = [exc[2, -1] / exc[1, -1], exc[3, -1] / exc[1, -1]]
    record("6.3-5", S, "Excess error ratios at N = 48 (theta = 1/2, 3/4 vs 1/4): 3.82 and 7.57, vs 4 and 9",
           "%.4f, %.4f; exact quadratic scaling (0.5/0.25)^2 = %g, (0.75/0.25)^2 = %g" % (q[0], q[1], (0.5 / 0.25)**2, (0.75 / 0.25)**2),
           "round to 3.82, 7.57", "PASS" if rounds_to(q[0], "3.82") and rounds_to(q[1], "7.57") else "FAIL")
    st = range_status(sX[1:], "-4.86", "-4.77")
    record("6.3-6", S, "Excess error slopes from -4.77 to -4.86, steeper than the limiting order N^{2-2m} = N^-4",
           ", ".join("%.4f" % s for s in sX[1:]), "-4.86 <= slope <= -4.77 and < -4",
           st if all(s < -4 for s in sX[1:]) else "FAIL",
           "extreme %.4f lies outside but rounds to the endpoint" % max(sX[1:]) if st == "PASS*" else "")
    stated = ["-2.96", "-2.85", "-2.71", "-2.65"]
    record("6.3-7", S, "L2 error slopes -2.96 (theta=0) and -2.85, -2.71, -2.65 (also in the alt text)",
           ", ".join("%.4f" % s for s in sL), "each rounds to the stated value",
           "PASS" if all(rounds_to(s, t) for s, t in zip(sL, stated)) else "FAIL")
    record("6.3-8", S, "L2 rates stay between N^-m (= -3) and N^{1-m} (= -2)",
           "range [%.3f, %.3f]" % (min(sL), max(sL)), "-3 < slope < -2",
           "PASS" if all(-3 < s < -2 for s in sL) else "FAIL")
    C, why = load("e4_sobolev_coarse.npz", "python e4_sobolev.py 48 --coarse")
    if C is None:
        skip("6.3-9", S, "A coarser rule changes the reported errors by less than 0.1%", why)
    elif list(C["Ns"]) != list(Ns) or list(C["thetas"]) != list(D["thetas"]):
        skip("6.3-9", S, "A coarser rule changes the reported errors by less than 0.1%",
             "e4_sobolev_coarse.npz uses N = %s; rerun: python e4_sobolev.py 48 --coarse" % list(map(int, C["Ns"])))
    else:
        rel = max(float(np.max(np.abs(C["E"] - E) / E)), float(np.max(np.abs(C["exc"][1:] - exc[1:]) / exc[1:])),
                  float(np.max(np.abs(C["L2"] - L2) / L2)))
        record("6.3-9", S, "A coarser rule (Gauss in |x|^2, fewer angular nodes) changes the reported errors by less than 0.1%",
               "max relative change %.2e (gradient error, excess error for theta > 0, L2 error)" % rel,
               "< 1e-3", "PASS" if rel < 1e-3 else "FAIL")


# ===================================================================== report
def write_reports():
    counts = {k: sum(r["status"] == k for r in RESULTS) for k in ORDER}
    with open(os.path.join(HERE, "verification_report.json"), "w", encoding="utf-8") as f:
        json.dump(dict(summary=counts, checks=RESULTS), f, indent=1)
    with open(os.path.join(HERE, "verification_report.md"), "w", encoding="utf-8") as f:
        f.write("# Verification of the numbers in Section 6\n\n")
        f.write("Generated %s. Summary: %s.\n\n" % (time.strftime("%Y-%m-%d %H:%M"),
                                                    ", ".join("%d %s" % (counts[k], k) for k in ORDER)))
        f.write("PASS: holds as written. PASS*: holds after rounding the computed value to the "
                "precision used in the paper. FAIL: does not hold. SKIP: not evaluated.\n\n")
        f.write("| Check | Status | Statement in the paper | Computed | Criterion | Note |\n|---|---|---|---|---|---|\n")
        for r in sorted(RESULTS, key=lambda r: ORDER[r["status"]]):
            cells = [r["id"], "**%s**" % r["status"], r["claim"], r["computed"], r["criterion"], r["note"]]
            f.write("| " + " | ".join(c.replace("|", "\\|") for c in cells) + " |\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="skip slow recomputations")
    ap.add_argument("--section", choices=["6.0", "6.1", "6.2", "6.3"], action="append")
    ap.add_argument("--no-files", action="store_true", help="do not write report files")
    a = ap.parse_args(argv)
    secs = a.section or ["6.0", "6.1", "6.2", "6.3"]
    RESULTS.clear()
    t0 = time.time()
    print("Verifying Section 6 numbers (%s)%s" % (", ".join(secs), " [quick]" if a.quick else ""))
    for s, fn in (("6.0", section_60), ("6.1", section_61), ("6.2", section_62), ("6.3", section_63)):
        if s in secs:
            print("\n--- Section %s " % s + "-" * 60)
            try:
                fn(a.quick)
            except Exception as e:          # report, do not crash the calling script
                record(s + "-ERR", s, "verification of Section %s could not complete" % s,
                       "%s: %s" % (type(e).__name__, e), "-", "FAIL",
                       "check that the data files come from the current scripts")
    counts = {k: sum(r["status"] == k for r in RESULTS) for k in ORDER}
    print("\nSummary: " + ", ".join("%d %s" % (counts[k], k) for k in ORDER) + "  (%.0f s)" % (time.time() - t0))
    for r in RESULTS:
        if r["status"] in ("FAIL", "PASS*"):
            print(_ascii("  %-5s %-7s %s" % (r["status"], r["id"], r["claim"])))
    if not a.no_files:
        write_reports()
        print("Reports written: verification_report.md, verification_report.json")
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
