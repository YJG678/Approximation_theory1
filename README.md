# Numerical experiments for Section 6

Requires Python 3 with numpy, scipy and matplotlib, plus a LaTeX installation
(the figures render labels with usetex, which needs cm-super and dvipng).
The exact-arithmetic check additionally needs sympy and mpmath
(`pip install sympy mpmath`).

## Producing the data and figures

| Step | Command | Output |
|---|---|---|
| Forward stability (Sec. 6.1) | `python run_forward_stability.py`, then `python plot_forward_stability.py` | `forward_stability_data.npz`, `fig_forward_stability.pdf` |
| Truncation, K_ref = 300 and 80 (Sec. 6.2) | `python kref_check.py` (or `plus` / `circle` separately) | `kref_circle_K300.npz`, `kref_plus_K80.npz` |
| Truncation, analytic density and w_+ with K_ref = 56 | `python e2_nonradial.py` (default 56) | `e2_nonradial_K56.npz` |
| Truncation figure | `python plot_truncation.py` | `fig_truncation.pdf`, `fig_truncation_data.npz`, `truncation_slopes.npy` |
| Sobolev rates (Sec. 6.3) | `python e4_sobolev.py 48`, then `python plot_sobolev.py` | `e4_sobolev.npz`, `fig_sobolev.pdf` |
| Sobolev robustness rerun (Sec. 6.3) | `python e4_sobolev.py 48 --coarse` | `e4_sobolev_coarse.npz` |
| Truncation with K_ref = 200 for w_o (cross-check only) | `python e2_radial_kink.py` | `e2_radial_kink.npy` |
| Inverse-modulus probe (not in the paper) | `python e3_run.py` | `e3_radial.npy` |

## Verifying every number in Section 6

```
python verify_numerics.py            # all 45 checks, incl. independent recomputations (2-6 min)
python verify_numerics.py --quick    # stored data only, seconds
python verify_numerics.py --section 6.2
```

Each number or numerical statement in Section 6 is one check. The console
lists each check with the computed value; `verification_report.md` and
`verification_report.json` hold the same table. The exit code is 1 if any
check fails.

- **PASS**: the statement holds for the computed values exactly as written.
- **PASS\***: the statement holds only after rounding the computed value to
  the number of decimals the paper uses (for example 0.2525 reported as "0.25").
- **FAIL**: the statement does not hold.
- **SKIP**: a data file or package is missing; the message names the script to run.

The checks use three kinds of evidence:
- the stored data files written by the scripts above;
- independent recomputations: an exact-rational value (`exact_defect.py`),
  quadrature refinements, a second quadrature to estimate the round-off floor,
  and a symmetry-split versus full-basis comparison;
- consistency checks, such as the figure using the same fitted exponents as
  the text.

`run_forward_stability.py`, `plot_truncation.py` and `e4_sobolev.py` also run
the stored-data checks of their own subsection at the end.

## Files

- `diskpoly.py`: the solver (generalized Zernike basis, quadrature, QR-based
  projections, nested computation of delta_{N,K}).
- `exact_defect.py`: exact rational computation of delta_{3,5}(rho_*)^2.
- `plotstyle.py`: shared figure style.
- `validate2.py`, `d4check.py`: the earlier stand-alone solver checks, now
  superseded by `verify_numerics.py`.
- `fig_projection_geometry.tex`: the TikZ figure that was drafted and then
  left out of the paper.
