# Verification of the numbers in Section 6

Generated 2026-10-01 12:11. Summary: 3 FAIL, 6 PASS*, 0 SKIP, 36 PASS.

PASS: holds as written. PASS*: holds after rounding the computed value to the precision used in the paper. FAIL: does not hold. SKIP: not evaluated.

| Check | Status | Statement in the paper | Computed | Criterion | Note |
|---|---|---|---|---|---|
| 6.2-13 | **FAIL** | ... and reaches round-off at K = 22 | remainder at K=22: 2.17e-15; numerical floor 8.3e-18 (N=2 vs N=3 spread for K >= 7 and a rerun with another quadrature); floor first reached at K = 27 | remainder(22) <= 10 x floor | the value at K = 22 is 260 times the floor and reproducible; the geometric decay continues to about K = 27 |
| 6.2-15 | **FAIL** | Caption of Fig. 2(b): the case N = 3 coincides (plotted range K = 3..26) | N=2 and N=3 differ at K = [4, 6] (e.g. K=4: 9.432e-04 vs 5.264e-04; K=6: 2.488e-04 vs 1.500e-04); identical for K >= 7 | relative difference <= 1e-6 or below the floor at every plotted K |  |
| 6.2-16 | **FAIL** | Refining the quadrature from 40 to 72 nodes per direction changes delta_{2,24}(w_+ nu_0)^2 by less than 1e-15 | 1.888866391130369e-03 vs 1.888866391129168e-03: change 1.20e-15 | < 1e-15 | the change is 1.2e-15, i.e. about 1e-15 (relative 6e-13) |
| 6.1-4 | **PASS*** | delta_{N,24} lies between 0.011 and 0.25 times the bound theta/(1-theta) | ratio min 0.01134, max 0.25245 | 0.011 <= ratio <= 0.25 | max ratio 0.2525 exceeds 0.25 but rounds to it at two decimals |
| 6.2-8 | **PASS*** | Decreasing K_ref (200 for w_o, 56 for w_+) or fitting to K_ref/2 changes the exponents by at most 0.09 | max change 0.0922; w_o, N=3: -3.0079 \| Kref/2 -3.0434 \| Kref=200 -3.0063 \| Kref=200,/2 -3.0454 ; w_+, N=2: -2.9189 \| Kref/2 -3.0112 \| Kref=56 -2.8692 \| Kref=56,/2 -3.0103 ; w_+, N=3: -3.0167 \| Kref/2 -3.0705 \| Kref=56 -3.0489 \| Kref=56,/2 -3.1024 | <= 0.09 | max change 0.0922 is above 0.09 but rounds to it |
| 6.2-9 | **PASS*** | Every fitted exponent lies between -3.10 and -2.87 | range [-3.1024, -2.8692] over 12 fits | -3.10 <= exponent <= -2.87 | extremes -3.1024 and -2.8692 lie outside but round to the endpoints |
| 6.3-3 | **PASS*** | Gradient error slopes between -2.01 and -2.05 for all four theta (N >= 12) | -2.0113, -2.0145, -2.0252, -2.0508 | -2.05 <= slope <= -2.01 | extreme -2.0508 lies outside but rounds to the endpoint |
| 6.3-4 | **PASS*** | At N = 48 the gradient error lies between 0.78 and 1 times its theta = 0 value | 1.0000, 0.9801, 0.9143, 0.7763 | 0.78 <= ratio <= 1 | smallest ratio 0.7763 is below 0.78 but rounds to it |
| 6.3-6 | **PASS*** | Excess error slopes from -4.77 to -4.86, steeper than the limiting order N^{2-2m} = N^-4 | -4.7661, -4.7961, -4.8638 | -4.86 <= slope <= -4.77 and < -4 | extreme -4.7661 lies outside but rounds to the endpoint |
| 6.0-1 | **PASS** | The solver reproduces (4.22) [\|\|C_3 u\|\|_nu = sqrt(6)/70 \|eps\|] to at least twelve digits | max relative error 9.5e-14 over eps in {0.5,-0.3,0.9} (13.0 digits) | relative error <= 1e-12 |  |
| 6.0-2 | **PASS** | Exact rational value delta_{3,5}(rho_*)^2 = 1.2152578150... x 10^-3 | exact (rational matrices, 40-digit eigenvalue) = 0.0012152578150466218048 | leading digits 1.2152578150 |  |
| 6.0-3 | **PASS** | Solver agrees with the exact value to at least twelve digits | solver 1.215257815046725e-03, relative error 8.5e-14 (13.1 digits) | relative error <= 1e-12 |  |
| 6.0-4 | **PASS** | Radial mu: restricting inputs to angular frequency <= 3 does not change delta_{3,K} | rho_*: max relative difference 5.9e-16 (K = 8, 16, 24); max \|\|C_3 u\|\|/\|\|grad u\|\| over 242 inputs with m >= 4: 1.1e-15 | relative difference <= 1e-12 and ratio <= 1e-12 |  |
| 6.0-5 | **PASS** | Radial reduction allows input degrees up to K = 300 | largest input degree in kref_circle_K300.npz: 300 | = 300 |  |
| 6.0-6 | **PASS** | Square symmetry splits the inputs into five classes and allows K = 80 | largest input degree 80; classes A1/A2/B1/B2/E; misassigned basis functions: 0 | K = 80, five classes, no misassignment |  |
| 6.0-7 | **PASS** | Five-class split reproduces the full defect (w_+, N = 2, 3, K <= 32) | max relative difference 2.7e-15 | <= 1e-12 |  |
| 6.0-8 | **PASS** | delta_{N,60}(nu_alpha) <= 2 x 10^-13 for alpha in {-1/2,0,2}, N in {2,3,6} | max 1.66e-13; a=-0.5,N=2:8.6e-15, a=-0.5,N=3:1.0e-13, a=-0.5,N=6:1.7e-13, a=0,N=2:7.7e-15, a=0,N=3:2.3e-14, a=0,N=6:3.4e-14, a=2,N=2:6.3e-15, a=2,N=3:2.3e-14, a=2,N=6:3.3e-14 | all <= 2e-13 | round-off level; the margin to 2e-13 depends on the BLAS/LAPACK build |
| 6.1-1 | **PASS** | Setting: nu = nu_0, K = 24, N = 2..10 (panels a,c,d), theta = 0.1 (panel b) | alpha = 0, K = 24, N = 2..10, theta_b = 0.1, l = [5, 9, 13, 17] | as stated |  |
| 6.1-2 | **PASS** | h is centred and scaled: nu_0(h) = 0 and \|\|h\|\|_inf = 1 | nu_0(h) = 0.0e+00; \|\|h\|\|_inf = 1 +4.6e-07 (grid normalisation vs refined maximisation) | \|nu_0(h)\| <= 1e-14 and \| \|\|h\|\|_inf - 1 \| <= 1e-6 |  |
| 6.1-3 | **PASS** | delta_{N,24}(mu_eps) has fitted slope 1.000 in theta for theta <= 1e-2, every N = 2..10 | slopes 0.99976, 0.99967, 0.99973, 0.99969, 0.99965, 0.99963, 0.99963, 0.99963, 0.99963 | each rounds to 1.000 |  |
| 6.1-5 | **PASS** | Panel (b): each ratio delta_{N,24}/theta rises to about 0.38 | maxima l=5: 0.3763, l=9: 0.3839, l=13: 0.3839, l=17: 0.3839; first values 0.307, 0.252, 0.218, 0.195 | each maximum rounds to 0.38 and exceeds the first value |  |
| 6.1-6 | **PASS** | Panel (b): stays near the plateau until N ~ l - 1, then decreases | l=5: plateau ends at N=4, then decreasing; l=9: plateau ends at N=8, then decreasing; l=13: plateau ends at N=12, then decreasing; l=17: plateau ends at N=16 (plateau = within 3% of the maximum) | plateau end within 1 of l-1; strictly decreasing afterwards | for l = 17 the decrease lies beyond the tested range N <= 16 |
| 6.1-7 | **PASS** | Panel (b): all ratios remain below the degree-independent envelope 1/(1-theta) | max 0.3839 vs 1/(1-0.1) = 1.1111 | max <= 1/(1-theta) |  |
| 6.1-8 | **PASS** | Excess squared error and linearization remainder both have fitted slope 2.00 | excess 1.9979..1.9989, remainder 1.9986..1.9999 (theta <= 1e-2) | each rounds to 2.00 |  |
| 6.1-9 | **PASS** | Direct difference of best-approximation errors agrees to a relative 7e-6 (theta >= 1e-2) | 1.76e-06 (stored check: 1.76e-06) | <= 7e-6 | cancellation-limited; the value varies with the BLAS build (observed 1.8e-6 to 6.5e-6) |
| 6.1-10 | **PASS** | All finite-input defects and fixed-input errors remain below their envelopes | max data/bound: (a) 0.252, (b) 0.346, (c) 1.58e-03, (d) 1.74e-02 | all <= 1 |  |
| 6.1-11 | **PASS** | ... in some panels by several orders of magnitude | min data/bound: (a) 4.8e-02 = 1.3 orders, (c) 2.4e-05 = 4.6 orders, (d) 9.1e-04 = 3.0 orders | some panel >= 2 orders below |  |
| 6.2-1 | **PASS** | The three measures are probability measures, centred, and isotropic after scalar rescaling | max \|mass-1\|, \|mean\|, \|m_xy\|, \|m_xx-m_yy\|: w_o 1e-16, w_+ 1e-15, (5.20) 3e-17 | all <= 1e-14 |  |
| 6.2-2 | **PASS** | Their densities lie between 0.57 and 1.29 | exact range [0.57559, 1.28269]; sampled range [0.57559, 1.28269] | 0.57 <= w <= 1.29 |  |
| 6.2-3 | **PASS** | w_o nu_0 is rotationally invariant, so delta_2 = 0 | delta_{2,K} = 1.4e-15, 1.4e-15, 2.8e-15 (K = 8, 16, 24) | <= 1e-12 |  |
| 6.2-4 | **PASS** | K_ref = 300 for w_o, 80 for w_+, 56 for the analytic density | K_ref in data: 300, 80, 56 | = 300, 80, 56 |  |
| 6.2-5 | **PASS** | Remainder plotted only for K <= K_ref/2 (panel a) | largest plotted K: w_o 150, w_+ N=2 40, w_+ N=3 40 | <= 150, 40, 40 |  |
| 6.2-6 | **PASS** | Exponents fitted over 10 <= K <= K_ref/3: -3.01 (w_o), -2.92 (w_+, N=2), -3.02 (w_+, N=3) | w_o, N=3 -3.0079, w_+, N=2 -2.9189, w_+, N=3 -3.0167 | each rounds to the stated value |  |
| 6.2-7 | **PASS** | The exponents used for the figure (plot_truncation.py) are the same numbers | max \|difference\| 0.0e+00 | <= 1e-10 |  |
| 6.2-10 | **PASS** | The remainder decays like K^-3 rather than K^-2 | all 12 exponents closer to -3 than to -2 | \|e+3\| < \|e+2\| |  |
| 6.2-11 | **PASS** | Analytic density: the remainder decays geometrically | log10-linear fit K = 5..21: -0.57 decades per degree, correlation -0.976, nonincreasing: True | correlation <= -0.95 and nonincreasing |  |
| 6.2-12 | **PASS** | ... falls below 1e-12 at K = 20 | first K with remainder < 1e-12: 20 (K=19: 1.31e-12, K=20: 4.11e-13) | = 20 |  |
| 6.2-14 | **PASS** | ... the two levels (N = 2, 3) give identical values there (K ~ 20-22) | max \|difference\| 2.2e-19 | <= numerical floor 8.3e-18 |  |
| 6.2-17 | **PASS** | w_+ runs with K_ref = 56 and 80 (different quadratures) agree to a relative 1e-12 for 5 <= K <= 56 | max relative difference 7.1e-13 | <= 1e-12 |  |
| 6.3-1 | **PASS** | Setting: h = sign(x1 x2) with nu_0(h) = 0, \|\|h\|\|_inf = 1; theta in {0,1/4,1/2,3/4}; N up to 48 | nu_0(h) = -5.9e-18; theta = [0.0, 0.25, 0.5, 0.75]; N = [4, 5, 6, 7, 9, 10, 13, 15, 18, 22, 27, 33, 40, 48] | as stated |  |
| 6.3-2 | **PASS** | nu_0 best-approximation errors decay at the H^{3-} rates (N^-2 gradient, N^-3 L2) | theta=0: gradient -2.011, L2 -2.956 | within 0.05 of -2 and -3 |  |
| 6.3-5 | **PASS** | Excess error ratios at N = 48 (theta = 1/2, 3/4 vs 1/4): 3.82 and 7.57, vs 4 and 9 | 3.8170, 7.5718; exact quadratic scaling (0.5/0.25)^2 = 4, (0.75/0.25)^2 = 9 | round to 3.82, 7.57 |  |
| 6.3-7 | **PASS** | L2 error slopes -2.96 (theta=0) and -2.85, -2.71, -2.65 (also in the alt text) | -2.9560, -2.8467, -2.7092, -2.6490 | each rounds to the stated value |  |
| 6.3-8 | **PASS** | L2 rates stay between N^-m (= -3) and N^{1-m} (= -2) | range [-2.956, -2.649] | -3 < slope < -2 |  |
| 6.3-9 | **PASS** | A coarser rule (Gauss in \|x\|^2, fewer angular nodes) changes the reported errors by less than 0.1% | max relative change 7.09e-04 (gradient error, excess error for theta > 0, L2 error) | < 1e-3 |  |
