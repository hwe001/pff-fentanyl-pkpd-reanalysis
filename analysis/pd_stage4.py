import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Stage 4: profile-likelihood 95% CIs for every PD parameter of the primary model (S1, pFF EC50 free), and derived quantities."""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pd_stage as P
R = RES
prim = json.load(open(R + os.sep + 'primary_fit.json')); kps = prim['kps']; x_hat = np.array(prim['x_potency'])
names = ['Emax_T (degC)', 'EC50_T fentanyl (ng/g)', 'ke0_T (1/min)', 'EC50_A fentanyl (ng/g)', 'ke0_A (1/min)', 'rho = EC50 pFF / fentanyl']


def fit_fixed(idx, val, nstart=6, seed=3):
    rng = np.random.default_rng(seed); best = None
    for i in range(nstart):
        x0 = x_hat.copy() if i == 0 else np.clip(x_hat + rng.normal(0, 0.6, 6), P.LO, P.HI); x0[idx] = np.log(val)
        lo, hi = P.LO.copy(), P.HI.copy(); lo[idx] = np.log(val) - 1e-9; hi[idx] = np.log(val) + 1e-9
        try: r = least_squares(P.resid, x0, args=(kps, None, True, None), bounds=(lo, hi), x_scale='jac')
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best


chi_hat = 2 * least_squares(P.resid, x_hat, args=(kps, None, True, None), bounds=(P.LO, P.HI), x_scale='jac').cost
rows = []
for i, nm in enumerate(names):
    centre = np.exp(x_hat[i]); grid = centre * np.exp(np.linspace(-2.3, 2.3, 41)); grid = grid[(np.log(grid) > P.LO[i] + 1e-6) & (np.log(grid) < P.HI[i] - 1e-6)]
    prof = np.array([2 * fit_fixed(i, v).cost for v in grid]) - chi_hat; ok = grid[prof < 3.84]
    lo_open = bool(ok.min() <= grid.min() * 1.0001); hi_open = bool(ok.max() >= grid.max() * 0.9999)
    rows.append(dict(parameter=nm, estimate=round(float(centre), 4), ci_low=round(float(ok.min()), 4), ci_high=round(float(ok.max()), 4), lower_unbounded=lo_open, upper_unbounded=hi_open))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(R + os.sep + 'pd_parameter_profile_CIs.csv', index=False)

# derived quantities from the primary model
Emax, EC, ke, ECA, keA, rho = np.exp(x_hat); out = {}
for g in ('fentanyl', 'pFF'):
    cb = P.brain(g, kps[g], None); pf = g == 'pFF'
    ce = P.filt(cb, ke); dT = -Emax * ce / (EC * (rho if pf else 1) + ce); ceA = P.filt(cb, keA); mp = 100 * ceA / (ECA * (rho if pf else 1) + ceA)
    out[g] = dict(nadir_dT=round(float(dT.min()), 2), t_nadir_min=round(float(P.T[dT.argmin()]), 0), peak_brain_ng_g=round(float(cb.max()), 0), t_peak_brain_min=round(float(P.T[cb.argmax()]), 0),
                  t_MPE_below_50_min=round(float(P.T[np.argmax((P.T > 60) & (mp < 50))]), 0), MPE_at_240=round(float(np.interp(240, P.T, mp)), 1), MPE_at_480=round(float(np.interp(480, P.T, mp)), 1),
                  dT_at_120=round(float(np.interp(120, P.T, dT)), 2), dT_at_240=round(float(np.interp(240, P.T, dT)), 2))
out['EC50_T_pFF'] = round(float(EC * rho), 0); out['EC50_A_pFF'] = round(float(ECA * rho), 1); out['chi2'] = round(float(chi_hat), 2); out['n_obs'] = len(P.OBS)
json.dump(out, open(R + os.sep + 'primary_derived.json', 'w'), indent=1); print(json.dumps(out, indent=1))
