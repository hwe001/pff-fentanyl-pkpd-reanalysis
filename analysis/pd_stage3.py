import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Stage 3: where does the exposure-only model fail, and is a potency (EC50) difference the best way to repair it? Alternatives: Emax or ke0 differing."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pd_stage as P
from scipy.optimize import least_squares

fit1 = json.load(open(os.path.join(RES, 'primary_fit.json'))); KPS = fit1['kps']


def predict2(par, extra):
    """par: log [EmaxT, EC50T, ke0T, EC50A, ke0A]; extra: dict of pFF multipliers on EC50 (rho), Emax_T (eps), ke0 (kap)"""
    Emax, EC, ke, ECA, keA = np.exp(par); rho, eps, kap = extra
    out = []
    for kind, g, t, y, s in P.OBS:
        cb = P.brain(g, KPS[g], None); pf = g == 'pFF'
        if kind == 'T':
            ce = P.filt(cb, ke * (kap if pf else 1)); c = np.interp(t, P.T, ce); out.append(-Emax * (eps if pf else 1) * c / (EC * (rho if pf else 1) + c))
        else:
            ce = P.filt(cb, keA * (kap if pf else 1)); c = np.interp(t, P.T, ce); out.append(100 * c / (ECA * (rho if pf else 1) + c))
    return np.array(out)


def resid2(x, which):
    par = x[:5]; ex = [1.0, 1.0, 1.0]
    for i, w in enumerate(which): ex[w] = np.exp(x[5 + i])
    pr = predict2(par, ex); r = []
    for (kind, g, t, y, s), p in zip(P.OBS, pr): r.append(max(0.0, 100 - p) / s if kind == 'Mc' else (p - y) / s)
    return np.array(r)


def fit2(which, nstart=40, seed=1):
    rng = np.random.default_rng(seed); lo = list(P.LO[:5]) + [np.log(0.05)] * len(which); hi = list(P.HI[:5]) + [np.log(20)] * len(which); best = None
    for i in range(nstart):
        x0 = rng.uniform(lo, hi) if i else np.log([6, 500, 0.03, 20, 0.07] + [1.0] * len(which))
        try: r = least_squares(resid2, x0, args=(which,), bounds=(lo, hi), x_scale='jac')
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best


N = len(P.OBS); rows = []
variants = {'exposure only': [], 'EC50 differs (rho)': [0], 'Emax_T differs': [1], 'ke0 differs': [2], 'EC50 + ke0 differ': [0, 2]}
res = {}
for name, w in variants.items():
    f = fit2(w); k = 5 + len(w); chi = 2 * f.cost; res[name] = f
    mult = {['rho', 'eps', 'kap'][j]: round(float(np.exp(f.x[5 + i])), 2) for i, j in enumerate(w)}
    rows.append(dict(model=name, n_params=k, chi2=round(chi, 2), AICc=round(P.aicc(chi, k, N), 2), pFF_multipliers=str(mult)))
    print(rows[-1], flush=True)
dfv = pd.DataFrame(rows); dfv['dAICc'] = (dfv.AICc - dfv.AICc.min()).round(2); dfv.to_csv(os.path.join(RES, 'pd_alternatives_S1.csv'), index=False)

# residuals, exposure-only vs potency model
rr = []
for name in ('exposure only', 'EC50 differs (rho)'):
    f = res[name]; w = variants[name]; ex = [1.0, 1.0, 1.0]
    for i, j in enumerate(w): ex[j] = np.exp(f.x[5 + i])
    pr = predict2(f.x[:5], ex)
    for (kind, g, t, y, s), p in zip(P.OBS, pr): rr.append(dict(model=name, endpoint={'T': 'temperature', 'M': 'tail flick', 'Mc': 'tail flick (ceiling)'}[kind], drug=g, t_min=t, observed=y, predicted=round(float(p), 2), z=round(float((p - y) / s) if kind != 'Mc' else max(0, 100 - p) / s, 2)))
pd.DataFrame(rr).to_csv(os.path.join(RES, 'pd_residuals_S1.csv'), index=False)
bad = pd.DataFrame(rr); bad = bad[(bad.model == 'exposure only') & (bad.z.abs() > 2)]
print('\nlargest misfits of the exposure-only model (|z| > 2):'); print(bad.to_string(index=False))
