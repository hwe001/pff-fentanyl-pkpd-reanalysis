import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Out-of-sample and joint test using the dose-ranging study (Canfield, Fry, Sprague 2025): 100 ug/kg singles with their own plasma PK and brain/plasma ratios."""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pk_stage as pk, pd_stage as P

R = RES
T, DT = P.T, P.DT
# ---------------- 100 ug/kg inputs (digitised; values marked * were read by eye where a marker was hidden) ----------------
W100 = 0.3465                                         # kg, mean of the 314-379 g range of the cannulated rats
DOSE100 = 100 * W100 * 1000                            # ng
NCA100 = {'fentanyl': dict(AUC=47.46, t_half=1.39), 'pFF': dict(AUC=33.55, t_half=1.51)}
PL100 = {'fentanyl': [34.08, 19.6, 6.0, 1.4, 0.5], 'pFF': [25.60, 13.8, 4.0, 1.2, 0.4]}      # 30,60,120*,240*,480* min; 30 min = published Cmax
BP100 = {'fentanyl': np.mean([21.65, 20.30, 27.66, 21.40]), 'pFF': np.mean([43.22, 43.35, 44.82, 39.67])}
CTRL_T = {30: 0.28, 60: -0.07, 120: -0.15, 240: -0.50, 480: -0.40}                        # saline group, Fig 3A (by eye)
# observations at 100 ug/kg: (kind, drug, t, y, sigma); T = temperature change, M = tail flick, Mc = ceiling constraint
OBS100 = [('T', 'fentanyl', 30, -0.94, 0.35), ('T', 'fentanyl', 60, -1.84, 0.7), ('T', 'fentanyl', 120, -1.00, 0.65), ('T', 'fentanyl', 240, -0.62, 0.4), ('T', 'fentanyl', 480, -0.62, 0.4),
          ('T', 'pFF', 30, -1.11, 0.35), ('T', 'pFF', 60, -1.30, 0.7), ('T', 'pFF', 120, -1.08, 0.8), ('T', 'pFF', 240, -0.17, 0.3), ('T', 'pFF', 480, 0.08, 0.3),
          ('Mc', 'fentanyl', 30, 100, 8), ('Mc', 'fentanyl', 60, 100, 8), ('M', 'fentanyl', 120, 58.15, 12), ('M', 'fentanyl', 240, 24.17, 12), ('M', 'fentanyl', 480, 25.0, 12),
          ('Mc', 'pFF', 30, 100, 8), ('M', 'pFF', 60, 86.28, 12), ('M', 'pFF', 120, 69.87, 14), ('M', 'pFF', 240, 14.5, 10), ('M', 'pFF', 480, 43.17, 11)]


def fit_pk100(drug):
    t = np.array([0, 30, 60, 120, 240, 480]) / 60.0; y = np.array([0] + PL100[drug]); CL = DOSE100 / NCA100[drug]['AUC']; beta = np.log(2) / NCA100[drug]['t_half']; sd = np.maximum(0.12 * y, 0.5)

    def model(par):
        ka, V1, Q = np.exp(par); V2 = pk.solve_V2(CL, V1, Q, beta); return None if not np.isfinite(V2) else (ka, V1, Q, V2)

    def res(par):
        m = model(par)
        if m is None: return np.full(5, 1e3)
        return (pk.conc(t[1:], m[0], CL, m[1], m[2], m[3], dose=DOSE100) - y[1:]) / sd[1:]
    best = None
    for x0 in ([2, 300, 600], [4, 500, 1000], [1.5, 200, 400], [3, 800, 1500]):
        try: r = least_squares(res, np.log(x0), bounds=(np.log([0.3, 30, 20]), np.log([30, 8000, 30000])))
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    ka, V1, Q, V2 = model(best.x); return dict(ka=ka, CL=CL, V1=V1, Q=Q, V2=V2, chi2=2 * best.cost)


pk100 = {d: fit_pk100(d) for d in ('fentanyl', 'pFF')}
CP100 = {d: pk.conc(T / 60.0, p['ka'], p['CL'], p['V1'], p['Q'], p['V2'], dose=DOSE100) for d, p in pk100.items()}
CP300 = P.CP; KP300 = {'fentanyl': P.kp_scenarios()['S1 brain/plasma ratio, mean of 4 regions'][0], 'pFF': P.kp_scenarios()['S1 brain/plasma ratio, mean of 4 regions'][1]}
GROUPS = {'f300': (CP300['fentanyl'], KP300['fentanyl'], 'fentanyl'), 'p300': (CP300['pFF'], KP300['pFF'], 'pFF'), 'f100': (CP100['fentanyl'], BP100['fentanyl'], 'fentanyl'), 'p100': (CP100['pFF'], BP100['pFF'], 'pFF')}


def effect(par, grp, rho_free, kind, t):
    lE, lEC, lk, lECA, lkA, lrho = par; Emax, EC, ke, ECA, keA = np.exp([lE, lEC, lk, lECA, lkA]); rho = np.exp(lrho) if rho_free else 1.0
    cp, kp, drug = GROUPS[grp]; cb = kp * cp; sc = rho if drug == 'pFF' else 1.0
    if kind == 'T': ce = np.interp(t, T, P.filt(cb, ke)); return -Emax * ce / (EC * sc + ce)
    ce = np.interp(t, T, P.filt(cb, keA)); return 100 * ce / (ECA * sc + ce)


def obs_set(use300=True, use100=True, ctrl=False):
    o = []
    if use300:
        for k, g, t, y, s in P.OBS: o.append((k, 'f300' if g == 'fentanyl' else 'p300', t, y, s))
    if use100:
        for k, g, t, y, s in OBS100: o.append((k, 'f100' if g == 'fentanyl' else 'p100', t, y - (CTRL_T[t] if (ctrl and k == 'T') else 0.0), s))
    return o


def resid(par, obs, rho_free):
    r = []
    for k, g, t, y, s in obs:
        p = effect(par, g, rho_free, 'T' if k == 'T' else 'M', t); r.append(max(0.0, 100 - p) / s if k == 'Mc' else (p - y) / s)
    return np.array(r)


def fit(obs, rho_free, nstart=30, seed=0, fixed_rho=None):
    rng = np.random.default_rng(seed); best = None
    for i in range(nstart):
        x0 = rng.uniform(P.LO, P.HI) if i else np.log([6.0, 500, 0.03, 20, 0.07, 1.0]); lo, hi = P.LO.copy(), P.HI.copy()
        if fixed_rho is not None: x0[5] = np.log(fixed_rho); lo[5] = x0[5] - 1e-9; hi[5] = x0[5] + 1e-9
        elif not rho_free: x0[5] = 0.0; lo[5] = -1e-9; hi[5] = 1e-9
        try: r = least_squares(resid, x0, args=(obs, rho_free), bounds=(lo, hi), x_scale='jac')
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best


if __name__ == '__main__':
    print('100 ug/kg PK fits:', {d: {k: round(float(v), 2) for k, v in p.items()} for d, p in pk100.items()})
    print('K_p at 100: fentanyl %.1f, pFF %.1f (ratio %.2f) | at 300: %.1f, %.1f (ratio %.2f)' % (BP100['fentanyl'], BP100['pFF'], BP100['pFF'] / BP100['fentanyl'], KP300['fentanyl'], KP300['pFF'], KP300['pFF'] / KP300['fentanyl']))
    prim = json.load(open(R + os.sep + 'primary_fit.json')); xs, xp = np.array(prim['x_exposure']), np.array(prim['x_potency'])
    res = {}
    for ctrl in (False, True):
        o100 = obs_set(False, True, ctrl)
        for name, x, rf in (('shared PD (fitted at 300)', xs, False), ('pFF EC50 free, rho=%.2f (fitted at 300)' % np.exp(xp[5]), xp, True)):
            r = resid(x, o100, rf); chi = float((r ** 2).sum())
            res[(ctrl, name)] = chi; print(f'OUT-OF-SAMPLE at 100 ug/kg | control-corrected={ctrl} | {name}: chi2 = {chi:.1f} on {len(o100)} observations')
    # joint fit
    rows = []
    for ctrl in (False, True):
        ob = obs_set(True, True, ctrl); N = len(ob)
        f0 = fit(ob, False); f1 = fit(ob, True); c0, c1 = 2 * f0.cost, 2 * f1.cost
        grid = np.exp(np.linspace(np.log(0.3), np.log(4), 33)); prof = np.array([2 * fit(ob, True, nstart=6, fixed_rho=g).cost for g in grid]); ok = grid[prof - prof.min() < 3.84]
        rows.append(dict(control_corrected=ctrl, n_obs=N, chi2_shared=round(c0, 1), chi2_ecc50_free=round(c1, 1), AICc_shared=round(P.aicc(c0, 5, N), 1), AICc_EC50free=round(P.aicc(c1, 6, N), 1), rho_hat=round(float(np.exp(f1.x[5])), 2), rho_lo=round(float(ok.min()), 2), rho_hi=round(float(ok.max()), 2),
                         Emax_T_shared=round(float(np.exp(f0.x[0])), 2)))
        print(rows[-1], flush=True)
        json.dump(dict(x_shared=list(f0.x), x_potency=list(f1.x), grid=list(grid), profile=list(prof)), open(R + f'\\dose_test_joint_ctrl{int(ctrl)}.json', 'w'))
    pd.DataFrame(rows).to_csv(R + os.sep + 'dose_test_joint.csv', index=False)
    # per-observation predictions at 100 ug/kg from the joint EC50-free and shared fits (raw data)
    j = json.load(open(R + os.sep + 'dose_test_joint_ctrl0.json')); o100 = obs_set(False, True, False); out = []
    for k, g, t, y, s in o100:
        out.append(dict(group=g, endpoint='temperature' if k == 'T' else 'tail flick', t_min=t, observed=y, pred_shared_joint=round(float(effect(np.array(j['x_shared']), g, False, 'T' if k == 'T' else 'M', t)), 2), pred_EC50free_joint=round(float(effect(np.array(j['x_potency']), g, True, 'T' if k == 'T' else 'M', t)), 2),
                        pred_shared_300only=round(float(effect(xs, g, False, 'T' if k == 'T' else 'M', t)), 2), pred_EC50free_300only=round(float(effect(xp, g, True, 'T' if k == 'T' else 'M', t)), 2)))
    pd.DataFrame(out).to_csv(R + os.sep + 'dose_test_predictions_100.csv', index=False); print(pd.DataFrame(out).to_string(index=False))
