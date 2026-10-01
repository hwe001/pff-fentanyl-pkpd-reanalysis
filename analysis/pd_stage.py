import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Stage 2: PD. Brain exposure -> effect compartment -> hypothermia and tail flick.
H_exposure: both drugs share every PD parameter; they differ only through their plasma curve and brain/plasma ratio.
H_potency : pFF has its own EC50 (a single fold-change rho on both endpoints, as a receptor-affinity difference would give)."""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.signal import lfilter
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pk_stage as pk

D = DATA
pkfit = json.load(open(os.path.join(RES, 'pk_fit.json')))
TMAX, DT = 480.0, 0.25
T = np.arange(0, TMAX + DT, DT)                        # minutes
CP = {d: pk.conc(T / 60.0, pkfit[d]['ka_per_h'], pkfit[d]['CL_F_mL_per_h'], pkfit[d]['V1_F_mL'], pkfit[d]['Q_F_mL_per_h'], pkfit[d]['V2_F_mL']) for d in ('fentanyl', 'pFF')}

temp = pd.read_csv(D + os.sep + 'fig4A_temperature_change.csv'); mpe = pd.read_csv(D + os.sep + 'fig3_tail_flick_MPE.csv')
bp = pd.read_csv(D + os.sep + 'fig5_brain_group_means.csv'); bp_ratio = bp[bp.panel.str.startswith('B')]; bp_conc = bp[bp.panel.str.startswith('A')]
REGIONS = ['Hippocampus', 'Medulla', 'Striatum', 'Frontal Cortex']


def kp_scenarios():
    r = {}
    for g in ('fentanyl', 'pFF'):
        r[g] = {reg: float(bp_ratio[(bp_ratio.group == g) & (bp_ratio.region == reg)]['mean'].iloc[0]) for reg in REGIONS}
    absb = {g: {reg: float(bp_conc[(bp_conc.group == g) & (bp_conc.region == reg)]['mean'].iloc[0]) for reg in REGIONS} for g in ('fentanyl', 'pFF')}
    cp60 = {g: float(np.interp(60, T, CP[g])) for g in CP}
    sc = {'S1 brain/plasma ratio, mean of 4 regions': (np.mean(list(r['fentanyl'].values())), np.mean(list(r['pFF'].values())), None),
          'S2 absolute brain / PK-cohort plasma at 60 min': (np.mean(list(absb['fentanyl'].values())) / cp60['fentanyl'], np.mean(list(absb['pFF'].values())) / cp60['pFF'], None),
          'S3 no brain difference (pFF given the fentanyl ratio)': (np.mean(list(r['fentanyl'].values())),) * 2 + (None,),
          'S4 as S1, finite brain equilibration (t1/2 10 min)': (np.mean(list(r['fentanyl'].values())), np.mean(list(r['pFF'].values())), 10.0),
          'S5 as S1, finite brain equilibration (t1/2 30 min)': (np.mean(list(r['fentanyl'].values())), np.mean(list(r['pFF'].values())), 30.0)}
    for reg in REGIONS: sc[f'S6 {reg} only'] = (r['fentanyl'][reg], r['pFF'][reg], None)
    return sc


def filt(x, k):
    a = np.exp(-k * DT); return lfilter([(1 - a) / 2, (1 - a) / 2], [1, -a], x)


def brain(drug, kp, tlag):
    cb = kp * CP[drug]
    return cb if tlag is None else filt(cb, np.log(2) / tlag)         # brain concentration (ng/g)


def sigma_T(row): return max(row.sem_estimate if row.sem_estimate == row.sem_estimate else 0.3, 0.25) * (1.6 if 'overlap' in str(row.note) else 1.0)


def build_data():
    obs = []
    for g in ('fentanyl', 'pFF'):
        for r in temp[(temp.group == g) & (temp.t_min > 0)].itertuples():
            obs.append(('T', g, r.t_min, r.delta_T_C, sigma_T(r)))
        for r in mpe[(mpe.group == g) & (mpe.t_min > 0)].itertuples():
            ceil = r.t_min <= 120
            obs.append(('Mc' if ceil else 'M', g, r.t_min, 100.0 if ceil else r.MPE_percent, 8.0 if ceil else max(r.sem_estimate, 8.0)))
    return obs


OBS = build_data()


def predict(par, kps, tlags, rho_free):
    lE, lEC, lk, lECA, lkA, lrho = par; Emax, EC, ke, ECA, keA = np.exp([lE, lEC, lk, lECA, lkA]); rho = np.exp(lrho) if rho_free else 1.0
    out = []
    for kind, g, t, y, s in OBS:
        cb = brain(g, kps[g], tlags); sc = rho if g == 'pFF' else 1.0
        if kind == 'T': ce = filt_ce(cb, ke); out.append(-Emax * np.interp(t, T, ce) / (EC * sc + np.interp(t, T, ce)))
        else: ce = filt_ce(cb, keA); c = np.interp(t, T, ce); out.append(100 * c / (ECA * sc + c))
    return np.array(out)


def filt_ce(cb, k): return filt(cb, k)


def resid(par, kps, tlags, rho_free, data=None):
    pr = predict(par, kps, tlags, rho_free); r = []
    for (kind, g, t, y, s), p in zip(data or OBS, pr):
        if kind == 'Mc': r.append(max(0.0, 100.0 - p) / s)
        else: r.append((p - y) / s)
    return np.array(r)


LO = np.log([2.0, 10, 0.002, 1, 0.002, 0.05]); HI = np.log([10.0, 1e5, 2.0, 1e4, 2.0, 20.0])


def fit(kps, tlags, rho_free, fixed_rho=None, nstart=24, seed=0, data=None):
    rng = np.random.default_rng(seed); best = None
    for i in range(nstart):
        x0 = rng.uniform(LO, HI) if i else np.log([6.0, 500, 0.03, 40, 0.1, 1.0])
        if fixed_rho is not None: x0[5] = np.log(fixed_rho)
        lo, hi = LO.copy(), HI.copy()
        if fixed_rho is not None: lo[5] = np.log(fixed_rho) - 1e-9; hi[5] = np.log(fixed_rho) + 1e-9
        elif not rho_free: lo[5] = -1e-9; hi[5] = 1e-9; x0[5] = 0.0
        try: r = least_squares(resid, x0, args=(kps, tlags, rho_free, data), bounds=(lo, hi), x_scale='jac')
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best


def aicc(chi2, k, n): return chi2 + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)


if __name__ == '__main__':
    rows = []; N = len(OBS)
    for name, (kf, kpf, tlag) in kp_scenarios().items():
        kps = {'fentanyl': kf, 'pFF': kpf}
        f0 = fit(kps, tlag, False); f1 = fit(kps, tlag, True)
        chi0, chi1 = 2 * f0.cost, 2 * f1.cost
        grid = np.exp(np.linspace(np.log(0.2), np.log(5), 41)); prof = []
        for rho in grid: prof.append(2 * fit(kps, tlag, True, fixed_rho=rho, nstart=8).cost)
        prof = np.array(prof); ok = grid[prof - prof.min() < 3.84]
        rows.append(dict(scenario=name, Kp_fentanyl=round(kf, 2), Kp_pFF=round(kpf, 2), chi2_exposure_only=round(chi0, 2), chi2_with_potency=round(chi1, 2),
                         AICc_exposure_only=round(aicc(chi0, 5, N), 2), AICc_with_potency=round(aicc(chi1, 6, N), 2), rho_hat=round(float(np.exp(f1.x[5])), 2),
                         rho_95CI_low=round(float(ok.min()), 2), rho_95CI_high=round(float(ok.max()), 2), n_obs=N,
                         Emax_T=round(float(np.exp(f0.x[0])), 2), EC50_T=round(float(np.exp(f0.x[1])), 0), ke0_T=round(float(np.exp(f0.x[2])), 4), EC50_A=round(float(np.exp(f0.x[3])), 1), ke0_A=round(float(np.exp(f0.x[4])), 3)))
        print(rows[-1], flush=True)
        if name.startswith('S1'): json.dump(dict(grid=list(grid), profile_chi2=list(prof), x_exposure=list(f0.x), x_potency=list(f1.x), kps=kps), open(os.path.join(RES, 'primary_fit.json'), 'w'))
    pd.DataFrame(rows).to_csv(os.path.join(RES, 'pd_scenarios.csv'), index=False)
