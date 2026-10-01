import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Stage 1: plasma PK (2-compartment, first-order absorption from the sc depot), anchored on the source's noncompartmental values."""
import numpy as np, pandas as pd, json
from scipy.optimize import least_squares, brentq

D = DATA
DOSE_UG = 300 * 0.33          # 300 ug/kg x ~0.33 kg rat (group weights 327-330 g) = 99 ug
DOSE_NG = DOSE_UG * 1000
plasma = pd.read_csv(D + os.sep + 'fig2_plasma_concentration.csv')
nca = pd.read_csv(D + os.sep + 'table1_noncompartmental_PK_exact_from_source.csv').set_index('parameter')
DRUGS = {'fentanyl': 'fentanyl_mean', 'pFF': 'pFF_mean'}


def eigen(CL, V1, Q, V2):
    k10, k12, k21 = CL / V1, Q / V1, Q / V2
    s = k10 + k12 + k21; disc = np.sqrt(s * s - 4 * k10 * k21)
    return (s + disc) / 2, (s - disc) / 2           # alpha (fast), beta (terminal)


def conc(t, ka, CL, V1, Q, V2, dose=DOSE_NG):
    """plasma concentration (ng/mL) of a 2-compartment model with first-order absorption (analytic, distinct rates)"""
    t = np.asarray(t, float); k10, k12, k21 = CL / V1, Q / V1, Q / V2
    al, be = eigen(CL, V1, Q, V2)
    A = ka * (k21 - al) / ((ka - al) * (be - al)); B = ka * (k21 - be) / ((ka - be) * (al - be)); C = ka * (k21 - ka) / ((al - ka) * (be - ka))
    return dose / V1 * (A * np.exp(-al * t) + B * np.exp(-be * t) + C * np.exp(-ka * t))


def solve_V2(CL, V1, Q, beta):
    """V2 such that the terminal rate equals beta (from the source's terminal half-life)"""
    f = lambda V2: eigen(CL, V1, Q, V2)[1] - beta
    lo, hi = 1.0, 1e6
    return brentq(f, lo, hi) if f(lo) * f(hi) < 0 else np.nan


def fit_drug(name, sigma_floor=1.0):
    p = plasma[plasma.drug == name]; t = p.t_min.values / 60.0; y = p.conc_ng_per_mL.values
    CL = DOSE_NG / (nca.loc['AUC_0_inf_ng_h_per_mL', DRUGS[name]])           # CL/F = dose/AUC (mL/h)
    CL_pub = nca.loc['Cl_over_F_mL_per_h', DRUGS[name]]; beta = np.log(2) / nca.loc['t_half_h', DRUGS[name]]
    sd = np.maximum(0.12 * y, sigma_floor)

    def model(par):
        ka, V1, Q = np.exp(par); V2 = solve_V2(CL, V1, Q, beta)
        return None if not np.isfinite(V2) else (ka, V1, Q, V2)

    def resid(par):
        m = model(par)
        if m is None: return np.full(len(t) - 1, 1e3)
        ka, V1, Q, V2 = m
        try: c = conc(t[1:], ka, CL, V1, Q, V2)
        except Exception: return np.full(len(t) - 1, 1e3)
        return (c - y[1:]) / sd[1:]
    best = None
    for x0 in ([2.0, 300, 600], [4.0, 500, 1000], [1.5, 200, 400], [3.0, 800, 1500]):
        try: r = least_squares(resid, np.log(x0), bounds=(np.log([0.3, 30, 20]), np.log([30, 4000, 20000])))
        except Exception: continue
        if np.isfinite(r.cost) and (best is None or r.cost < best.cost): best = r
    ka, V1, Q, V2 = model(best.x); al, be = eigen(CL, V1, Q, V2)
    tt = np.linspace(0, 8, 4001); c = conc(tt, ka, CL, V1, Q, V2); auc = np.trapezoid(c, tt) + c[-1] / be
    return dict(drug=name, ka_per_h=ka, CL_F_mL_per_h=CL, CL_F_published=CL_pub, V1_F_mL=V1, Q_F_mL_per_h=Q, V2_F_mL=V2, Vz_F_mL=CL / be, Vz_published=nca.loc['Vd_over_F_mL', DRUGS[name]],
                t_half_h=np.log(2) / be, Cmax_model=float(c.max()), Tmax_h=float(tt[c.argmax()]), AUC_model=float(auc), chi2=float(2 * best.cost), n=len(t) - 1)


if __name__ == '__main__':
    out = {n: fit_drug(n) for n in DRUGS}
    for n, r in out.items():
        print(n, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
        p = plasma[plasma.drug == n]; ka, V1, Q, V2 = r['ka_per_h'], r['V1_F_mL'], r['Q_F_mL_per_h'], r['V2_F_mL']
        pred = conc(p.t_min.values / 60, ka, r['CL_F_mL_per_h'], V1, Q, V2)
        print('   t(min):', list(p.t_min), '\n   data  :', [round(v, 1) for v in p.conc_ng_per_mL], '\n   model :', [round(float(v), 1) for v in pred])
    json.dump(out, open(os.path.join(RES, 'pk_fit.json'), 'w'), indent=1)
