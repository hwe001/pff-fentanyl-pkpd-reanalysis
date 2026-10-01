import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Summary numbers of the two-dose analysis, saved for the manuscript: out-of-sample chi-squared, per-dose split of the joint fits, joint parameters, 100 ug/kg PK."""
import json, numpy as np
import dose_test as D
R = D.R
prim = json.load(open(R + os.sep + 'primary_fit.json')); xs, xp = np.array(prim['x_exposure']), np.array(prim['x_potency'])
out = {'oos': {}, 'joint': {}, 'pk100': {d: {k: float(v) for k, v in p.items()} for d, p in D.pk100.items()}, 'Kp100': {k: float(v) for k, v in D.BP100.items()}, 'Kp300': {k: float(v) for k, v in D.KP300.items()}}
for ctrl in (False, True):
    o100 = D.obs_set(False, True, ctrl)
    out['oos']['ctrl' + str(int(ctrl))] = {'shared': float((D.resid(xs, o100, False) ** 2).sum()), 'EC50free': float((D.resid(xp, o100, True) ** 2).sum()), 'n': len(o100)}
    j = json.load(open(R + f'\\dose_test_joint_ctrl{int(ctrl)}.json')); d = {}
    for name, x, rf in (('shared', np.array(j['x_shared']), False), ('EC50free', np.array(j['x_potency']), True)):
        d[name] = {'chi2_300': float((D.resid(x, D.obs_set(True, False, ctrl), rf) ** 2).sum()), 'chi2_100': float((D.resid(x, D.obs_set(False, True, ctrl), rf) ** 2).sum()),
                   'params': dict(zip(['Emax_T', 'EC50_T', 'ke0_T', 'EC50_A', 'ke0_A', 'rho'], [float(v) for v in np.exp(x)]))}
    out['joint']['ctrl' + str(int(ctrl))] = d
json.dump(out, open(R + os.sep + 'dose_test_summary.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('oos',)}, indent=1)); print(json.dumps(out['joint']['ctrl0']['EC50free'], indent=1)); print('Kp100', out['Kp100'], '| PK100 ka/V1', {d: (round(p['ka'], 2), round(p['V1'])) for d, p in out['pk100'].items()})
