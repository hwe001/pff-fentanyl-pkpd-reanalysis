import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
"""Figure 4: the two-dose test. (a,b) 100 ug/kg data against the joint fits; (c) profile likelihood for rho from the joint fit."""
import json, os, numpy as np, pandas as pd
import matplotlib, matplotlib.ticker; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import dose_test as D
R = D.R; FD = FIG; T = D.T
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.linewidth': 0.8})
COL = {'fentanyl': '#8a2be2', 'pFF': '#119911'}; MK = {'fentanyl': 'o', 'pFF': '^'}
j = json.load(open(R + os.sep + 'dose_test_joint_ctrl0.json')); xs, xp = np.array(j['x_shared']), np.array(j['x_potency'])
obs = D.OBS100
fig, a = plt.subplots(1, 3, figsize=(10.5, 3.3))
for drug, grp in (('fentanyl', 'f100'), ('pFF', 'p100')):
    for ax, kind, key in ((a[0], 'T', 'T'), (a[1], 'M', 'M')):
        pts = [(t, y, s) for k, g, t, y, s in obs if g == drug and (k == key or (key == 'M' and k == 'Mc'))]
        t_, y_, s_ = zip(*pts); ax.errorbar(t_, y_, yerr=s_, fmt=MK[drug], color=COL[drug], ms=5, capsize=2, lw=0.8)
        ax.plot(T, [D.effect(xs, grp, False, 'T' if key == 'T' else 'M', t) for t in T], '--', color=COL[drug], lw=1.2)
        ax.plot(T, [D.effect(xp, grp, True, 'T' if key == 'T' else 'M', t) for t in T], '-', color=COL[drug], lw=1.6)
a[0].set_ylabel('Change in core temperature (°C)'); a[0].axhline(0, color='k', lw=0.5, ls=':'); a[1].set_ylabel('Tail flick (%MPE)'); a[1].set_ylim(-5, 120)
a[0].plot([], [], '--', color='0.4', label='shared PD parameters'); a[0].plot([], [], '-', color='0.4', label='pFF EC50 free'); a[0].legend(frameon=False, fontsize=7, loc='lower right')
for ax in a[:2]: ax.set_xlabel('Time (min)')
g_, pr_ = np.array(j['grid']), np.array(j['profile']); a[2].plot(g_, pr_ - pr_.min(), color='k'); a[2].axhline(3.84, color='0.5', ls='--', lw=0.8); a[2].axvline(1, color='r', lw=0.8, ls=':')
a[2].set_xscale('log'); a[2].set_xticks([0.5, 1, 2, 4]); a[2].set_xticklabels(['0.5', '1', '2', '4']); a[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a[2].set_xlim(0.3, 4); a[2].set_ylim(0, 40)
a[2].set_xlabel('pFF / fentanyl EC50 ratio (ρ), both doses'); a[2].set_ylabel('Δχ² (profile likelihood)'); a[2].text(1.05, 30, 'equal potency', color='r', fontsize=7)
for i, t in enumerate('abc'): a[i].text(-0.2, 1.05, t, transform=a[i].transAxes, fontsize=11, fontweight='bold')
fig.tight_layout(); fig.savefig(FD + os.sep + 'fig4_dose_test.png', dpi=300); fig.savefig(FD + os.sep + 'fig4_dose_test.pdf'); fig.savefig(FD + os.sep + 'fig4_dose_test.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'}); print('saved Figure 4')
