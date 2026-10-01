import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
import sys, json, numpy as np, pandas as pd, os
import matplotlib, matplotlib.ticker; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pd_stage as P, pd_stage3 as P3
R = RES; FD = FIG
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.linewidth': 0.8})
COL = {'fentanyl': '#8a2be2', 'pFF': '#119911'}
plasma = pd.read_csv(P.D + os.sep + 'fig2_plasma_concentration.csv'); temp, mpe = P.temp, P.mpe
prof = json.load(open(R + os.sep + 'primary_fit.json')); scen = pd.read_csv(R + os.sep + 'pd_scenarios.csv'); alt = pd.read_csv(R + os.sep + 'pd_alternatives_S1.csv')
fx, fp = P3.res['exposure only'], P3.res['EC50 differs (rho)']


def curves(f, which):
    ex = [1.0, 1.0, 1.0]
    for i, j in enumerate(which): ex[j] = np.exp(f.x[5 + i])
    Emax, EC, ke, ECA, keA = np.exp(f.x[:5]); rho, eps, kap = ex; out = {}
    for g in ('fentanyl', 'pFF'):
        cb = P.brain(g, P3.KPS[g], None); pf = g == 'pFF'
        ce = P.filt(cb, ke); out[(g, 'T')] = -Emax * ce / (EC * (rho if pf else 1) + ce)
        ce = P.filt(cb, keA); out[(g, 'M')] = 100 * ce / (ECA * (rho if pf else 1) + ce)
    return out


cx, cp = curves(fx, []), curves(fp, [0])
mk = {'fentanyl': 'o', 'pFF': '^'}
# ---- Figure 2: data and fits ----
f2, a = plt.subplots(1, 3, figsize=(10.5, 3.3))
for g in ('fentanyl', 'pFF'):
    p = plasma[plasma.drug == g]; a[0].plot(p.t_min, p.conc_ng_per_mL, mk[g], color=COL[g], ms=5, label=f'{g}'); a[0].plot(P.T, P.CP[g], color=COL[g], lw=1.5)
    d = temp[temp.group == g].dropna(subset=['delta_T_C']); a[1].errorbar(d.t_min, d.delta_T_C, yerr=d.sem_estimate.clip(lower=0.2), fmt=mk[g], color=COL[g], ms=5, capsize=2)
    a[1].plot(P.T, cx[(g, 'T')], '--', color=COL[g], lw=1.2); a[1].plot(P.T, cp[(g, 'T')], '-', color=COL[g], lw=1.6)
    d = mpe[mpe.group == g].dropna(subset=['MPE_percent']); a[2].errorbar(d.t_min, d.MPE_percent, yerr=d.sem_estimate.clip(lower=0), fmt=mk[g], color=COL[g], ms=5, capsize=2)
    a[2].plot(P.T, cx[(g, 'M')], '--', color=COL[g], lw=1.2); a[2].plot(P.T, cp[(g, 'M')], '-', color=COL[g], lw=1.6)
a[0].set_ylabel('Plasma concentration (ng/mL)'); a[1].set_ylabel('Change in core temperature (\u00b0C)'); a[1].axhline(0, color='k', lw=0.5, ls=':'); a[2].set_ylabel('Tail flick (%MPE)'); a[2].set_ylim(-5, 115)
a[0].legend(frameon=False, fontsize=8); a[1].plot([], [], '--', color='0.4', label='shared PD parameters'); a[1].plot([], [], '-', color='0.4', label='pFF EC50 free'); a[1].legend(frameon=False, fontsize=7, loc='lower right')
for i, t in enumerate('abc'): a[i].set_xlabel('Time (min)'); a[i].text(-0.2, 1.05, t, transform=a[i].transAxes, fontsize=11, fontweight='bold')
f2.tight_layout(); [f2.savefig(FD + os.sep + 'fig2_fits.' + e, dpi=600 if e == 'tif' else 300, **({'pil_kwargs': {'compression': 'tiff_lzw'}} if e == 'tif' else {})) for e in ('png', 'tif')]; f2.savefig(FD + os.sep + 'fig2_fits.pdf')
# ---- Figure 3: inference ----
f3, b = plt.subplots(1, 3, figsize=(10.5, 3.3))
g_, pr_ = np.array(prof['grid']), np.array(prof['profile_chi2']); b[0].plot(g_, pr_ - pr_.min(), color='k'); b[0].axhline(3.84, color='0.5', ls='--', lw=0.8); b[0].axvline(1, color='r', lw=0.8, ls=':')
b[0].set_xscale('log'); b[0].set_xticks([0.25, 0.5, 1, 2, 4]); b[0].set_xticklabels(['0.25', '0.5', '1', '2', '4']); b[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); b[0].set_xlim(0.2, 5); b[0].set_ylim(0, 40)
b[0].set_xlabel('pFF / fentanyl EC50 ratio (\u03c1)'); b[0].set_ylabel('\u0394\u03c7\u00b2 (profile likelihood)'); b[0].text(1.06, 30, 'equal potency', color='r', fontsize=7)
SHORT = ['S1 brain/plasma ratio', 'S2 absolute brain level', 'S3 no brain difference', 'S4 brain lag 10 min', 'S5 brain lag 30 min', 'S6 hippocampus', 'S6 medulla', 'S6 striatum', 'S6 frontal cortex']
ys = np.arange(len(scen))[::-1]; b[1].errorbar(scen.rho_hat, ys, xerr=[scen.rho_hat - scen.rho_95CI_low, scen.rho_95CI_high - scen.rho_hat], fmt='s', color='k', ms=4, capsize=2); b[1].axvline(1, color='r', lw=0.8, ls=':')
b[1].set_yticks(ys); b[1].set_yticklabels(SHORT, fontsize=7); b[1].set_xlabel('pFF / fentanyl EC50 ratio (95% CI)')
names = ['Shared PD parameters', 'pFF EC50 free', 'pFF Emax free', 'pFF ke0 free', 'pFF EC50 + ke0 free']
b[2].barh(np.arange(len(alt))[::-1], alt.dAICc, color='0.6'); b[2].set_yticks(np.arange(len(alt))[::-1]); b[2].set_yticklabels(names, fontsize=7); b[2].set_xlabel('\u0394AICc (0 = best)')
for i, t in enumerate('abc'): b[i].text(-0.2 if i != 1 else -0.55, 1.05, t, transform=b[i].transAxes, fontsize=11, fontweight='bold')
f3.tight_layout(); [f3.savefig(FD + os.sep + 'fig3_inference.' + e, dpi=600 if e == 'tif' else 300, **({'pil_kwargs': {'compression': 'tiff_lzw'}} if e == 'tif' else {})) for e in ('png', 'tif')]; f3.savefig(FD + os.sep + 'fig3_inference.pdf')
# ---- graphical contents entry (50 x 60 mm) ----
kp_r = P3.KPS['pFF'] / P3.KPS['fentanyl']; rho = float(np.exp(prof['x_potency'][5])); eff = kp_r / rho
g, ax = plt.subplots(figsize=(50 / 25.4, 60 / 25.4)); ax.bar([0, 1], [kp_r, eff], color=['#8a2be2', '#119911'], width=0.6); ax.axhline(1, color='k', lw=0.6, ls=':')
ax.set_xticks([0, 1]); ax.set_xticklabels(['measured\nbrain ratio', 'effective\n(fits effect)'], fontsize=6); ax.set_ylabel('pFF / fentanyl exposure', fontsize=6); ax.tick_params(labelsize=6)
for x, v in zip([0, 1], [kp_r, eff]): ax.text(x, v + 0.05, f'{v:.1f}', ha='center', fontsize=7)
ax.set_ylim(0, 2.8); g.tight_layout(pad=0.4); g.savefig(FD + os.sep + 'graphical_contents.png', dpi=600); g.savefig(FD + os.sep + 'graphical_contents.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
print('saved figures; measured brain ratio %.2f, rho %.2f, effective %.2f' % (kp_r, rho, eff))
