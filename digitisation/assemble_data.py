import sys, json, csv, shutil, os
import pandas as pd, numpy as np
S = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(S)), 'data')

# ---- Fig 2: plasma ----
raw = json.load(open(S + '/dig_fig2_raw.json')); rows = []
TRI = 0.74        # ng/mL: pFF triangle offset, calibrated on the exact Cmax (112.57) of Table 1
for g in ('fentanyl', 'pFF'):
    for r in raw[g]:
        t = r['t_min']; v = 0.0 if t == 0 else r['y_bbox'] - (TRI if g == 'pFF' else 0)
        bar = max(v - r['bottom'], r['top'] - v) if t else float('nan')
        note = 'baseline (zero by design)' if t == 0 else ('partly hidden by the other marker; accuracy about +/-1 ng/mL' if r['how'] == 'window' else '')
        rows.append(dict(drug=g, t_min=t, conc_ng_per_mL=round(max(v, 0), 2), error_bar_extent=round(bar, 2) if bar == bar else '', note=note))
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'fig2_plasma_concentration.csv'), index=False)

# ---- Fig 3 and 4A ----
pdx = pd.read_csv(S + '/dig_pd_raw.csv')
pdx['sem_estimate'] = np.maximum(pdx.value - pdx.bottom, pdx.top - pdx.value).round(2)
base = []
for fig in pdx.figure.unique():
    for g in pdx.group.unique(): base.append(dict(figure=fig, group=g, t_min=0, value=0.0, top=0, bottom=0, h_px='', note='baseline (zero by design)', sem_estimate=0.0))
pdx = pd.concat([pd.DataFrame(base), pdx], ignore_index=True)
tf = pdx[pdx.figure.str.startswith('Fig3')].copy()
tf.loc[(tf.group == 'pFF') & tf.value.isna(), ['value', 'note']] = [100.0, 'hidden behind the fentanyl marker; both groups at the 100% ceiling at 30-120 min']
tf.loc[(tf.group == 'fentanyl') & (tf.t_min.isin([30, 60, 120])), 'value'] = 100.0
tf.loc[(tf.group == 'fentanyl') & (tf.t_min.isin([30, 60, 120])), 'note'] = 'at the 100% ceiling (digitised 100.2)'
tf.rename(columns={'value': 'MPE_percent'}).drop(columns=['figure', 'h_px']).sort_values(['group', 't_min']).to_csv(os.path.join(OUT, 'fig3_tail_flick_MPE.csv'), index=False)
te = pdx[pdx.figure.str.startswith('Fig4A')].copy()
te.rename(columns={'value': 'delta_T_C'}).drop(columns=['figure', 'h_px']).sort_values(['group', 't_min']).to_csv(os.path.join(OUT, 'fig4A_temperature_change.csv'), index=False)

# ---- Fig 4B and Fig 5 ----
shutil.copy(S + '/dig_fig4B_points.csv', os.path.join(OUT, 'fig4B_max_temperature_change_individual_animals.csv'))
shutil.copy(S + '/dig_brain_points.csv', os.path.join(OUT, 'fig5_brain_individual_animals.csv'))
b = pd.read_csv(S + '/dig_brain_bars.csv'); b.to_csv(os.path.join(OUT, 'fig5_brain_group_means.csv'), index=False)

# ---- Table 1 of the source (exact values, from the article text) ----
t1 = pd.DataFrame([['t_half_h', 1.10, 0.01, 1.12, 0.03], ['AUC_0_inf_ng_h_per_mL', 112.83, 10.82, 158.77, 21.60], ['Cmax_ng_per_mL', 82.12, 11.78, 112.57, 16.26],
                   ['Tmax_h', 0.5, '', 0.5, ''], ['Cl_over_F_mL_per_h', 917.60, 110.11, 697.02, 109.87], ['Vd_over_F_mL', 1524.90, 368.52, 1121.11, 168.08]],
                  columns=['parameter', 'fentanyl_mean', 'fentanyl_SEM', 'pFF_mean', 'pFF_SEM']); t1.to_csv(os.path.join(OUT, 'table1_noncompartmental_PK_exact_from_source.csv'), index=False)
print('files:', sorted(os.listdir(OUT)))
