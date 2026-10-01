import os
"""Fig 5 of the dose-ranging paper: brain/plasma ratio bars (single-drug 100 ug/kg groups only) by region."""
import sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dose_ranging_figs'); OUT = os.path.dirname(os.path.abspath(__file__))
im = load(F + os.sep + 'p8_img0.jpeg')[0:838, 0:1350]; dm = dark_mask(im, 100)
xa = int(np.argmax(dm[:, 100:200].sum(0))) + 100
yt = sorted(y_ticks(dm, xa, depth=(6, 14), y_lo=60))[:4]; print('y-axis col', xa, 'y ticks', [round(v) for v in yt])
fy, ey = linmap(yt, [60, 40, 20, 0]); y0 = int(round(max(yt))); print('calibration residual', round(ey, 3))
R, G, B = im[..., 0], im[..., 1], im[..., 2]
purple = (abs(R - 160) < 60) & (G < 90) & (B > 150); green = (G > 140) & (R < 90) & (B < 90)
rows = []


def runs(v, mn):
    out = []; i = 0
    while i < len(v):
        if v[i]:
            j = i
            while j + 1 < len(v) and v[j + 1]: j += 1
            if j - i + 1 >= mn: out.append((i, j))
            i = j + 1
        else: i += 1
    return out


for grp, m in (('fentanyl', purple), ('pFF', green)):
    m = m.copy(); m[:, :xa + 4] = False; m[y0 + 2:, :] = False; m[:130, 640:] = False       # drop the legend
    # solid bars only (the hatched co-administration bars have gaps): require solid colour over a band just above the axis
    base = m[y0 - 40:y0 - 12, :].mean(0) > 0.85; bars = runs(base, 15)
    print(grp, 'solid bars found at columns', bars)
    for reg, (b0, b1) in zip(['Hippocampus', 'Medulla', 'Striatum', 'Frontal cortex'], bars):
        w = b1 - b0; side = list(range(b0 + 3, b0 + w // 3)) + list(range(b1 - w // 3, b1 - 2)); tops = [np.where(m[:y0, c])[0].min() for c in side if m[:y0, c].any()]
        cx = (b0 + b1) // 2; ce = np.where(m[:y0, cx - 1:cx + 2].any(1))[0]
        mean = float(fy(np.median(tops))); cap = float(fy(ce.min())); rows.append(dict(group=grp, region=reg, BP_ratio=round(mean, 2), BP_plus_sem=round(cap, 2), sem=round(cap - mean, 2)))
        print(f'  {grp:9s} {reg:15s} B:P {mean:6.2f}  +SEM {cap:6.2f}')
with open(OUT + os.sep + 'dad_brain_plasma_ratio.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
