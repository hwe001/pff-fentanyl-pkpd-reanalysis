import os
import sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
from scipy import ndimage as ndi
from PIL import Image, ImageDraw
S = sys.argv[1]
im_full = load(S + '/src_figs/Fig5.png'); H = im_full.shape[0]
REG = ['Hippocampus', 'Medulla', 'Striatum', 'Frontal Cortex']
panels = {'A_conc_ng_per_g': dict(rows=(0, 880), yvals=[0, 1000, 2000, 3000, 4000], legend=(520, 0, 970, 140)),
          'B_brain_plasma_ratio': dict(rows=(880, H), yvals=[0, 20, 40, 60, 80], legend=(0, 0, 0, 0))}
out_bars, out_dots = [], []
img = Image.open(S + '/src_figs/Fig5.png').convert('RGB'); d = ImageDraw.Draw(img)


def runs(v, min_len):
    out = []; i = 0
    while i < len(v):
        if v[i]:
            j = i
            while j + 1 < len(v) and v[j + 1]: j += 1
            if j - i + 1 >= min_len: out.append((i, j))
            i = j + 1
        else: i += 1
    return out


for pname, p in panels.items():
    r0, r1 = p['rows']; im = im_full[r0:r1]; dm = dark_mask(im)
    xa = int(np.argmax(dm[:, 120:260].sum(0))) + 120
    yt = sorted(y_ticks(dm, xa, depth=(5, 12), y_lo=20))[:len(p['yvals'])]
    ya = int(round(max(yt))); fy, ey = linmap(yt, sorted(p['yvals'], reverse=True)); print(f'\n== {pname}: y-axis col {xa}, zero row {ya}, calibration residual {ey:.3f}')
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    black = (im.max(2) < 120); black[:, :xa + 4] = False; black[ya - 1:, :] = False
    for grp, m in (('fentanyl', (R > 110) & (B > 150) & (G < 90) & (R < 235)), ('pFF', (G > 130) & (R < 90) & (B < 90))):
        m = m.copy(); m[:, :xa + 4] = False; m[ya + 1:, :] = False
        lx0, ly0, lx1, ly1 = p['legend']; m[ly0:ly1, lx0:lx1] = False
        base = m[ya - 18:ya - 8, :].all(0)          # bar columns: coloured down to just above the (dark) axis line
        bars = runs(base, 35)
        print('  ', grp, 'bars', bars)
        for reg, (bx0, bx1) in zip(REG, bars):
            w = bx1 - bx0; side = list(range(bx0 + 3, bx0 + w // 4)) + list(range(bx1 - w // 4, bx1 - 2))      # columns away from the centre (error bar) and the dots
            tops = [np.where(m[:ya, c])[0].min() for c in side if m[:ya, c].any()]
            top = float(np.median(tops)); mean = float(fy(top))
            cx = (bx0 + bx1) // 2; ce = np.where(m[:ya, cx - 1:cx + 2].any(1))[0]
            cap = float(fy(ce.min())); out_bars.append(dict(panel=pname, group=grp, region=reg, mean=round(mean, 2), mean_plus_sem=round(cap, 2), sem=round(cap - mean, 2)))
            sub = np.zeros_like(black); sub[:, bx0 - 8:bx1 + 8] = black[:, bx0 - 8:bx1 + 8]
            sub = ndi.binary_opening(sub, structure=np.ones((7, 7))); lb, nb = ndi.label(sub); vals = []
            for kk in range(1, nb + 1):
                yy, xx = np.where(lb == kk)
                if 50 < len(xx) < 450: vals.append(float(fy(yy.mean()))); d.ellipse((xx.mean() - 5, yy.mean() + r0 - 5, xx.mean() + 5, yy.mean() + r0 + 5), outline=(255, 0, 0), width=2)
            for v in sorted(vals, reverse=True): out_dots.append(dict(panel=pname, group=grp, region=reg, value=round(v, 2)))
            print(f'   {reg:14s} mean {mean:8.2f}  +SEM {cap:8.2f}  points {len(vals)}  mean of points {np.mean(vals) if vals else float("nan"):8.2f}')
            d.rectangle((bx0, top + r0, bx1, ya + r0), outline=(0, 0, 255), width=2)
img.save(S + '/dig_Fig5_check.png')
for fn, rows in (('dig_brain_bars.csv', out_bars), ('dig_brain_points.csv', out_dots)):
    with open(S + '/' + fn, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('\nsaved bars', len(out_bars), 'points', len(out_dots))
