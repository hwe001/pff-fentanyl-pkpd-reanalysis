import os
"""Digitise Canfield, Fry & Sprague (Drug Alcohol Depend 2025;272:112710; CC BY-NC-ND) Figs 2-5: plasma, temperature, tail flick at 50/100/200 ug/kg."""
import sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
from scipy import ndimage as ndi
from PIL import Image, ImageDraw
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dose_ranging_figs'); OUT = os.path.dirname(os.path.abspath(__file__))
PAL = {'purple': (150, 30, 185), 'blue': (20, 20, 175), 'green': (35, 170, 35), 'red': (205, 75, 75), 'orange': (205, 130, 30)}
TIMES = [30, 60, 120, 240, 480]
# panel: file, crop, y-axis search window, x-axis search window, y tick values (top to bottom), x tick values, series colour -> label, legend box to blank (relative to crop)
P = {
 'Fig3A_temp_fentanyl': dict(f='p5_img0.jpeg', crop=(0, 0, 860, 640), xw=(100, 230), yw=(480, 560), yv=[2, 0, -2, -4], xv=list(range(0, 481, 60)), ser={'purple': 'fentanyl_200', 'blue': 'fentanyl_100', 'green': 'fentanyl_50', 'red': 'control'}, leg=(380, 30, 860, 190)),
 'Fig3B_temp_pFF': dict(f='p5_img0.jpeg', crop=(860, 0, 1800, 640), xw=(100, 230), yw=(480, 560), yv=[2, 0, -2, -4], xv=list(range(0, 481, 60)), ser={'orange': 'pFF_200', 'blue': 'pFF_100', 'green': 'pFF_50', 'red': 'control'}, leg=(380, 30, 940, 190)),
 'Fig4A_TF_fentanyl': dict(f='p6_img0.jpeg', crop=(0, 0, 900, 790), xw=(100, 230), yw=(520, 620), yv=[125, 75, 25, 0, -25], xv=list(range(0, 481, 60)), ser={'purple': 'fentanyl_200', 'blue': 'fentanyl_100', 'green': 'fentanyl_50', 'red': 'control'}, leg=(380, 30, 900, 230)),
 'Fig4B_TF_pFF': dict(f='p6_img0.jpeg', crop=(900, 0, 1800, 790), xw=(100, 230), yw=(520, 620), yv=[125, 75, 25, 0, -25], xv=list(range(0, 481, 60)), ser={'orange': 'pFF_200', 'blue': 'pFF_100', 'green': 'pFF_50', 'red': 'control'}, leg=(380, 30, 900, 230)),
}


def classify(im):
    sat = (im.max(2) - im.min(2)) > 55; cols = np.array(list(PAL.values())); names = list(PAL)
    d = np.linalg.norm(im[..., None, :] - cols[None, None], axis=3); k = d.argmin(2); ok = sat & (d.min(2) < 120)
    return {n: (ok & (k == i)) for i, n in enumerate(names)}


rows = []
for pname, c in P.items():
    x0, y0, x1, y1 = c['crop']; img = load(F + os.sep + c['f'])[y0:y1, x0:x1]; dm = dark_mask(img, 100)
    xa, ya = find_axes(dm, x_range=c['xw'], y_range=c['yw'])
    xt = x_ticks(dm, ya, depth=(6, 14), x_lo=xa - 5); yt = y_ticks(dm, xa, depth=(6, 14), y_hi=ya + 4)
    fx, ex = linmap(xt, c['xv'][:len(xt)]); fy, ey = linmap(yt, c['yv'][:len(yt)]); ax_, bx_ = np.polyfit(xt, c['xv'][:len(xt)], 1); ay_, by_ = np.polyfit(yt, c['yv'][:len(yt)], 1)
    masks = classify(img); pil = Image.open(F + os.sep + c['f']).convert('RGB').crop(c['crop']); dr = ImageDraw.Draw(pil)
    print(f'{pname}: calibration residual x {ex:.2f} min, y {ey:.3f}')
    for col, label in c['ser'].items():
        m = masks[col].copy(); m[:, :xa + 4] = False; m[ya - 2:, :] = False; lx0, ly0, lx1, ly1 = c['leg']; m[ly0:ly1, lx0:lx1] = False
        f = ndi.binary_fill_holes(ndi.binary_closing(m, structure=np.ones((3, 3))))
        for t in TIMES:
            xp = (t - bx_) / ax_; best = None; used = None
            for op in (13, 11, 9, 7, 5):
                o = ndi.binary_opening(f, structure=np.ones((op, op))); lab, n = ndi.label(o)
                for k in range(1, n + 1):
                    ys, xs = np.where(lab == k)
                    if abs(xs.mean() - xp) < 11 and len(xs) > 35 and (best is None or len(xs) > len(best[0])): best = (xs, ys)
                if best is not None: used = op; break
            if best is None: rows.append(dict(panel=pname, series=label, t_min=t, value='', note='hidden')); continue
            xs, ys = best; cy = (ys.min() + ys.max()) / 2 + (0.16 * (ys.max() - ys.min()) if col == 'green' else 0)
            strip = m[:, int(xp) - 1:int(xp) + 2].any(1); yy = np.where(strip)[0]
            rows.append(dict(panel=pname, series=label, t_min=t, value=round(float(fy(cy)), 2), lo=round(float(fy(yy.max())), 2), hi=round(float(fy(yy.min())), 2), note='partly overlapped' if used < 11 else ''))
            dr.ellipse((xp - 8, cy - 8, xp + 8, cy + 8), outline=(0, 0, 0), width=3)

# ---- Fig 2: filled symbols only (single-drug groups) ----
img = load(F + os.sep + 'p4_img0.jpeg')[0:745, 0:1050]; dm = dark_mask(img, 100); xa, ya = find_axes(dm, x_range=(120, 260), y_range=(560, 640))
xt = x_ticks(dm, ya, depth=(6, 14), x_lo=xa - 5); yt = y_ticks(dm, xa, depth=(6, 14), y_hi=ya + 4)
fy, ey = linmap(yt, [40, 30, 20, 10, 0][:len(yt)]); ax_, bx_ = np.polyfit(xt, list(range(0, 481, 30))[:len(xt)], 1); print(f'Fig2: calibration residual y {ey:.3f}')
masks = classify(img); pil = Image.open(F + os.sep + 'p4_img0.jpeg').convert('RGB'); dr = ImageDraw.Draw(pil)
for col, label in (('purple', 'fentanyl_100'), ('green', 'pFF_100')):
    m = masks[col].copy(); m[:, :xa + 4] = False; m[ya - 2:, :] = False; m[0:200, 330:] = False
    f = ndi.binary_fill_holes(ndi.binary_closing(m, structure=np.ones((3, 3))))
    for t in TIMES:
        xp = (t - bx_) / ax_; cand = []
        for op in (11, 9, 7):
            o = ndi.binary_opening(f, structure=np.ones((op, op))); lab, n = ndi.label(o)
            for k in range(1, n + 1):
                ys, xs = np.where(lab == k)
                if abs(xs.mean() - xp) < 12 and len(xs) > 60:
                    bb = m[ys.min():ys.max() + 1, xs.min():xs.max() + 1]; cand.append((bb.mean(), (ys.min() + ys.max()) / 2, op))     # fill ratio: filled markers are solid
            if cand: break
        if not cand: rows.append(dict(panel='Fig2_plasma', series=label, t_min=t, value='', note='hidden')); continue
        cand.sort(reverse=True); fr, cy, op = cand[0]
        rows.append(dict(panel='Fig2_plasma', series=label, t_min=t, value=round(float(fy(cy)), 2), note=('' if fr > 0.55 else 'hollow marker chosen (filled one hidden)')))
        dr.ellipse((xp - 9, cy - 9, xp + 9, cy + 9), outline=(0, 0, 0), width=3)
with open(OUT + os.sep + 'dad_dose_data_raw.csv', 'w', newline='') as fh:
    keys = ['panel', 'series', 't_min', 'value', 'lo', 'hi', 'note']; w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
import collections
for r in rows:
    if r['series'] in ('fentanyl_100', 'pFF_100'): print(r['panel'].ljust(20), r['series'].ljust(13), str(r['t_min']).rjust(4), r['value'], r['note'])
