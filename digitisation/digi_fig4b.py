import os
import sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
from scipy import ndimage as ndi
from PIL import Image, ImageDraw
S = sys.argv[1]
full = load(S + '/src_figs/Fig4.png'); r0 = 960; im = full[r0:]
dm = dark_mask(im)
xa = int(np.argmax(dm[:, 150:260].sum(0))) + 150
yt = sorted(y_ticks(dm, xa, depth=(5, 12), y_lo=20))
print('y-axis col', xa, 'y ticks', [round(v, 1) for v in yt])
yt = yt[:5]; fy, ey = linmap(yt, [0, -2, -4, -6, -8]); print('calibration residual', round(ey, 3))
y0 = int(round(yt[0]))
R, G, B = im[..., 0], im[..., 1], im[..., 2]
black = (im.max(2) < 120); black[:y0 + 3, :] = False; black[:, :xa + 4] = False; black[:, 700:] = False
out_b, out_p = [], []
img = Image.open(S + '/src_figs/Fig4.png').convert('RGB'); d = ImageDraw.Draw(img)
for grp, m in (('fentanyl', (R > 110) & (B > 150) & (G < 90) & (R < 235)), ('pFF', (G > 130) & (R < 90) & (B < 90))):
    m = m.copy(); m[:y0 + 3, :] = False; m[:, :xa + 4] = False; m[:, 560:] = False
    band = m[y0 + 12:y0 + 40, :].all(0); cols = np.where(band)[0]
    if len(cols) == 0: print(grp, 'no bar'); continue
    bx0, bx1 = cols.min(), cols.max(); w = bx1 - bx0
    side = list(range(bx0 + 3, bx0 + w // 4)) + list(range(bx1 - w // 4, bx1 - 2))
    bottoms = [np.where(m[:, c])[0].max() for c in side if m[:, c].any()]
    bot = float(np.median(bottoms)); mean = float(fy(bot))
    cx = (bx0 + bx1) // 2; ce = np.where(m[:, cx - 1:cx + 2].any(1))[0]; cap = float(fy(ce.max()))
    sub = np.zeros_like(black); sub[:, bx0 - 8:bx1 + 8] = black[:, bx0 - 8:bx1 + 8]
    sub = ndi.binary_opening(sub, structure=np.ones((7, 7))); lb, nb = ndi.label(sub); vals = []
    for kk in range(1, nb + 1):
        yy, xx = np.where(lb == kk)
        if 50 < len(xx) < 450: vals.append(float(fy(yy.mean()))); d.ellipse((xx.mean() - 5, yy.mean() + r0 - 5, xx.mean() + 5, yy.mean() + r0 + 5), outline=(255, 0, 0), width=2)
    out_b.append(dict(group=grp, mean_max_dT=round(mean, 2), mean_minus_sem=round(cap, 2), sem=round(mean - cap, 2)))
    for v in sorted(vals): out_p.append(dict(group=grp, max_dT=round(v, 2)))
    print(f'{grp}: bar columns {bx0}-{bx1}; mean max dT {mean:.2f} C; mean-SEM {cap:.2f}; {len(vals)} points, mean of points {np.mean(vals):.2f}: {[round(v, 2) for v in sorted(vals)]}')
img.save(S + '/dig_Fig4B_check.png')
for fn, rows in (('dig_fig4B_bars.csv', out_b), ('dig_fig4B_points.csv', out_p)):
    with open(S + '/' + fn, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
