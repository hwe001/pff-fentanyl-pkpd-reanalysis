import os
import sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
from scipy import ndimage as ndi
S = sys.argv[1]
im = load(S + '/src_figs/Fig2.png'); dm = dark_mask(im)
xa, ya = find_axes(dm, x_range=(100, 260), y_range=(500, 620))
xt = x_ticks(dm, ya, depth=(5, 12), x_lo=xa - 5); yt = y_ticks(dm, xa, depth=(5, 12), y_hi=ya)
fx, ex = linmap(xt, np.arange(0, 481, 30)[:len(xt)]); fy, ey = linmap(yt, np.arange(120, -1, -20)[:len(yt)])
print('calibration residual: x %.3f min, y %.3f ng/mL' % (ex, ey))
R, G, B = im[..., 0], im[..., 1], im[..., 2]
purple = (R > 110) & (B > 150) & (G < 90) & (R < 235); green = (G > 130) & (R < 90) & (B < 90)
TIMES = [0, 30, 60, 120, 240, 480]; XPIX = [(t - np.polyfit(xt, np.arange(0, 481, 30)[:len(xt)], 1)[1]) / np.polyfit(xt, np.arange(0, 481, 30)[:len(xt)], 1)[0] for t in TIMES]


def filled(m):
    m = m.copy(); m[:, :160] = False; m[540:, :] = False; m[60:250, 410:] = False
    return ndi.binary_fill_holes(ndi.binary_closing(m, structure=np.ones((3, 3))))


out = {}
for name, m, op in (('fentanyl', purple, 11), ('pFF', green, 9)):
    f = filled(m); o = ndi.binary_opening(f, structure=np.ones((op, op))); lab, n = ndi.label(o)
    res = []
    for t, xp in zip(TIMES, XPIX):
        # prefer a clean opened blob near the expected x; otherwise use the coloured pixels in a local window
        best = None
        for k in range(1, n + 1):
            ys, xs = np.where(lab == k)
            if abs(xs.mean() - xp) < 8 and len(xs) > 80: best = (xs, ys, 'blob')
        if best is None:
            win = f[:, int(xp) - 14:int(xp) + 15]; ys, xs = np.where(win); xs = xs + int(xp) - 14
            best = (xs, ys, 'window')
        xs, ys, how = best
        cy_c, cy_b = ys.mean(), (ys.min() + ys.max()) / 2
        lo, hi = errorbar_extent(m, xp, 1)
        res.append(dict(t_min=t, y_centroid=float(fy(cy_c)), y_bbox=float(fy(cy_b)), top=float(fy(lo)), bottom=float(fy(hi)), how=how, h_px=int(ys.max() - ys.min())))
    out[name] = res
    print(name)
    for r in res: print('   t=%3d  centroid %.2f  bbox-centre %.2f  extent[%.1f, %.1f]  (%s, h=%d px)' % (r['t_min'], r['y_centroid'], r['y_bbox'], r['bottom'], r['top'], r['how'], r['h_px']))
json.dump(out, open(S + '/dig_fig2_raw.json', 'w'), indent=1)
print('\nSOURCE Table 1 Cmax: fentanyl 82.12, pFF 112.57 ng/mL at 30 min')
