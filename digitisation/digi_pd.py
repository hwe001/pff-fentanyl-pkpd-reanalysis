import os
import sys, json, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi import *
from scipy import ndimage as ndi
from PIL import Image, ImageDraw
S = sys.argv[1]
TRI_OFFSET = 0.16          # fraction of marker height: data point sits this far BELOW the bounding-box centre of a down-triangle (calibrated on the exact Cmax of Fig 2)
TIMES = [30, 60, 120, 240, 480]

CFG = {
    'Fig3_tailflick': dict(file='Fig3.png', crop=None, ax_x=(100, 230), ax_y=(520, 640), xr=(0, None), legend=(1000, 0, 1501, 330),
                           ytick_vals={56.5: 125, 210.5: 75, 365.5: 25, 442.5: 0, 519.5: -25}, xtick_vals='0-480/60'),
    'Fig4A_temperature': dict(file='Fig4.png', crop=810, ax_x=(150, 260), ax_y=(630, 700), legend=(1040, 80, 1502, 330),
                              ytick_vals=None, xtick_vals='0-480/30'),
}
rows_out = []
for name, c in CFG.items():
    im = load(f"{S}/src_figs/{c['file']}")
    if c['crop']: im = im[:c['crop']]
    dm = dark_mask(im)
    xa, ya = find_axes(dm, x_range=c['ax_x'], y_range=c['ax_y'])
    xt = x_ticks(dm, ya, depth=(5, 12), x_lo=xa - 5)
    yt = y_ticks(dm, xa, depth=(5, 12), y_hi=ya + 3 if name.startswith('Fig3') else ya - 20)
    print(f'\n=== {name}: y-axis col {xa}, x-axis row {ya}, {len(xt)} x ticks, y ticks at {[round(v, 1) for v in yt]}')
    if name.startswith('Fig3'):
        xv = np.arange(0, 481, 60)[:len(xt)]
        pv = [(p, c['ytick_vals'][min(c['ytick_vals'], key=lambda k: abs(k - p))]) for p in yt]
    else:
        xv = np.arange(0, 481, 30)[:len(xt)] if len(xt) > 10 else np.arange(0, 481, 60)[:len(xt)]
        pv = list(zip(yt, [2, 0, -2, -4, -6, -8][:len(yt)]))
    fx, ex = linmap(xt, xv); fy, ey = linmap([p for p, _ in pv], [v for _, v in pv])
    print(f'  calibration residuals: x {ex:.3f} min, y {ey:.3f}')
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    masks = {'control': (R > 180) & (G < 120) & (B < 120), 'fentanyl': (R > 110) & (B > 150) & (G < 90) & (R < 235), 'pFF': (G > 130) & (R < 90) & (B < 90)}
    lx0, ly0, lx1, ly1 = c['legend']
    img = Image.open(f"{S}/src_figs/{c['file']}").convert('RGB'); d = ImageDraw.Draw(img)
    a_x, b_x = np.polyfit(xt, xv, 1); a_y, b_y = np.polyfit([p for p, _ in pv], [v for _, v in pv], 1)
    for g, m in masks.items():
        m = m.copy(); m[:, :xa + 4] = False; m[ya - 3:, :] = False; m[ly0:ly1, lx0:lx1] = False
        f = ndi.binary_fill_holes(ndi.binary_closing(m, structure=np.ones((3, 3))))
        for t in TIMES:
            xp = (t - b_x) / a_x; best = None; used = None
            for op in ((11, 9, 7, 5) if g == 'pFF' else (13, 11, 9, 7, 5)):
                o = ndi.binary_opening(f, structure=np.ones((op, op))); lab, n = ndi.label(o)
                for k in range(1, n + 1):
                    ys, xs = np.where(lab == k)
                    if abs(xs.mean() - xp) < 9 and len(xs) > 40 and (best is None or len(xs) > len(best[0])): best = (xs, ys)
                if best is not None: used = op; break
            if best is None:
                rows_out.append(dict(figure=name, group=g, t_min=t, value=None, top=None, bottom=None, h_px=None, note='hidden behind another marker')); continue
            xs, ys = best; h = ys.max() - ys.min(); cy = (ys.min() + ys.max()) / 2 + (TRI_OFFSET * h if g == 'pFF' else 0)
            lo, hi = errorbar_extent(m, xp, 1)
            small = used < (11 if g == 'pFF' else 13)
            rows_out.append(dict(figure=name, group=g, t_min=t, value=round(float(fy(cy)), 2), top=round(float(fy(lo)), 2), bottom=round(float(fy(hi)), 2), h_px=int(h), note='partly overlapped (lower accuracy)' if small else ''))
            col = {'control': (0, 0, 0), 'fentanyl': (255, 140, 0), 'pFF': (0, 0, 255)}[g]
            d.ellipse((xp - 7, cy - 7, xp + 7, cy + 7), outline=col, width=3)
    img.save(f'{S}/dig_{name}_check.png')
with open(S + '/dig_pd_raw.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0])); w.writeheader(); w.writerows(rows_out)
for r in rows_out: print(r['figure'][:6], r['group'].ljust(8), str(r['t_min']).rjust(4), r['value'], '[%s, %s]' % (r['bottom'], r['top']), r['note'])
