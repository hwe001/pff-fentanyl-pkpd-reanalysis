import os
"""Small figure digitiser for GraphPad-style plots: axis calibration from tick marks, marker centres from colour segmentation."""
import numpy as np
from PIL import Image
from scipy import ndimage as ndi


def load(path): return np.asarray(Image.open(path).convert('RGB')).astype(int)


def dark_mask(im, thr=90): return (im.max(2) < thr)


def find_axes(dm, x_range=None, y_range=None):
    """Return (x_of_yaxis, y_of_xaxis): the column / row with the longest dark run in the allowed window."""
    H, W = dm.shape
    x0, x1 = x_range or (0, W); y0, y1 = y_range or (0, H)
    col = dm[:, x0:x1].sum(0); row = dm[y0:y1, :].sum(1)
    return x0 + int(np.argmax(col)), y0 + int(np.argmax(row))


def clusters(v):
    """centres of runs of True in a 1-D boolean array"""
    out = []; i = 0
    while i < len(v):
        if v[i]:
            j = i
            while j + 1 < len(v) and v[j + 1]: j += 1
            out.append((i + j) / 2); i = j + 1
        else: i += 1
    return out


def x_ticks(dm, yaxis_row, depth=(6, 14), x_lo=0, x_hi=None, below=True):
    sgn = 1 if below else -1
    rows = [yaxis_row + sgn * d for d in range(depth[0], depth[1])]
    v = dm[rows, :].all(0)
    if x_hi is None: x_hi = dm.shape[1]
    v[:x_lo] = False; v[x_hi:] = False
    return clusters(v)


def y_ticks(dm, xaxis_col, depth=(6, 14), y_lo=0, y_hi=None):
    cols = [xaxis_col - d for d in range(depth[0], depth[1])]
    v = dm[:, cols].all(1)
    if y_hi is None: y_hi = dm.shape[0]
    v[:y_lo] = False; v[y_hi:] = False
    return clusters(v)


def linmap(pix, val):
    a, b = np.polyfit(pix, val, 1); res = np.array(val) - (a * np.array(pix) + b)
    return (lambda p: a * np.asarray(p) + b), float(np.abs(res).max())


def colour_mask(im, rgb, tol=60):
    return (np.abs(im - np.array(rgb)).sum(2) < tol)


def marker_centres(mask, open_px=7, min_area=60):
    """remove thin lines by opening, label what is left, return (cx, cy, area, bbox)"""
    op = ndi.binary_opening(mask, structure=np.ones((open_px, open_px)))
    lab, n = ndi.label(op); out = []
    for k in range(1, n + 1):
        ys, xs = np.where(lab == k)
        if len(xs) < min_area: continue
        out.append(dict(cx=xs.mean(), cy=ys.mean(), area=len(xs), x0=xs.min(), x1=xs.max(), y0=ys.min(), y1=ys.max(),
                        bx=(xs.min() + xs.max()) / 2, by=(ys.min() + ys.max()) / 2))
    return sorted(out, key=lambda d: d['cx'])


def errorbar_extent(mask, cx, half_w=2):
    """vertical extent (rows) of coloured pixels in the column strip at x = cx: includes marker and any visible error bar"""
    strip = mask[:, int(cx) - half_w:int(cx) + half_w + 1].any(1)
    ys = np.where(strip)[0]
    return (int(ys.min()), int(ys.max())) if len(ys) else (None, None)
