import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = str(ROOT / 'data'); RES = str(ROOT / 'results'); FIG = str(ROOT / 'figures')
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
R = os.path.join(RES, 'figures')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9})
fig, ax = plt.subplots(figsize=(7.2, 3.0)); ax.set_xlim(0, 100); ax.set_ylim(0, 44); ax.axis('off')


def box(x, y, w, h, text, fc='#f2f2f2'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.0', fc=fc, ec='k', lw=0.9)); ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=8)


def arrow(x0, y0, x1, y1, both=False):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='<|-|>' if both else '-|>', mutation_scale=9, lw=0.9, color='k', shrinkA=0, shrinkB=0))


box(1, 26, 12, 9, 'sc depot'); box(21, 26, 15, 9, 'Plasma\n(central)'); box(21, 5, 15, 9, 'Peripheral', fc='#f7f7f7')
box(49, 26, 15, 9, 'Brain\n(tissue)', fc='#e8dcf5'); box(76, 26, 18, 9, 'Effect site\n(Ce)', fc='#e8dcf5')
box(40, 5, 24, 9, 'Tail flick (%MPE)', fc='#fde9d9'); box(72, 5, 25, 9, 'Hypothermia (\u0394T)', fc='#fde9d9')
arrow(13.5, 30.5, 20.5, 30.5); ax.text(17, 32, 'ka', ha='center', fontsize=7.5)
arrow(28.5, 25.5, 28.5, 14.5, both=True); ax.text(30.5, 20, 'Q', fontsize=7.5)
arrow(28.5, 36.5, 28.5, 42); ax.text(30.5, 39.5, 'CL', fontsize=7.5)
arrow(36.5, 30.5, 48.5, 30.5); ax.text(42.5, 32, 'brain/plasma\nratio (measured)', ha='center', va='bottom', fontsize=6.6)
arrow(64.5, 30.5, 75.5, 30.5); ax.text(70, 32, 'ke0', ha='center', fontsize=7.5)
arrow(85, 25.5, 85, 14.5); arrow(79, 25.5, 62, 14.5)
ax.text(1, 1.5, 'Hypothermia: \u0394T = \u2212Emax\u00b7Ce/(EC50 + Ce).   Tail flick: %MPE = 100\u00b7Ce/(EC50 + Ce).', fontsize=7, va='center')
fig.savefig(R + os.sep + 'fig1_schematic.png', dpi=300, bbox_inches='tight'); fig.savefig(R + os.sep + 'fig1_schematic.tif', dpi=600, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(R + os.sep + 'fig1_schematic.pdf', bbox_inches='tight'); print('saved')
