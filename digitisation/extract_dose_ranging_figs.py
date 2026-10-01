"""Extract the figure images from the dose-ranging article's PDF (Canfield JR, Fry KK, Sprague JE. Drug Alcohol Depend 2025;272:112710, CC BY-NC-ND).
Download the open-access PDF yourself and run:  python digitisation/extract_dose_ranging_figs.py "path/to/article.pdf"
The images are written to digitisation/dose_ranging_figs/ (not committed: the article's licence does not allow redistribution of altered or adapted material)."""
import sys, os, fitz
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dose_ranging_figs'); os.makedirs(out, exist_ok=True)
d = fitz.open(sys.argv[1])
for i, p in enumerate(d):
    for k, img in enumerate(p.get_images(full=True)):
        x = d.extract_image(img[0])
        if x['width'] >= 800: open(os.path.join(out, f"p{i + 1}_img{k}.{x['ext']}"), 'wb').write(x['image']); print('wrote', f"p{i + 1}_img{k}.{x['ext']}", x['width'], 'x', x['height'])
