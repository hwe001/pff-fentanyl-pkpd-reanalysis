# pFF versus fentanyl: model-based re-analysis of published rat data

Data and code for the manuscript *Brain exposure does not account proportionally for the greater effects of para-fluorofentanyl than fentanyl in rats: a model-based re-analysis of published data* (in preparation for *Drug Testing and Analysis*).

All data come from Canfield JR, Sprague JE. *In vivo pharmacokinetic, pharmacodynamic and brain concentration comparison of fentanyl and para-fluorofentanyl in rats.* Arch Toxicol 2025;99:287-297, [doi:10.1007/s00204-024-03887-z](https://doi.org/10.1007/s00204-024-03887-z) (open access, CC BY 4.0). **No new animal experiments were performed, and nothing here is unpublished data**: the values were read from the figures of that article and Table 1.

## Question

The source study found that para-fluorofentanyl (pFF) reached more than twice the brain concentration of fentanyl at a similar plasma half-life, with deeper hypothermia and similar analgesia. Does the higher brain exposure alone explain the greater effects? A pharmacokinetic-pharmacodynamic model is fitted twice: with all pharmacodynamic parameters shared by the two drugs (differences arise only from exposure), and with a separate pFF EC50.

## Model

Two-compartment plasma model with first-order absorption (clearance and half-life fixed to the published noncompartmental values) -> brain concentration = plasma concentration x measured brain/plasma ratio -> effect compartment (ke0) -> hypothermia and tail-flick response (saturable). Maximal hypothermia is bounded at 10 degC (a run with a 20 degC bound is in `results/Emax_bound_20C/`). Models are compared by chi-squared and AICc, with profile-likelihood intervals, and the result is tested against six assumptions about brain exposure (scenarios S1-S6).

## Contents

| Folder | Content |
|---|---|
| `data/` | Digitised group means and individual animals (CSV) and `README.md` with provenance and accuracy checks |
| `digitisation/` | Scripts that read the figures, the four source figure images of the first article (Figs 2-5, CC BY; the second article's figures are not included), and overlay images for visual checking |
| `analysis/` | Model fitting and figure scripts (`run_all.py` runs everything) |
| `results/` | Fitted PK parameters, scenario table, model comparison, parameter intervals, residuals |
| `figures/` | Figures 1-4 of the manuscript and the graphical-contents image |

## Reproducing

```bash
pip install -r requirements.txt
python analysis/run_all.py        # about 15 min; rewrites results/ and figures/
```

To redo the digitisation from the source figure images:

```bash
python digitisation/digi_fig2.py digitisation
python digitisation/digi_pd.py digitisation
python digitisation/digi_brain.py digitisation
python digitisation/digi_fig4b.py digitisation
python digitisation/assemble_data.py digitisation      # writes data/
```

## Second dataset: dose-ranging study at 100 ug/kg

The result is tested against a later study from the same laboratory: Canfield JR, Fry KK, Sprague JE. *Enhanced pharmacokinetic and pharmacodynamic effects of combining p-fluorofentanyl and fentanyl.* Drug Alcohol Depend 2025;272:112710, [doi:10.1016/j.drugalcdep.2025.112710](https://doi.org/10.1016/j.drugalcdep.2025.112710) (open access, **CC BY-NC-ND**). It gives plasma PK, brain/plasma ratios and effects for fentanyl and pFF at 100 ug/kg.
Only digitised **numbers** are included (`data/dose_ranging_*`); the article's figures are not redistributed. To redo the digitisation, download the article PDF and run `digitisation/extract_dose_ranging_figs.py`, then `digi_dad.py` and `digi_dad_brain.py`. Values that had to be read by eye where a marker was hidden are flagged in `data/dose_ranging_100ugkg_effect_observations.csv`.
`analysis/dose_test.py` predicts the 100 ug/kg time courses from the 300 ug/kg fit (out-of-sample) and fits both doses together; `results/dose_test_*` hold the output.

## Key result

With the measured brain/plasma ratios (pFF 2.3 times fentanyl), the shared-parameter model fits poorly (chi-squared 36.0 vs 23.6; dAICc 8.2) and pFF needs an EC50 1.85 times that of fentanyl (95% CI 1.47-2.24). The data identify an effective exposure ratio of about 1.2 (brain/plasma ratio divided by the EC50 ratio); how it splits between exposure and potency depends on the assumed brain exposure (Table 3 of the manuscript).

## Two-dose result

Fitted to both doses (40 observations) the shared-parameter model fits poorly (chi-squared 94.0 vs 69.8; dAICc 21.4) and the pFF EC50 ratio is 1.77 (95% CI 1.51-2.09), against 1.85 from the 300 ug/kg data alone. Tail-flick parameters are identified only with both doses.

## Limitations

Digitised group means only (no individual animals for plasma or effects); brain data are from separate animals at one time point; Emax for hypothermia is not identified; a few markers hidden behind others are less accurate (see `data/README.md`). Individual-animal data would allow a population analysis.

## Licences and citation

Code: MIT (`LICENSE`). Data: CC BY 4.0, attributed first to the source article (`LICENSE-DATA.md`). See `CITATION.cff`.
