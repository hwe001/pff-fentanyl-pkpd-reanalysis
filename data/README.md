# Digitised data from Canfield and Sprague (2025)

Source: Canfield JR, Sprague JE. *In vivo pharmacokinetic, pharmacodynamic and brain concentration comparison of fentanyl and para-fluorofentanyl in rats.* Arch Toxicol 2025;99:287-297. https://doi.org/10.1007/s00204-024-03887-z (open access, **CC BY**; PMC11748481). Figures 2-5 were taken from the publisher's open-access image files (copies in `source_figures_CC-BY/`) and digitised again at higher resolution. The original data are "available upon request" from the authors; none of the values here are unpublished.

Study design (from the paper): male Sprague-Dawley rats, 300 ug/kg sc; PK/PD groups n = 6 per drug and 5 saline controls, sampled at 0, 30, 60, 120, 240 and 480 min; **brain concentrations were measured at 60 min only, in a separate group of 12 rats (n = 6 per drug)**.

## Files

| File | Source figure | Content |
|---|---|---|
| `fig2_plasma_concentration.csv` | Fig. 2 | Mean plasma concentration, ng/mL, six times, both drugs |
| `fig3_tail_flick_MPE.csv` | Fig. 3 | Mean %MPE, three groups |
| `fig4A_temperature_change.csv` | Fig. 4A | Mean change in core temperature, degC |
| `fig4B_max_temperature_change_individual_animals.csv` | Fig. 4B | Maximum change in temperature for each of 6 animals per drug |
| `fig5_brain_individual_animals.csv`, `fig5_brain_group_means.csv` | Fig. 5A, B | Brain concentration (ng/g) and brain/plasma ratio at 60 min, four regions, individual animals and bar means |
| `table1_noncompartmental_PK_exact_from_source.csv` | Table 1 | Exact published noncompartmental values (typed from the article) |
| `overlay_check_*.png` | | Digitised points drawn over the figure, for visual checking |
| `digitiser_code/` | | The scripts that produced these files |

## Method

Axes were calibrated from the detected tick marks (residuals 0.0-0.26 min on time and 0.0-0.1 units on the y axes; 4 ng/g on the 4,000 ng/g brain axis). Markers were found by colour segmentation, filling the outlines of hollow symbols, and removing the connecting lines. Down-pointing triangles (pFF) sit about 16% of the marker height below the centre of their bounding box; this was calibrated on the exact Cmax and applied to all triangles.

## Accuracy checks against numbers the paper states exactly

| Check | Digitised | Published |
|---|---|---|
| Fentanyl Cmax (ng/mL) | 82.13 | 82.12 |
| pFF Cmax (ng/mL), before the triangle correction | 113.3 | 112.57 |
| Fentanyl SEM at Cmax, read from the error bar | 11.5 | 11.78 |
| Mean maximum temperature drop, fentanyl / pFF (degC) | -4.10 / -5.55 | abstract: -5.6 for pFF |
| Mean of the 6 individual points vs bar height (Figs 4B, 5) | within about 1% | n = 6 |
| Tail flick at 240 and 480 min vs the earlier CSVs | within 0.9 %MPE | |

## Known limitations

- Markers hidden behind another marker are marked in the `note` column. The pFF tail-flick points at 30, 60 and 120 min are hidden behind the fentanyl squares; both groups are at the 100% ceiling then.
- Temperature at 240 and 480 min (fentanyl, control) and plasma pFF at 240 and 480 min are partly overlapped: accuracy about +/-0.3 degC and +/-1 ng/mL.
- `error_bar_extent` and `sem_estimate` are the visible vertical extent of the error bar. Where the SEM is smaller than the marker, the value reflects the marker size, not the SEM, and should not be used as one.
- Fig. 5 pFF striatum (ng/g) has 5 visible points instead of 6 (two points overlap).
