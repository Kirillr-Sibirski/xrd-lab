# XRD Lab Task Summary

## Assignment Summary

The lab report was built around the analysis of an unknown powder sample measured by X-ray diffraction. The key assumptions and tasks were:

- the unknown is an inorganic compound with a cubic unit cell
- plot the experimental XRD pattern clearly
- calculate the cubic unit-cell parameter `a`
- assign `hkl` values to the observed peaks
- explain one precision-improving analysis choice
- simulate two candidate crystal models from an online crystal database
- compare the simulated and experimental diffractograms
- conclude with the most likely crystal identity and explain the uncertainty

## What Was Done Here

This repository contains the analysis workflow and supporting files used for that assignment:

- experimental data copied into `data/gr3.xy`
- analyzer code for cubic indexing and lattice-parameter estimation
- local CIF-based simulation code for candidate comparison
- plots of the full experimental pattern and zoomed weak-feature regions
- simulated comparisons for `ZnS` and `CuCl`
- pyrite-check zooms to inspect competing reflections
- draft text for the report questions, including the discussion and conclusion

## Main Technical Outcome

For the Group 3 dataset, the dominant reflections were identified near `28.48°`, `47.31°`, `56.16°`, and `69.17°`, with tentative weak features near `33.32°` and `59.16°`. Using the strongest reflections, the best-fit cubic lattice parameter was:

- `a = 5.4276 Å = 0.54276 nm`

The pattern was most consistent with an FCC, zinc-blende-like diffraction pattern. Among the candidates checked in detail, `ZnS` was the best current match and `CuCl` remained a plausible alternative.
