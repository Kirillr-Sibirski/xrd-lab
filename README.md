# XRD Lab Project

This repository contains the code, data, figures, and report notes used for a cubic powder-XRD lab analysis.

## Lab Instructions Summary

The lab task was to analyze an unknown powder diffractogram under the assumption that the unknown is an inorganic compound with a cubic unit cell. The main required deliverables were:

- plot the experimental powder XRD pattern with all visible peaks
- calculate the cubic unit-cell parameter `a`
- assign `hkl` indices to the observed reflections
- explain one choice that improves the precision of the `a` value
- simulate two candidate crystal models from an online database
- compare the simulations with the experimental pattern
- conclude with the most likely crystal identity and discuss uncertainty

This repository includes the data analysis, figures, and draft text used to answer those tasks.

## What We Did

We built a small Python workflow to:

- read powder-XRD data from `.int` and two-column `.xy` files
- detect the main diffraction peaks
- calculate `d` spacings from Bragg's law
- compare `sin^2(theta)` ratios against allowed cubic sequences
- estimate the cubic unit-cell parameter `a`
- simulate powder diffractograms from AMCSD CIF files
- compare candidate crystal models with the experimental pattern

## Tasks Completed In This Repo

- prepared experimental plots for the Group 3 dataset
- checked weak-feature regions near `33.32°` and `59.16°`
- indexed the dominant reflections with an FCC-like cubic sequence
- estimated the lattice parameter from the strongest reflections
- compared FCC and BCC indexing logic
- simulated candidate diffractograms for `ZnS` and `CuCl`
- checked pyrite as a competing candidate and compared its extra predicted peaks
- drafted report-ready wording for the main lab questions

## Sample Analyzed Here

The included dataset is:

- `data/gr3.xy`

For the Group 3 sample, the main experimental reflections were identified near:

- `28.4794°`
- `47.3083°`
- `56.1585°`
- `69.1674°`

We also observed very weak additional features near:

- `33.3245°`
- `59.1649°`

Using the strong reflections, we obtained:

- `a = 5.4276 Å = 0.54276 nm`

The dominant pattern is most consistent with an FCC, zinc-blende-like diffraction pattern. From the simulated candidate comparisons, `ZnS` was the best current match and `CuCl` remained a plausible alternative.

## Key Files

- `xrd_analyzer.py`
  Main CLI for peak detection, cubic indexing, lattice-parameter estimation, and candidate ranking.
- `cif_matcher.py`
  CIF parsing, symmetry expansion, powder-pattern simulation, and candidate scoring.
- `amcsd_cache.py`
  Helpers for downloading and indexing the AMCSD CIF archive locally.
- `data/gr3.xy`
  Experimental Group 3 XRD data used in this analysis.
- `gr3_xrd_plot.png`
  Main experimental plot.
- `gr3_raw_full_linear.png`
  Raw-data full-range plot on linear scale.
- `gr3_raw_full_log.png`
  Raw-data full-range plot on log scale.
- `gr3_weak_peaks_simple.png`
  Simple zoomed plot showing the weak features near `33.3°` and `59.2°`.
- `gr3_vs_zns.png`
  Experimental vs simulated comparison for `ZnS`.
- `gr3_vs_cucl.png`
  Experimental vs simulated comparison for `CuCl`.
- `xrd_lab_submission_draft.md`
  Draft report text and supporting notes.
- `docs/lab_task_summary.md`
  Short summary of the lab instructions and the work completed here.

## Reproducing The Analysis

Run the analyzer on the included dataset:

```bash
python3 xrd_analyzer.py data/gr3.xy
```

JSON output:

```bash
python3 xrd_analyzer.py data/gr3.xy --json
```

To scan the full filtered AMCSD candidate set:

```bash
python3 xrd_analyzer.py data/gr3.xy --cache-dir . --amcsd-limit 0 --top-candidates 10
```

## Notes On Candidate Selection

The AMCSD search returned several candidates with lattice parameters close to the measured `a`, including silicon, nitrogen, `NpO2`, pyrite, and `AlP`. We did not choose candidates by `a` alone. We also filtered by plausibility and by the full diffraction pattern:

- silicon was excluded because it is a pure element rather than the unknown inorganic compound described in the lab slides
- nitrogen was excluded because it is not a solid powder at room temperature
- `NpO2` was excluded because neptunium compounds are radioactive
- pyrite was rejected as the dominant phase because it predicts extra reflections that are not clearly present in the measured pattern
- `AlP` was treated as less likely because it can be hazardous

That left `ZnS` and `CuCl` as the most reasonable candidates to compare in detail.

## Weak Peaks

The weak features near `33.32°` and `59.16°` are only marginal in the raw data. In the report, they should be described as weak or tentative features rather than strong resolved peaks. The plot `gr3_weak_peaks_simple.png` was made specifically to show these regions more clearly without smoothing away the raw signal.

## Tests

```bash
python3 -m unittest discover -s tests
```

## Large Local Cache Files

This repository does not commit the full local AMCSD cache or extracted CIF tree, because those files are large and can be regenerated locally by the scripts if needed.
