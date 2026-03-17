# XRD Lab Technical Report Draft

Status: this draft already covers `1a`, `1b`, `2`, and `3c`, and it sets up placeholders for `3a`, `3b`, `3d`, and `4` once the actual lab data arrives. Replace any text in square brackets after the lab session.

## 1a. Final overall goal / research question

The overall analytical goal of this experiment is to determine the cubic unit-cell parameter of an unknown powder from its X-ray diffractogram and, by comparison with reference crystal structures, identify the most likely crystalline phase or phases present in the sample, including any impurity phase if visible.

If the full identification turns out to be uncertain, the intermediate goal is still certainly reachable: determine the best possible cubic lattice parameter `a` and assign Miller indices to the main peaks that belong to the dominant phase.

Citations: [S1], [S3]

## 1b. Own schematic drawing of Bragg's law

Figure: [bragg_law_schematic.svg](/Users/kirillrybkov/Desktop/xrd/assets/bragg_law_schematic.svg)

Suggested submission text:

Bragg diffraction can be understood by treating the crystal as a set of parallel atomic planes with spacing `d`. X-rays with wavelength `lambda` are scattered by atoms on successive planes, and the path difference between these scattered rays is `Delta = 2 d sin(theta)`. When this path difference equals an integer multiple of the wavelength, `n lambda`, the waves return in phase and constructive interference produces a diffraction peak. If the phase difference is not an integer multiple of the wavelength, the scattered waves partly or fully cancel and the corresponding peak is weak or absent.

The drawing also shows which factors affect peak intensity. The angle `theta`, wavelength `lambda`, and spacing `d` determine whether the phase relation is constructive at all, so they fix the peak position condition. The peak intensity is then governed mainly by how strongly the atoms scatter and how their arrangement inside the unit cell adds or cancels the scattered amplitudes, which is summarized by the structure factor. In addition, a larger number of well-ordered coherently diffracting planes gives a stronger and sharper peak than a poorly ordered or weakly crystalline material.

Citations: [S1], [S2], [S3]

## 2. Schematic drawing of the experimental setup

Figure: [powder_diffractometer_schematic.svg](/Users/kirillrybkov/Desktop/xrd/assets/powder_diffractometer_schematic.svg)

Suggested submission text:

The experiment is a powder X-ray diffraction measurement in Bragg-Brentano-type reflection geometry. The X-ray tube generates characteristic radiation from a metal anode; in the final report the anode material, wavelength, and tube power should be reported explicitly because they determine the diffraction condition and the achievable intensity. If a filter, mirror, or monochromator is used, that should also be stated because it affects background and whether unwanted spectral lines are suppressed.

The incident-beam optics, such as divergence slits, Soller slits, masks, or other conditioning optics, shape the beam before it reaches the sample. These settings matter because they control beam divergence, illuminated sample length, resolution, and counting statistics, so a reader needs them to reproduce the trade-off between intensity and angular resolution.

The sample itself is a powdered solid mounted in a holder, usually as a flat plate in reflection geometry. The mounting method and whether the sample was spun during the scan should be reported, because preferred orientation, surface flatness, and sample displacement can shift peak positions or distort relative intensities.

The goniometer controls the angular relation between source, sample, and detector, and in a standard coupled scan the diffracted beam is measured as a function of `2theta` while the sample angle changes accordingly. To reproduce the measurement, the final report should give the scan range, angular step size or scan speed, and counting time per step or equivalent dwell time.

On the diffracted-beam side, a receiving slit and the detector collect the scattered X-rays that satisfy the diffraction condition. The detector type should be reported, because point detectors, strip detectors, and area-derived detector modes differ in angular acceptance, count-rate behavior, and data-collection speed.

The acquisition computer stores intensity versus `2theta` and is used for background treatment, peak finding, and export of the diffractogram. If calibration standards or alignment checks were used, those should also be mentioned because they affect the confidence in the final peak positions and therefore in the derived lattice parameter.

Replace the bracketed items below after the lab:

- X-ray anode / wavelength: `[fill in after lab]`
- Geometry: `[fill in after lab]`
- Scan range in `2theta`: `[fill in after lab]`
- Step size or scan speed: `[fill in after lab]`
- Counting time / dwell time: `[fill in after lab]`
- Detector type: `[fill in after lab]`
- Slit settings / optics: `[fill in after lab]`

Citations: [S1], [S3], [S4]

## 3a. Experimental plot

Placeholder until lab data arrives.

Planned output:

- properly formatted intensity versus `2theta` plot
- axis labels with units
- all peaks visible
- logarithmic `y` axis if weak peaks or impurity peaks would otherwise disappear

## 3b. Calculation of the cubic unit-cell parameter and peak indexing

Placeholder until lab data arrives.

Planned content:

1. Show the measured peak positions in a table.
2. State the starting assumption for the first assigned peak, for example `h^2 + k^2 + l^2 = 1`, `2`, `3`, or `4`.
3. Calculate `sin^2(theta)` ratios and compare them with allowed cubic sequences.
4. Assign `(hkl)` values to the main phase peaks.
5. Mark impurity peaks explicitly if they do not fit the dominant cubic sequence.
6. Calculate `a = d * sqrt(h^2 + k^2 + l^2)` for each main-phase peak and report mean, lowest, and highest `a`.

## 3c. One precision-improving choice in the data analysis

Suggested submission text:

To obtain the most precise possible value of the cubic lattice parameter `a`, I would not base the final result on only one low-angle peak. Instead, I would first identify the peaks that belong to the dominant phase and then calculate `a` from several non-overlapping peaks, preferably including higher-angle peaks. The reason is that for a fixed angular uncertainty, the relative uncertainty in `d` and therefore in `a` becomes smaller at higher diffraction angles, so high-angle peaks usually constrain the lattice parameter more accurately than the first low-angle peak. I would also avoid weak impurity peaks or strongly overlapping peaks, because errors in peak position assignment would propagate directly into the calculated `a` value.

Citations: [S1], [S3]

## 3d. Simulation of two candidate crystal models

Placeholder until lab data arrives.

Planned content once data is available:

1. Use the measured `a` value to search the database for cubic structures with three or fewer elements.
2. Pick two plausible candidates and cite the database entries.
3. Simulate both powder patterns.
4. Compare simulated and experimental peak positions and relative intensities in one combined plot.
5. Discuss whether mismatches are better explained by a wrong candidate, impurity peaks, missing weak peaks, or an incorrect starting indexing assumption.

## 4. Discussion and conclusion

Placeholder until lab data arrives.

Planned structure:

- answer the research question from `1a` directly
- state the best candidate phase
- state how certain that identification is
- mention the strongest evidence and the strongest remaining doubt

## Sources

- [S1] Arnoud Onnink, *Intro to the X-ray diffraction lab*, course slides, 12 March 2026. File provided by user: [2026-03-12 Intro XRD.pdf](/Users/kirillrybkov/Downloads/2026-03-12%20Intro%20XRD.pdf)
- [S2] W. H. Bragg, *The Reflection of X-Rays by Crystals*, Nature 91, 477 (10 June 1913). DOI: [10.1038/091477b0](https://doi.org/10.1038/091477b0)
- [S3] Rigaku, *X-ray diffraction (XRD)*, official resource page: [rigaku.com/resources/techniques/x-ray-diffraction-xrd](https://rigaku.com/resources/techniques/x-ray-diffraction-xrd)
- [S4] Malvern Panalytical, *XRD analysis*, official overview page: [malvernpanalytical.com/en/learn/knowledge-center/analysis-methods-and-applications/xrd-analysis](https://www.malvernpanalytical.com/en/learn/knowledge-center/analysis-methods-and-applications/xrd-analysis)
