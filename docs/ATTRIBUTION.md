# Attribution and provenance

## The method

The rolling Radon transform approach to englacial layer slopes is Nick
Holschuh's:

> N. Holschuh, B. R. Parizek, R. B. Alley, S. Anandakrishnan (2017),
> *Decoding ice sheet behavior using englacial layer slopes*,
> Geophysical Research Letters 44, 5561-5570.
> <https://doi.org/10.1002/2017GL073417>

**Please cite that paper.** The idea implemented here is his: project a
window of the radargram along a family of directions, and take the
direction whose projection is most coherent as the layer orientation.

Reference implementations, worth reading:

- <https://github.com/nholschuh/SlopeExtraction_Radar>
- <https://github.com/nholschuh/NDH_MatlabTools>

## The code

**This repository contains no code from either of those repositories.**
Everything in `src/` and `opr/` was written here. There is no dependency on
`SlopeExtraction_Radar` or `NDH_MatlabTools`, and neither needs to be on the
MATLAB path.

That was a deliberate decision, for three reasons:

1. `SlopeExtraction_Radar` carries no licence file, so no terms of
   redistribution attach to it. Writing our own avoids the question.
2. Its published release and `NDH_MatlabTools` disagree with each other on
   the SIGN of the result: the public release returns rising-positive, while
   the full toolbox applies `slopegrid*-1` at the end. Inheriting from
   either invites the sign to drift. Ours is stated, tested and fixed
   (see below).
3. Several implementation choices needed to change for this application -
   a variance criterion rather than a peak amplitude, a rectangular window,
   an explicit angular step - and it is cleaner to own the estimator than to
   carry a patched copy of someone else's.

`test_holschuh_regime.m` checks that the reimplementation still reproduces
the method in the regime the paper targeted: RDS/impulse radar, 2.8 m range
resolution, layers 45 m apart, folds giving +/-15 degree slopes. It recovers
the slope with an RMS error of 0.04 degrees and r = 1.000 against truth.

## Sign convention

    POSITIVE slope = the layer RISES (gets shallower) with increasing x.

This is `d(elevation)/dx`, the standard glaciological sense, and it matches
the public `SlopeExtraction_Radar` release. `NDH_MatlabTools` produces the
opposite. `tests/test_sign.m` pins ours down against synthetics so it cannot
drift; check which one your downstream analysis expects before comparing
numbers between codebases.

## Files

| Path | Origin |
|---|---|
| `src/ls_radon_dip.m` | Written here. Radon slope estimator. |
| `src/ls_rolling_radon.m` | Written here. Rolling-window driver. |
| `src/ls_polarstereo_fwd.m` | Written here, from Snyder (1987). |
| `src/ls_cice.m`, `src/ls_slope_colormap.m` | Written here. |
| `opr/*` | Written here. OPR/CReSIS front end, figures, diagnostics. |
| `examples/*`, `tests/*` | Written here. |

The only external dependencies are MATLAB's own `radon` (Image Processing
Toolbox) and `prctile` (Statistics Toolbox).
