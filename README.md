# layer_slopes

Englacial layer slopes from Open Polar Radar / CReSIS echograms, by rolling
Radon transform.

The eventual target is mapping **mega-dune** position and migration through a
survey domain: mega-dunes leave a characteristic dipping signature in the
englacial stratigraphy, so a reliable layer-slope field across many profiles
is the raw material for locating them and tracking how they move.

The Radon method, and most of the code under `src/`, is Nick Holschuh's:

> N. Holschuh, B. R. Parizek, R. B. Alley, S. Anandakrishnan (2017),
> *Decoding ice sheet behavior using englacial layer slopes*,
> Geophysical Research Letters 44, 5561-5570.

- <https://github.com/nholschuh/SlopeExtraction_Radar>
- <https://github.com/nholschuh/NDH_MatlabTools>

See `docs/ATTRIBUTION.md` for exactly what came from where and what changed.

---


## What the OPR front end does

- **Reads strictly read-only.** Results always go to a separate file.
- **Finds the bed.** The standard product ships an all-NaN `Bottom` for this
  season; the bed is pulled from `CSARP_layer` (layer id 2) and interpolated
  onto the frame's `GPS_time`.
- **Crops to the ice column.** The example frame carries 32066 samples over
  ~4.5 km of range for ~1 km of ice.
- **Grids isotropically, once.** Radon needs square pixels. Native sampling is
  0.14 m vertically against 5.9 m along track, so `radon_ndh` would otherwise
  upsample every window ~42x.
- **Conditions for layering, not gain.** A depth low-pass (`smooth_len`)
  suppresses structure finer than the layering; per-trace balancing removes
  profile-to-profile gain swings that show up as vertical striping and that
  the Radon will otherwise lock onto. An optional depth high-pass
  (`detrend_len`, off by default) removes the long-wavelength power envelope
  when it dominates a window.
- **Exaggerates vertically.** Interior layers dip ~0.1 deg, which is below
  one range cell of displacement across any tractable window. See
  `vert_exag` under *Choosing the parameters*.
- **Removes the vertical striping.** Trace-to-trace gain changes draw
  vertical stripes, and a stripe is a strong linear feature the Radon will
  fit. `smooth_x` averages them away, and at these dips it is nearly free: a
  0.1 deg layer moves 0.4 m across 250 m of track, inside one range cell.
  This is the single biggest lever on how coherent the field looks -
  along-track continuity improves from 0.130 to 0.010 deg between adjacent
  cells, and sign consistency from 0.75 to 0.89.
- **Gates honestly.** The whole window must sit inside the ice column, not
  just its centre, and rejection reasons are returned in `R.status`.
- **Masks the merged pulse return.** An accumulation radar's first and
  second pulse returns merge at a fixed range (~75 m here), and that band is
  a strong horizontal feature with no stratigraphic meaning. `exclude_z`
  blanks it, so windows overlapping it abstain while the layering above and
  below is still solved. The shallow section above the band carries some of
  the clearest layering in the profile.
- **Plots a continuous raster.** Overlapping windows drawn as raw cells make
  a staircase of rectangles whose edges come from the window spacing, not
  the ice. `plot_slope_field` interpolates onto a fine grid for the smooth
  slope raster of Holschuh et al. (2017, fig. 3).

---

## Requirements

MATLAB with the Image Processing Toolbox (`radon`) and Statistics Toolbox
(`prctile`). Tested on R2024b (`/opt/sw/matlab/2024b/bin/matlab` on the
CReSIS machines).

**No other dependency.** This repository contains no code from
`SlopeExtraction_Radar` or `NDH_MatlabTools` - neither needs to be on the
path. The method is Holschuh's and is cited; the implementation is our own.
See `docs/ATTRIBUTION.md`.

## Install and run on the CReSIS servers

Clone into your own scratch space and run in place. Echograms are opened
read-only, and nothing is written outside the output directory you choose.

```bash
ssh mem1                                    # or any CReSIS compute node
cd /kucresis/scratch/$USER/scripts          # your scratch, not someone else's
git clone https://github.com/hoffmaao/layer_slopes.git
cd layer_slopes
```

Check the install before pointing it at real data. The suite needs no
external files and runs in a couple of minutes:

```bash
cd tests
/opt/sw/matlab/2024b/bin/matlab -batch "run_tests"
```

Then run a profile. The scripts resolve their own paths, so they work from
any checkout location:

```bash
cd ../examples
nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \\
    "maxNumCompThreads(8); run_slopes; make_slope_figure"
```

There are three scripts and one config:

| file | does |
|---|---|
| `ls_config.m` | the target frame and every solver setting, in one place |
| `run_slopes.m` | computes the slope field, saves a `.mat` |
| `make_slope_figure.m` | slope field over the power image, saves a `.png` |
| `make_window_movie.m` | animated GIF of the windows the Radon is fed |

All three read `ls_config.m`, so they cannot drift apart. The default target
is the 20 km mega-dune profile `20250108_02_005` from
`CSARP_post/CSARP_standard_HH` - the SAR-focused product, so the layers are
migrated. Output goes under
`/kucresis/scratch/$USER/layer_slopes/products/`.

`make_window_movie` steps the solver's own window in raster order, left to
right along a row then down to the next, and shows for each one the window
exactly as the estimator receives it, the fitted slope drawn on it, and the
criterion against candidate slope. It is the tool for answering "why did it
return that?" - a flat criterion means the window had nothing to lock onto
and should be abstaining.

The frame can be overridden from the command line; anything else is edited
in `ls_config.m`:

```bash
matlab -batch "frame='20250108_02_001'; run_slopes"
matlab -batch "make_window_movie"
```

Long lines are worth detaching, since nothing needs a terminal:

```bash
nohup nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \\
    "run_slopes" > run.log 2>&1 &
```

`maxNumCompThreads(8)` and `nice` keep a shared node usable; neither is
required for correctness.

Directly:

```matlab
addpath src; addpath opr;
R = RollingRadon_OPR(data_file, ...
        'grid_spacing', 0.25, ...   % m; must resolve the range resolution
        'window_x', 120, ...        % m along track
        'window_z', 30, ...         % m vertical
        'dip_max', 12, ...
        'detrend_len', 10, ...
        'out_file', 'slopes.mat');
plot_slope_field(R, 'slopes.png');
```

### Choosing the parameters

Interior englacial layers dip by **~0.1 degrees**. Almost every setting
follows from that number, and getting one wrong produces a slope field that
looks plausible and means nothing.

- **`vert_exag`** is what makes the measurement possible at all. A 0.1 deg dip
  displaces a layer by 0.2 m across a 120 m window, less than the 0.53 m range
  resolution, so it simply is not in the data. Sampling along track at
  `vert_exag` times the vertical spacing and presenting the grid to the Radon
  as isotropic multiplies the apparent dip by exactly `vert_exag`
  (`tan(apparent) = vert_exag * tan(true)`), and the driver inverts it on
  output. At 20x a 0.1 deg dip presents as 2 deg, for 20x fewer pixels per
  window. `tests/test_vert_exag.m` shows the plain grid returning a single
  quantised value where the exaggerated grid resolves thirteen.
- **`window_x`** must be long enough that the dip displaces a layer by more
  than one range cell. The solver prints that floor on every run
  (`one range cell across the window = 0.030 deg` for a 1 km window); keep the
  expected dip well above it.
- **`dip_step`** was effectively 0.1 deg upstream, so the entire signal fit
  inside a single search increment.
- **`grid_spacing`** must sample the system's range resolution, 0.53 m in ice
  for this radar (~159 MHz bandwidth). At 1-2 m the anti-alias average removes
  the layering before the Radon sees it: low coverage, implausibly large dips.
- **`window_z`** wants many range cells, but short enough that system power
  drift is not the strongest linear feature in the window.
- **`z_max` / `z_pad_surface`** decide whether you measure stratigraphy or
  speckle. On `20250104_01_002` the 20-110 m band has 5.9 dB of band-passed
  contrast and discrete spectral peaks at 1-11 m; the 120-600 m band has
  1.5 dB, sits 5 dB above the noise floor at -122 dB, and has a FLAT spectrum
  from 1 to 50 m, which is white noise. Run `tests/diag_layering.m` on a new
  site before choosing the depth window.

`R.status` reports why windows were rejected (1 outside ice, 2 low SNR,
3 slope gate), which is the fastest way to tell whether a setting is too
strict.

## Test

```bash
cd tests && matlab -batch "run_tests"
```

`test_sign` pins the dip sign convention against synthetics, `test_regrid`
covers the unit bug, and `test_opr_units` drives the whole OPR path over a
synthetic echogram with a known dip - including the all-NaN `Bottom` case.

## Choosing the window

The two window dimensions have independent limits and pull opposite ways.
Measured on `20250108_02_005` against an independently tracked horizon
(`tests/validate_horizon.m`), where `gain` is the regression slope of solver
against truth and 1.00 would be unbiased:

| `window_x` | smallest measurable slope | gain | r |
|---|---|---|---|
| 750 m | 0.040 deg | 0.99 | 0.98 |
| **1000 m** | **0.030 deg** | **0.96** | **0.98** |
| 2000 m | 0.015 deg | 0.85 | 0.97 |
| 3000 m | 0.010 deg | 0.71 | 0.91 |

A longer window sees smaller slopes but averages over a range of true ones
and regresses toward their mean. 1000 m balances the two on this data.

`window_z` behaves differently: it is cheap, because the system resolves
0.53 m and layers sit ~8 m apart, so a 20 m window already holds several
cycles. Sampling MORE layers actively hurts - gain falls from 0.86 at 20 m
to 0.74 at 100 m, because dip varies with depth and a tall window averages
across it.

Combining several window sizes was tried and **did not help**. Merging
scales on mutual agreement gave gain 0.86 and r 0.91 against 0.96 and 0.98
for the best single scale. The reason is worth recording: the coarse scale
is systematically biased low, so accepting cells where two scales agree
preferentially selects cells that agree with that bias. Corroboration
between estimators with different biases selects for the bias, not for
truth. One calibrated window is better.


## Validation

`tests/` covers the sign convention, the `regrid` unit bug, the whole OPR
path over a synthetic echogram, and sub-degree recovery under vertical
exaggeration.

**The method, as published, still works.** `test_holschuh_regime` builds a
synthetic in the regime Holschuh et al. (2017) targeted - RDS/impulse radar,
~2.8 m range resolution, layers 45 m apart, folds giving +/-15 deg reflector
slopes - and runs it with Nick's own published defaults and every addition
here switched off (no vertical exaggeration, no along-track smoothing, no
trace balancing, `o_f` 2/6, `snr_thresh` 2, `vr` 3, Radon method 0). It
recovers the dip with an **RMS error of 0.07 deg** and correlation
**r = 1.000** against truth, at 72% coverage. Nothing in this repo has
broken the method; the difficulty on accumulation-radar data is that 0.1 deg
dips in the top 200 m of a 0.53 m-resolution image is a different and much
harder problem than the one it was built for.

The check that matters is against real data, and it does not use the Radon
at all. `tests/validate_horizon.m` seeds on the bright reflector, follows it
trace by trace under a continuity constraint, writes the pick out as a
figure so it can be confirmed by eye, then differentiates it. On
`20250108_02_005` that horizon rises from 158 m to 121 m over 20 km.

Against that truth the calibrated settings give a **regression gain of
0.96** and **r = 0.98** - the magnitude is right, not just the sign. Two
things were needed to get there, and both were found this way rather than
by guessing:

- Without the depth high-pass the gain is 0.70 and r is 0.78, because a
  window spans part of the power-vs-depth decay and that gradient is
  horizontal.
- A 2000 m window gives gain 0.86 and a 3000 m window 0.71, because a long
  window averages across varying dip.

## Output

`RollingRadon_OPR` returns, and saves, a struct:

| field | meaning |
|---|---|
| `slope_x` | along-track distance of each window centre (m) |
| `slope_z` | depth below the ice surface of each window centre (m) |
| `slopes` | layer dip (deg). **Positive = the layer deepens with increasing x.** |
| `status` | 0 solved, 1 outside ice, 2 low SNR, 3 slope gate |
| `lat`, `lon`, `x`, `y` | geolocation of each slope column |
| `bed_z` | bed depth below surface at each column (m) |
| `param` | every setting, plus bed/surface provenance and runtime |

## Layout

```
src/      Holschuh's functions, with the fixes above
opr/      OPR/CReSIS front end and the standard figure
examples/ runnable drivers, one per season
tests/    regression suite, plus diag_layering.m for siting the depth window
docs/     attribution and the upstream README
```
