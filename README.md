# layer_slopes

Englacial layer slopes from Open Polar Radar / CReSIS echograms, by rolling
Radon transform.

The eventual target is mapping **mega-dune** position and migration through a
survey domain: mega-dunes leave a characteristic dipping signature in the
englacial stratigraphy, so a reliable layer-slope field across many profiles
is the raw material for locating them and tracking how they move.

The rolling-Radon method is Nick Holschuh's:

> N. Holschuh, B. R. Parizek, R. B. Alley, S. Anandakrishnan (2017),
> *Decoding ice sheet behavior using englacial layer slopes*,
> Geophysical Research Letters 44, 5561-5570.

- <https://github.com/nholschuh/SlopeExtraction_Radar>
- <https://github.com/nholschuh/NDH_MatlabTools>

See `docs/ATTRIBUTION.md` for exactly what came from where and what changed.

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
        'grid_spacing', 0.25, ...   % m; must resolve the 0.53 m range cell
        'window_x', 1000, ...       % m along track
        'window_z', 20, ...         % m vertical
        'dip_max', 1, ...           % deg; interior layers dip ~0.1 deg
        'exclude_z', [70 88], ...   % merged pulse return
        'out_file', 'slopes.mat');
plot_slope_field(R, 'slopes.png', 'clim_dip', [-0.3 0.3]);
```

## Test

```bash
cd tests && matlab -batch "run_tests"
```

`test_sign` pins the dip sign convention against synthetics,
`test_opr_units` drives the whole OPR path over a synthetic echogram with a
known dip - including the all-NaN `Bottom` case - `test_vert_exag` covers
sub-degree recovery, `test_holschuh_regime` runs Nick's published settings
in the regime they were built for, `test_null_gate` checks that the quality
gate holds pure speckle to the requested false-alarm rate while admitting
layering, and `test_search_edge` checks that a slope at the edge of the
search is withheld rather than reported.

## Validation

**Synthetic.** `test_holschuh_regime` builds an echogram in the regime
Holschuh et al. (2017) targeted and runs it with Nick's published settings:
RMS dip error 0.04 deg, r = 1.000 against truth.

**Real data.** `tests/validate_horizon.m` tracks a bright reflector trace
by trace, with no Radon involved, differentiates it, and compares every
window on the reflector with the reflector's own dip across that window.
It saves the pick over the echogram so it can be checked by eye first.
With the settings in `ls_config.m`:

| frame | window_x | windows on horizon | solved | gain | r |
|---|---|---|---|---|---|
| 20250108_02_005 | 500 m | 157 | 100% | 1.04 | 0.96 |
| 20250108_02_005 | 1000 m | 77 | 100% | 1.04 | 0.98 |
| 20250108_02_005 | 2000 m | 37 | 100% | 1.09 | 0.96 |
| 20250112_01_008 | 500 m | 157 | 100% | 1.01 | 0.97 |
| 20250112_01_008 | 1000 m | 77 | 100% | 1.02 | 0.995 |
| 20250112_01_008 | 2000 m | 37 | 100% | 1.04 | 0.995 |

Gain is the regression slope of the solver on the tracked dip, so 1 means
the magnitude is right as well as the sign.

```matlab
validate_horizon('20250112_01_008', 10000, 130)   % frame, seed x (m), seed z (m)
```

## Quality gate

A window is kept when its **semblance** - the fraction of its energy that
stacks coherently when every column is shifted along the fitted slope - is
higher than noise would reach. "Noise" is measured, not assumed: the same
echogram with every trace shifted in depth by up to +/-`null_jitter` (15 m,
about two layer spacings), put through the identical conditioning and
windows. The layers no longer line up, but each trace keeps its own
statistics and the power envelope stays where it was. The threshold is the
99th percentile of that noise at each depth (`false_alarm` 0.01), so a
horizontal feature that is not layering, such as the firn power envelope,
raises the bar where it occurs rather than passing as layering.

Semblance replaced the Radon peak/median ratio `q`. `q` measures how much
better the best slope is than the other candidate slopes, and when a window
cannot resolve slopes much finer than the search range, that ratio stays
small even over clean layering. On 20250112_01_008, `q` told real windows
from noise no better than chance at `window_x` 500 m (AUC 0.55), while
semblance did so at AUC 0.94-0.99 at every window size.

Two things to know when reading a field:

- A sharp horizontal amplitude boundary is genuinely coherent, so windows on
  it are kept at the boundary's own dip. On accumulation radar, watch the
  bright-to-dark transition around 150-200 m.
- Calibration scores the echogram twice, so a run takes about twice as long.
  Pass a fixed `semb_thresh` to skip it.

## Output

`RollingRadon_OPR` returns, and saves, a struct:

| field | meaning |
|---|---|
| `slope_x` | along-track distance of each window centre (m) |
| `slope_z` | depth below the ice surface of each window centre (m) |
| `slopes` | layer dip (deg). **Positive = the layer rises (gets shallower) with increasing x.** |
| `semb` | semblance along the fitted slope, the quality measure (0-1) |
| `semb_thresh` | the gate applied at each window row, calibrated against noise |
| `q` | Radon criterion peak/median, kept as a diagnostic |
| `status` | 0 solved, 1 outside ice, 2 below the semblance gate, 3 slope gate, 4 best slope at the edge of the search |
| `lat`, `lon`, `x`, `y` | geolocation of each slope column |
| `bed_z` | bed depth below surface at each column (m) |
| `param` | every setting, plus bed/surface provenance and runtime |

## Layout

```
src/      the solver: Radon slope estimator, rolling window, helpers
opr/      OPR/CReSIS front end, the standard figure, the window movie
examples/ ls_config and the drivers that read it
tests/    regression suite, validate_horizon.m, and diag_layering.m for
          siting the depth window
docs/     attribution and the upstream README
```
