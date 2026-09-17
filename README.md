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

## Test

```bash
cd tests && matlab -batch "run_tests"
```

`test_sign` pins the dip sign convention against synthetics, `test_regrid`
covers the unit bug, and `test_opr_units` drives the whole OPR path over a
synthetic echogram with a known dip - including the all-NaN `Bottom` case.

## Validation

`tests/` covers the sign convention, the `regrid` unit bug, the whole OPR
path over a synthetic echogram, and sub-degree recovery under vertical
exaggeration.

`test_holschuh_regime` builds a synthetic test case and runs it with Nick's own 
published defaults. It recovers the dip with an **RMS error of 0.07 deg** and correlation
**r = 1.000** against truth. 

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
src/      Holschuh functions, with small fixes
opr/      OPR/CReSIS front end and the standard figure
examples/ runnable drivers, one per season
tests/    regression suite, plus diag_layering.m for siting the depth window
docs/     attribution and the upstream README
```
