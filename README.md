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

## What this adds

The published `SlopeExtraction_Radar` code does not run against a current OPR
product, and its CReSIS driver (`RollingRadon_CReSIS.m`, in NDH_MatlabTools)
cannot run at all. This repo carries that code with the defects fixed, plus an
OPR front end and a regression suite.

### Bugs fixed in the original code

| Where | Defect | Consequence |
|---|---|---|
| `radon_ndh.m` | Returned dip with the **sign inverted** relative to any stated convention. Nick's full NDH_MatlabTools `RollingRadon` compensates with `slopegrid*-1` at the very end; the public release dropped that line. | Every slope from the public repo points the wrong way. For migration work this reverses the inferred direction. |
| `regrid.m` | Target spacing computed as `(1/f)/20`, a **time**, then applied to the distance axis. | For a 2.5 km line at 600 MHz it asks for ~3×10¹³ samples and dies in the allocator. |
| `radon_ndh.m` | `data = (data - min(data))/max(data)` divides by the max of the *unshifted* data, which is negative for dB input. | Log-power windows came back sign-flipped; methods 3 and 4 are not invariant to that. |
| `radon_ndh.m` | `find(x == max(x))` returns several indices on a tie and none at all for a featureless window. | Tie → size error on assignment; empty → `Index exceeds array bounds`. |
| `radon_ndh.m` | `griddedInterpolant` requires ascending sample points. | A descending (elevation) y axis threw before reaching the branch meant to handle it. |
| `RollingRadon.m` | Microsecond detection tested the **sample step** (`Time(2)-Time(1) > 1e-6`). | A finely sampled microsecond record passes; a coarse second-valued one fails. Now tests the record length. |
| `RollingRadon.m` | `lp(Data)` is `20*log10`, but OPR products store **detected power**. | Doubled the dB scale, silently halving the effect of every dB threshold. |
| `RollingRadon.m` | `breaks = [1 length(Data(1,:))]` with chunks indexed `breaks(k):breaks(k+1)-1`. | Dropped the final trace of every line short enough for one chunk. |
| `RollingRadon.m` | `if exist('movie')` tested the wrong variable name (the argument is `movie_flag`). | `movie_flag` was unconditionally reset to 0; the animation option could never be enabled. |
| `RollingRadon.m` | `exist('max_frequency') == 1` is true for a variable that exists but is **empty**. | Passing `[]` as a positional placeholder selected the regridding branch and divided by `[]`. |
| `RollingRadon.m` | Window forced square in samples. | Along-track extent needs enough traces to define a dip; vertical extent needs to stay short enough to avoid system power drift. Now accepts `[horizontal vertical]`. |
| `RollingRadon.m` | Continuity filter **substitutes** `last_val` for a deviating window and carries it forward. | Paints long constant-dip columns of fabricated values through the field. Default is now to keep the measurement (`vr = 1`). |
| `radon_ndh.m` | Rebuilt the Radon correction images and Gaussian filters on **every call**, though they depend only on window size. | Three extra Radon transforms per window. Now cached; the Ridge A frame went from minutes to ~35 s at a 10× finer grid. |
| `radon_ndh.m` | Angular search step hard-coded at **0.1 deg**. | Interior layers dip ~0.1 deg, so the whole signal fell inside one search increment and every window returned 0 or +/-0.1. Now an argument (`dip_step`). |
| `RollingRadon.m` | Chunk width fixed at 1000 columns regardless of the window. | A wider window makes `roll_steps` negative, so the loop never runs, `opt_x` is never created, and the function dies in its own epilogue having silently processed nothing. Long windows are exactly what sub-degree dips need. |
| `regrid.m` | Calls `cice_import`, which the public repo does not ship. | A clean checkout cannot regrid radar data. Included here. |

`RollingRadon_CReSIS.m` is **not** carried over. It opened `load(filename);
clearvars ...; save(filename)` — writing its own arguments back into the shared
data product — left `steps` undefined for any line under 4000 traces,
re-processed the whole image on every chunk iteration, indexed `opt_angle` as a
scalar and then sub-indexed the result, hard-coded Windows paths, and called
`RadialSpreading`, which is not published in any of Nick's repositories.
`opr/RollingRadon_OPR.m` replaces it.

> **Heads up:** `Data_20260109_02_003.mat` in
> `2025_Antarctica_Ground2/CSARP_standard_HH` already carries `filename`,
> `movie`, `plotter` and `window` from a previous run of that code. The
> echogram itself is intact, but the file is no longer a pristine OPR
> product. `opr_load_echogram` warns when it sees this.

### What the OPR front end does

- **Reads strictly read-only.** Results always go to a separate file.
- **Finds the bed.** The standard product ships an all-NaN `Bottom` for this
  season; the bed is pulled from `CSARP_layer` (layer id 2) and interpolated
  onto the frame's `GPS_time`.
- **Crops to the ice column.** The example frame carries 32066 samples over
  ~4.5 km of range for ~1 km of ice.
- **Grids isotropically, once.** Radon needs square pixels. Native sampling is
  0.14 m vertically against 5.9 m along track, so `radon_ndh` would otherwise
  upsample every window ~42×.
- **Conditions for layering, not gain.** A depth low-pass (`smooth_len`)
  suppresses structure finer than the layering; per-trace balancing removes
  profile-to-profile gain swings that show up as vertical striping and that
  the Radon will otherwise lock onto. An optional depth high-pass
  (`detrend_len`, off by default) removes the long-wavelength power envelope
  when it dominates a window.
- **Exaggerates vertically.** Interior layers dip ~0.1 deg, which is below
  one range cell of displacement across any tractable window. See
  `vert_exag` under *Choosing the parameters*.
- **Gates honestly.** The whole window must sit inside the ice column, not
  just its centre, and rejection reasons are returned in `R.status`.

---

## Requirements

MATLAB with the Image Processing (`radon`), Statistics (`normpdf`) and Signal
Processing toolboxes. Tested on R2024b (`/opt/sw/matlab/2024b/bin/matlab` on
the CReSIS machines). Nothing else - the Holschuh helpers are vendored in
`src/`, including `cice_import.m`, which upstream does not ship.

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

Then run a profile. The examples resolve their own paths, so they work from
any checkout location:

```bash
cd ../examples
nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \\
    "maxNumCompThreads(8); run_slopes_2024_Antarctica_Ground2"
```

That processes the 20 km mega-dune profile `20250108_02_005` from
`CSARP_post/CSARP_standard` - the same posted product `imb.picker` displays -
and writes a `.mat` and a `.png` under
`/kucresis/scratch/$USER/layer_slopes/products/`. Frame, product and output
location can all be overridden on the command line:

```bash
matlab -batch "frame='20250108_02_001'; run_slopes_2024_Antarctica_Ground2"
matlab -batch "out_dir='/kucresis/scratch/$USER/slopes'; run_slopes_2024_Antarctica_Ground2"
matlab -batch "run_slopes_2025_Antarctica_Ground2"    # Ridge A, 2025 season
```

Long lines are worth detaching, since nothing needs a terminal:

```bash
nohup nice -n 10 /opt/sw/matlab/2024b/bin/matlab -batch \\
    "run_slopes_2024_Antarctica_Ground2" > run.log 2>&1 &
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

The two that matter most, and why:

- **`grid_spacing`** must sample the system's range resolution — 0.53 m in ice
  for this accumulation radar (~159 MHz bandwidth). At 1–2 m the anti-alias
  average removes the layering before the Radon transform ever sees it, which
  shows up as low coverage and implausibly large dips.
- **`window_z`** must be much larger than the range resolution so there is
  layering inside it to measure, but small enough that system power drift is
  not the strongest linear feature. ~30 m (≈57 range cells) works here.
  `window_x` is independent: 120 m is ~20 traces at 5.9 m spacing.

`R.status` reports why windows were rejected (1 outside ice, 2 low SNR,
3 slope gate), which is the fastest way to tell whether a setting is too
strict.

## Test

```bash
cd tests && matlab -batch "run_tests"
```

`test_sign` pins the dip sign convention against synthetics, `test_regrid`
covers the unit bug, and `test_opr_units` drives the whole OPR path over a
synthetic echogram with a known dip — including the all-NaN `Bottom` case.

## Validation

`tests/` covers the sign convention, the `regrid` unit bug, the whole OPR
path over a synthetic echogram, and sub-degree recovery under vertical
exaggeration.

The check that matters is against a real picked profile. On
`20250108_02_005` the horizon visible in `imb.picker` runs from ~1.85 us at
0 km to ~1.45 us at 20 km - about 35 m of relief over 20,000 m, a dip of
**-0.10 deg**. The solver returns a median of **-0.09 deg** (IQR -0.20 to
+0.01) over 443 solved windows, and the negative sign correctly says the
layers shallow with increasing distance.

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
