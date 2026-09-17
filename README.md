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

`SlopeExtraction_Radar` was published in 2017 alongside the paper, as a
deliberate first step toward making the method usable by other people. Its
README says so plainly:

> I am still making sure this has all the dependencies required, but this is
> my first step in trying to make the code public and useful for folks who
> may want to do similar analysis. Apologies for any errors that come up --
> feel free to contact me if so!

This repo takes that invitation up. The method is Nick's and is unchanged;
what follows is the work of getting a 2017 research release to run against
2025-era OPR products, plus a handful of genuine fixes that are worth
sending back upstream. Most of the list is format drift and buried
constants rather than anything wrong with the science.

### Changes from the published release

Every change is marked in place in the source with a `%%% FIX:` or
`%%% ADD:` comment, so the diff against upstream is readable without a diff
tool.

**Worth folding back upstream.** These are real and would affect anyone
using the public release today:

| Where | Change | Why it matters |
|---|---|---|
| `radon_ndh.m` | Negate `dip_angles` so a layer deepening with increasing x returns a positive dip. | Nick's full NDH_MatlabTools `RollingRadon` applies `slopegrid*-1` at the end; the public release does not carry that line, so its sign convention is the opposite of his own pipeline's. For migration work the sign is the answer. |
| `regrid.m` | Target spacing in mode 1 is a length, `(c_ice/f)/20`, and is converted back to seconds for a travel-time axis. | As released it is `(1/f)/20`, a time, applied to both axes. On a 2.5 km line at 600 MHz the distance axis then asks for ~3e13 samples. |
| `radon_ndh.m` | Normalise to `[0,1]` with `(data-min)/(max-min)`. | The released form divides by the max of the unshifted data, which is negative for dB input. Methods 3 and 4 are not invariant to the resulting sign flip. |
| `radon_ndh.m` | Take the first candidate angle, and return NaN when there is none. | `find(x == max(x))` yields several indices on a tie and none for a featureless window, and callers assign into a scalar. |
| `RollingRadon.m` | `breaks` terminates one past the last column. | Chunks are indexed `breaks(k):breaks(k+1)-1`, so the released value drops the final trace of any line processed in one chunk. |
| `RollingRadon.m` | Test `exist('movie_flag')` rather than `exist('movie')`. | The argument is `movie_flag`, so the animation option could not be switched on. |

**Adaptations for current data and for sub-degree dips.** These are choices
this application needs, not corrections:

| Where | Change | Why it matters here |
|---|---|---|
| `radon_ndh.m` | Angular step `d_theta` is an argument (was 0.1 deg). | Interior layers dip ~0.1 deg, so the whole signal fits inside one search increment at the original step. |
| `RollingRadon.m` | Window may be `[horizontal vertical]` rather than square in samples. | Along-track extent needs enough traces to define a dip; vertical extent wants to stay short enough that system power drift is not the strongest feature in the window. |
| `RollingRadon.m` | Chunk width scales with the window (was a flat 1000 columns). | A window wider than the chunk makes `roll_steps` negative, so the loop does not execute. Long windows are what sub-degree dips require. |
| `RollingRadon.m` | `lp(Data, 1)` for the filename path, and an explicit `data_is_power` flag. | OPR standard and qlook products store detected power, so `10*log10` is the right conversion; the default `20*log10` doubles the dB scale and with it every threshold expressed in dB. |
| `RollingRadon.m` | Microsecond detection tests the record length, not the sample step. | A finely sampled microsecond record has a small step too. |
| `RollingRadon.m` | Treat an empty optional argument as not supplied. | `exist('max_frequency') == 1` is true for a variable that exists but is `[]`, which is the natural placeholder when reaching a later positional argument. |
| `RollingRadon.m` | Gating constants (`snr_thresh`, `vr`, `radon_method`, ...) exposed through a `params` struct; `status_flag` returned. | Tuning previously meant editing the file, and there was no way to see why a window was rejected. |
| `RollingRadon.m` | Optional gate on the Radon peak prominence that `radon_ndh` already returns. | Without it a window containing no reflector still reports whichever angle won, rather than abstaining. |
| `RollingRadon.m` | `vr = 1` by default, so a window that deviates from the one above it keeps its own measurement. | The continuity filter otherwise substitutes the previous value and carries it forward, which draws constant-dip columns through the field. |
| `radon_ndh.m` | Handle a descending y axis by solving in the ascending frame. | `griddedInterpolant` requires ascending sample points, so an elevation axis errors before reaching the branch written for it. |
| `radon_ndh.m` | Cache the Radon correction images and window filters against window size. | They do not depend on the data, but were rebuilt on every call: three extra Radon transforms per window. Caching is what made a 10x finer grid affordable. |
| `regrid.m` | `cice_import.m` included. | It lives in NDH_MatlabTools, so a clean checkout of the slope repo alone cannot regrid radar data. |

`RollingRadon_CReSIS.m` (from NDH_MatlabTools) is superseded rather than
carried over. It is an internal 2016 driver, not part of the published
release: it writes back to its own input file, its chunk loop re-processes
the whole image each iteration, it assumes Windows paths, and it calls
`RadialSpreading`, which is not in either public repository.
`opr/RollingRadon_OPR.m` replaces it.

> **Note on one data file:** `Data_20260109_02_003.mat` in
> `2025_Antarctica_Ground2/CSARP_standard_HH` carries `filename`, `movie`,
> `plotter` and `window` variables, which is the signature of that old
> driver's `save(filename)`. The echogram itself is intact, but the file is
> no longer a pristine OPR product. `opr_load_echogram` warns when it sees
> this so nobody is surprised by it later.

### What the OPR front end does

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
