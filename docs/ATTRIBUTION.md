# Attribution and provenance

## Origin of the method

The rolling Radon transform approach to englacial layer slopes is Nick
Holschuh's:

> N. Holschuh, B. R. Parizek, R. B. Alley, S. Anandakrishnan (2017),
> *Decoding ice sheet behavior using englacial layer slopes*,
> Geophysical Research Letters 44, 5561-5570.
> <https://doi.org/10.1002/2017GL073417>

Please cite that paper for the method.

## Origin of the code

Everything in `src/` other than `cice_import.m` was taken from

> <https://github.com/nholschuh/SlopeExtraction_Radar>
> at commit `df8a2d1b548570c1279de3a68c7a552a4de00310`
> ("Added depth_shift and elevation_shift to the repo")

`cice_import.m` is a two-line reimplementation of a script `regrid.m` calls
but the public repository does not ship (it lives in
<https://github.com/nholschuh/NDH_MatlabTools>). It is written here rather
than copied.

The upstream README is preserved verbatim at `docs/UPSTREAM_README`.

### Licensing

`nholschuh/SlopeExtraction_Radar` carries no licence file, so no explicit
terms of redistribution are attached to it. The code is republished here with
attribution, in the spirit of the upstream README's stated intent to make it
"public and useful for folks who may want to do similar analysis", and with
every modification marked. If you plan to build on this, contact Nick
(Nick.Holschuh@gmail.com) - both to let him know and because several of the
fixes below are worth folding back upstream.

## Files and their status

| File | Status |
|---|---|
| `src/RollingRadon.m` | Holschuh, modified - see `%%% FIX:` / `%%% ADD:` comments |
| `src/radon_ndh.m` | Holschuh, modified - sign convention, normalisation, tie/empty handling, descending axis, filter caching |
| `src/regrid.m` | Holschuh, modified - mode-1 target spacing is a length, not a time |
| `src/cice_import.m` | Written here (upstream ships it in a different repo) |
| `src/b2r2.m`, `combvec.m`, `depth_shift.m`, `distance_vector.m`, `elevation_shift.m`, `exclude.m`, `find_nearest.m`, `interpNaN.m`, `lp.m`, `matrix_to_vector.m`, `plot_indicator_lines.m`, `pointdistance.m`, `polarstereo_fwd.m`, `rad2deg.m`, `value2value.m` | Holschuh, unmodified |
| `opr/*` | Written here |
| `examples/*`, `tests/*` | Written here |

Every edit to an upstream file is marked in place with a `%%% FIX:` or
`%%% ADD:` comment explaining what was wrong and why the change is correct, so
the diff against upstream is readable without a diff tool.

## Not carried over

`RollingRadon_CReSIS.m` (from NDH_MatlabTools) is deliberately absent. It is
superseded by `opr/RollingRadon_OPR.m`. The reasons it was not repaired in
place are listed in the README and in the header of its replacement; the
short version is that it writes to its own input file, its chunking loop does
not actually chunk, and it depends on `RadialSpreading`, which is not
published in any of Nick's public repositories.
