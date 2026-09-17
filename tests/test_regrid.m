function n_fail = test_regrid()
% TEST_REGRID  regrid mode 1 sizes its target grid in the right units.
%
% The released regrid.m computed its target spacing as (1/f)/20, a TIME,
% and then applied it to the distance axis as well. For a 2.5 km line at
% 600 MHz that asks for ~3e13 samples, so the call died in the allocator.
% The step is a LENGTH, (c_ice/f)/20, converted back to seconds only for a
% travel-time axis.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

n_fail = 0;
cice_import

% Geometry matching the Ridge A example frame.
nt = 4000;
dt = 1.6667e-9;
twtt = (0:nt-1)'*dt;
dist = (0:423)*5.9;
data = repmat(sin(2*pi*(1:nt)'/50), 1, numel(dist));
f = 600e6;

t0 = tic;
[nx, ny, out] = regrid(dist, twtt, data, 1, f);
elapsed = toc(t0);

% 1. it returns at all, and quickly
if elapsed > 60
    fprintf('  FAIL regrid took %.1f s\n', elapsed); n_fail = n_fail+1;
else
    fprintf('  ok   regrid returned in %.2f s\n', elapsed);
end

% 2. the grids are physically sized. The broken version asked for ~3e13
%    samples along x alone; anything in that neighbourhood is the bug.
if numel(nx)*numel(ny) > 1e8
    fprintf('  FAIL output grid is %d x %d (%.3g elements)\n', ...
        numel(ny), numel(nx), numel(nx)*numel(ny));
    n_fail = n_fail+1;
else
    fprintf('  ok   output grid %d x %d (%.3g elements)\n', ...
        numel(ny), numel(nx), numel(nx)*numel(ny));
end

% 2b. the x axis is stepped by a LENGTH, not by a travel time
dxo = median(diff(nx));
if ~(dxo > 1e-3 && dxo < 100)
    fprintf('  FAIL x step is %.3g - not a plausible distance in metres\n', dxo);
    n_fail = n_fail+1;
else
    fprintf('  ok   x step is %.3f m\n', dxo);
end

% 3. the x axis still spans the input track, in metres
if abs(nx(end) - dist(end)) > 2*median(diff(nx))
    fprintf('  FAIL x axis ends at %g, track ends at %g\n', nx(end), dist(end));
    n_fail = n_fail+1;
else
    fprintf('  ok   x axis spans the track (%.0f m)\n', nx(end));
end

% 4. the y axis is still two-way travel time, not metres
if ny(end) > 1e-3
    fprintf('  FAIL y axis looks like metres (ends at %g)\n', ny(end));
    n_fail = n_fail+1;
else
    fprintf('  ok   y axis is still twtt (ends at %.3g s)\n', ny(end));
end

% 5. the output matches the axes
if ~isequal(size(out), [numel(ny) numel(nx)])
    fprintf('  FAIL output is %s, axes imply %s\n', ...
        mat2str(size(out)), mat2str([numel(ny) numel(nx)]));
    n_fail = n_fail+1;
else
    fprintf('  ok   output shape agrees with the axes\n');
end

% 6. an already-isotropic grid is passed through untouched, which is what
%    lets the OPR driver pre-grid once instead of per window
gs = 2;
xi = (0:99)*gs; zi = (0:79)*gs;
d2 = rand(80,100);
[nx2, ny2, out2] = regrid(xi, zi, d2, 0, 0);
if isequal(out2, d2) && isequal(nx2(:).', xi) && isequal(ny2(:).', zi)
    fprintf('  ok   isotropic input passes through unchanged\n');
else
    fprintf('  FAIL isotropic input was resampled\n');
    n_fail = n_fail+1;
end

if n_fail == 0
    fprintf('test_regrid: PASS\n');
else
    fprintf('test_regrid: %d FAILURE(S)\n', n_fail);
end
end
