function n_fail = test_search_edge()
% TEST_SEARCH_EDGE  A best slope at the edge of the search is not reported.
%
% When the criterion peaks at either end of the slope search, the real peak
% lies beyond it, and returning the end value would report a clipped slope
% as a measurement. Layers rising 3 deg, searched over only +/-1 deg, must
% abstain with status 4; searched over +/-6 deg they must come back at 3.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));
n_fail = 0;

rng(3);
nx = 400; nz = 200; d = 1;
[X, Z] = meshgrid((0:nx-1)*d, (0:nz-1)*d);
img = sin(2*pi*(Z + tand(3)*X)/9) + 0.3*randn(size(X));   % rises 3 deg
G = struct('x', X(1,:), 'z', Z(:,1).', 'img', img, 'dx', d, 'dz', d, ...
    'vert_exag', 1);
base = struct('window_x', 100, 'window_z', 60, 'slope_step', 0.1, ...
    'semb_thresh', 0.2, 'verbose', false);

p = base; p.slope_max = 1;
S = ls_rolling_radon(G, p);
if any(isfinite(S.slopes(:))) || ~all(S.status(:) == 4)
    fprintf('  FAIL search +/-1 deg: %d windows reported, %d of %d at status 4\n', ...
        nnz(isfinite(S.slopes)), nnz(S.status == 4), numel(S.status));
    n_fail = n_fail + 1;
else
    fprintf('  ok   3 deg layers searched over +/-1 deg abstain (status 4)\n');
end

p = base; p.slope_max = 6;
S = ls_rolling_radon(G, p);
v = S.slopes(isfinite(S.slopes));
if isempty(v) || abs(median(v) - 3) > 0.05
    fprintf('  FAIL search +/-6 deg: median %.3f, want 3\n', median(v));
    n_fail = n_fail + 1;
else
    fprintf('  ok   searched over +/-6 deg they come back at %.3f deg\n', median(v));
end

if n_fail == 0
    fprintf('test_search_edge: PASS\n');
else
    fprintf('test_search_edge: %d FAILURE(S)\n', n_fail);
end
end
