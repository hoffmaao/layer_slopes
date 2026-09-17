function M = slope_multiscale(data_file, varargin)
% SLOPE_MULTISCALE  Layer slope field at the finest scale the data supports.
%
%   M = SLOPE_MULTISCALE(data_file, 'name', value, ...)
%
% Runs ROLLINGRADON_OPR at several window scales and merges them, keeping
% the finest scale whose answer is corroborated by the next coarser one.
%
% WHY MULTIPLE SCALES. The two window dimensions have independent limits on
% accumulation-radar data, and they pull in opposite directions:
%
%   window_z  sets depth detail, and is cheap. The radar resolves 0.53 m
%             and layers sit ~8 m apart, so a 10-15 m window still holds
%             several layer cycles. Small features in depth are resolvable.
%   window_x  sets the smallest measurable dip, and is expensive. A dip of
%             theta displaces a layer by window_x*tan(theta), which has to
%             exceed a useful fraction of one range cell. At 500 m the
%             floor is 0.061 deg - the size of the signal itself, and the
%             field falls apart (sign consistency 0.61). At 2000 m the
%             floor is 0.015 deg and it holds together (0.82).
%
% A single window therefore either resolves small features and cannot see
% small dips, or sees small dips and smooths the features away. Running
% several and merging on agreement gets both: a fine scale is trusted only
% where a coarse scale independently says the same thing, so small features
% survive where they are real and are rejected where they are noise.
%
% Options
%   scales      n x 2 array of [window_x window_z] in metres, coarse first.
%               Default [4000 30; 2000 20; 1000 10].
%   tol         agreement tolerance (deg) between neighbouring scales.
%               Default 0.05.
%   tol_frac    alternative tolerance as a fraction of the coarse |dip|;
%               the looser of the two is used. Default 0.5.
%   out_file    save the merged result here, '' = do not save
%   (any ROLLINGRADON_OPR option is accepted and passed through)
%
% Returns
%   .x, .z          merged grid axes (m)
%   .slopes         merged dip field (deg), NaN where unsupported
%   .scale          [n x m] index of the scale each cell came from (1 = coarsest)
%   .window_x       [n x m] along-track window length behind each cell (m)
%   .per_scale      struct array of the individual ROLLINGRADON_OPR results
%   .param          settings and per-scale statistics
%
% See also ROLLINGRADON_OPR, PLOT_SLOPE_FIELD

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% --- split our options from the solver's --------------------------------
o.scales = [4000 30; 2000 20; 1000 10];
o.tol = 0.05;
o.tol_frac = 0.5;
o.out_file = '';
o.verbose = true;

passthrough = {};
i = 1;
while i <= numel(varargin)
    name = varargin{i};
    if isfield(o, name) && ~strcmp(name,'verbose')
        o.(name) = varargin{i+1};
    elseif strcmp(name,'verbose')
        o.verbose = varargin{i+1};
        passthrough(end+1:end+2) = {name, varargin{i+1}};
    else
        passthrough(end+1:end+2) = varargin(i:i+1); %#ok<AGROW>
    end
    i = i + 2;
end

if size(o.scales,2) ~= 2 || size(o.scales,1) < 1
    error('slope_multiscale:badScales', ...
        'scales must be an n x 2 array of [window_x window_z] in metres.');
end
% Coarsest first, so index 1 is always the most trustworthy fallback.
[~, ord] = sort(o.scales(:,1), 'descend');
o.scales = o.scales(ord,:);
ns = size(o.scales,1);

% --- run each scale -----------------------------------------------------
R = cell(1,ns);
for k = 1:ns
    wx = o.scales(k,1);
    wz = o.scales(k,2);
    if o.verbose
        fprintf('scale %d/%d: %g x %g m (dip floor %.3f deg)\n', ...
            k, ns, wx, wz, atand(0.53/wx));
    end
    args = [passthrough, {'window_x', wx, 'window_z', wz}];
    R{k} = RollingRadon_OPR(data_file, args{:});
end

% --- common grid: the finest scale's sampling ---------------------------
xq = R{ns}.slope_x;
zq = R{ns}.slope_z;
[XQ, ZQ] = meshgrid(xq, zq);

S = nan(numel(zq), numel(xq), ns);
for k = 1:ns
    S(:,:,k) = local_regrid(R{k}, XQ, ZQ);
end

% --- merge: finest scale corroborated by the next coarser one -----------
merged = S(:,:,1);                       % coarsest as the base
scale_idx = ones(size(merged));
scale_idx(~isfinite(merged)) = NaN;

for k = 2:ns
    fine = S(:,:,k);
    coarse = S(:,:,k-1);
    tol = max(o.tol, o.tol_frac*abs(coarse));
    agree = isfinite(fine) & isfinite(coarse) & abs(fine-coarse) <= tol;
    % A finer answer is adopted only where the coarser scale corroborates
    % it. Where they disagree the coarse value stands, because the finer
    % window is the one that cannot see small dips reliably.
    merged(agree) = fine(agree);
    scale_idx(agree) = k;
    % Where the coarse scale had nothing but the fine scale does, there is
    % nothing to corroborate against, so leave it out rather than trust it.
end

M = struct();
M.x = xq;
M.z = zq;
M.slopes = merged;
M.scale = scale_idx;
M.window_x = nan(size(merged));
for k = 1:ns
    M.window_x(scale_idx == k) = o.scales(k,1);
end
M.per_scale = [R{:}];
M.grid = R{ns}.grid;
M.lat = R{ns}.lat; M.lon = R{ns}.lon;
M.x_geo = R{ns}.x;  M.y_geo = R{ns}.y;
M.bed_z = R{ns}.bed_z;

dh = diff(merged, 1, 2);
M.continuity = median(abs(dh(isfinite(dh))));

M.param = R{ns}.param;
M.param.scales = o.scales;
M.param.tol = o.tol;
M.param.tol_frac = o.tol_frac;
M.param.window_x = o.scales(1,1);   % so plot_slope_field can size its ticks

if o.verbose
    v = merged(isfinite(merged));
    fprintf('merged: %d/%d cells (%.1f%%), median %.3f deg, continuity %.3f deg\n', ...
        numel(v), numel(merged), 100*numel(v)/numel(merged), ...
        median(v), M.continuity);
    if ~isempty(v)
        fprintf('  sign consistency %.2f\n', mean(sign(v)==sign(median(v))));
    end
    for k = 1:ns
        n = nnz(scale_idx == k);
        fprintf('  from scale %d (%g x %g m): %d cells (%.0f%%)\n', ...
            k, o.scales(k,1), o.scales(k,2), n, 100*n/max(1,numel(v)));
    end
end

if ~isempty(o.out_file)
    d = fileparts(o.out_file);
    if ~isempty(d) && exist(d,'dir') ~= 7, mkdir(d); end
    save(o.out_file, '-struct', 'M', '-v7.3');
    if o.verbose, fprintf('  wrote %s\n', o.out_file); end
end
end

% ------------------------------------------------------------------------
function Sq = local_regrid(R, XQ, ZQ)
% Put one scale's cells onto the common grid. Natural-neighbour inside the
% solved footprint, nothing outside it.
ok = isfinite(R.slopes);
if nnz(ok) < 4
    Sq = nan(size(XQ));
    return
end
[SX, SZ] = meshgrid(R.slope_x, R.slope_z);
F = scatteredInterpolant(SX(ok), SZ(ok), double(R.slopes(ok)), 'natural', 'none');
Sq = F(XQ, ZQ);
end
