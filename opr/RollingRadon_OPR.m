function R = RollingRadon_OPR(data_file, varargin)
% ROLLINGRADON_OPR  Englacial layer slopes from an OPR / CReSIS echogram.
%
%   R = ROLLINGRADON_OPR(data_file)
%   R = ROLLINGRADON_OPR(data_file, 'name', value, ...)
%
% Replacement for RollingRadon_CReSIS.m, which could not run: it opened
% with load/save on its own input (writing the function's arguments back
% into the shared data product), left `steps` undefined for any line
% shorter than 4000 traces, indexed opt_angle as a scalar and then
% sub-indexed the result, re-processed the whole image on every chunk
% iteration, hard-coded Windows output paths, and called RadialSpreading,
% which is not published in any of Nick's repositories.
%
% Options (name/value)
%   grid_spacing   isotropic working grid (m), default 2
%   window_x       along-track window length (m), default 120
%   window_z       vertical window height (m), default 30
%   window         legacy square window (m); sets both when given
%   range_resolution  vertical resolution of the system (m), default 0.53
%   dip_max        maximum dip searched (deg), default 20
%   dip_accept     dips beyond this are discarded (deg), default dip_max-5
%   z_pad_surface  ignore this far below the surface (m), default 30
%   z_pad_bed      stop this far above the bed (m), default 25
%   z_max          hard depth cap (m), [] = from the bed pick
%   smooth_len     depth low-pass length (m), 0 = off, default 1.5
%   detrend_len    depth high-pass length (m), 0 = off, default 15
%   trace_balance  equalise traces against each other, default true
%   agc_len        running-RMS normalisation (m), 0 = off, default 0
%   layer_file     explicit CSARP_layer path, '' = derive it
%   out_file       save the result here, '' = do not save
%   plotter        1 = debug figure (needs a display), default 0
%   full_window_in_ice  require the whole window inside the ice column,
%                  not just its centre (default true)
%   solver_params  struct forwarded to RollingRadon (snr_thresh, ...)
%   verbose        default true
%
% Returns a struct R with
%   .slope_x    [1 x m] along-track distance of each window centre (m)
%   .slope_z    [1 x n] depth below surface of each window centre (m)
%   .slopes     [n x m] layer dip (deg); positive = deepening with +x
%   .lat,.lon,.x,.y   geolocation of each slope column
%   .bed_z      bed depth at each slope column (m)
%   .param      everything needed to reproduce the run
%
% The input file is opened read-only and is never written to.
%
% See also OPR_LOAD_ECHOGRAM, OPR_FLATTEN_GRID, ROLLINGRADON

p = local_options(varargin{:});

if p.verbose
    fprintf('RollingRadon_OPR\n');
end

% --- read ---------------------------------------------------------------
D = opr_load_echogram(data_file, struct( ...
    'layer_file', p.layer_file, 'verbose', p.verbose));

if ~any(isfinite(D.bed_twtt)) && isempty(p.z_max)
    error('RollingRadon_OPR:noBed', ...
        ['No bed pick is available for "%s" and z_max was not set. ' ...
         'Either supply layer_file or set z_max so the solver knows ' ...
         'where the ice column ends.'], data_file);
end

% --- grid ---------------------------------------------------------------
G = opr_flatten_grid(D, struct( ...
    'grid_spacing', p.grid_spacing, 'z_pad_bed', p.z_pad_bed, ...
    'z_max', p.z_max, 'detrend_len', p.detrend_len, ...
    'smooth_len', p.smooth_len, 'agc_len', p.agc_len, ...
    'trace_balance', p.trace_balance, 'verbose', p.verbose));

% --- window sizing ------------------------------------------------------
window_samples = round([p.window_x p.window_z]/p.grid_spacing);
window_samples = window_samples + (mod(window_samples,2) == 0);
if any(window_samples < 9)
    error('RollingRadon_OPR:windowTooSmall', ...
        ['A %g x %g m window at %g m grid spacing is only %d x %d samples. ' ...
         'Use a larger window or a finer grid.'], p.window_x, p.window_z, ...
        p.grid_spacing, window_samples(1), window_samples(2));
end
if window_samples(1) > size(G.img,2) || window_samples(2) > size(G.img,1)
    error('RollingRadon_OPR:windowTooLarge', ...
        ['A %g x %g m window is %d x %d samples, larger than the %d x %d ' ...
         'grid. Use a smaller window.'], p.window_x, p.window_z, ...
        window_samples(1), window_samples(2), size(G.img,2), size(G.img,1));
end

% The vertical window has to be much larger than the range resolution or
% there is no layering inside it to measure.
if p.window_z < 10*p.range_resolution
    warning('RollingRadon_OPR:thinWindow', ...
        ['A %g m vertical window is only %.0f range resolution cells ' ...
         '(%.2f m). The dip estimate will be poorly constrained.'], ...
        p.window_z, p.window_z/p.range_resolution, p.range_resolution);
end

% Surface/bed gate, in the same depth-below-surface metres as the y axis.
% Where the bed is unpicked, gate on the bottom of the gridded column
% instead of NaN - a NaN bound makes the ice-column test false everywhere,
% which is exactly how the original code produced an all-NaN slope field
% without ever reporting a problem.
bed_gate = G.bed_z;
nanbed = ~isfinite(bed_gate);
if any(nanbed)
    bed_gate(nanbed) = G.z(end);
    if p.verbose
        fprintf('    %d/%d columns have no bed pick; gating those at %.0f m\n', ...
            sum(nanbed), numel(bed_gate), G.z(end));
    end
end
surf_gate = repmat(p.z_pad_surface, 1, numel(G.x));

% RollingRadon tests only the window CENTRE against these bounds, so a
% window whose centre sits just above the bed still has half its height in
% the dead zone below it - and a Radon transform is perfectly happy to
% return a confident dip for noise. Pull the bounds in by half a window so
% the whole window has to be inside the ice.
if p.full_window_in_ice
    surf_gate = surf_gate + p.window_z/2;
    bed_gate = bed_gate - p.window_z/2;
end
surface_bottom = [surf_gate; bed_gate];

if all(bed_gate <= surf_gate)
    error('RollingRadon_OPR:noRoom', ...
        ['A %g m vertical window leaves no depth range fully inside the ' ...
         'ice (usable column is %.0f m). Use a smaller window, or set ' ...
         'full_window_in_ice to false.'], p.window_z, ...
        max(G.bed_z) - p.z_pad_surface);
end

if p.verbose
    fprintf('    window %g x %g m = %d x %d samples (%.0f traces, %.0f range cells)\n', ...
        p.window_x, p.window_z, window_samples(1), window_samples(2), ...
        p.window_x/median(diff(D.dist)), p.window_z/p.range_resolution);
    fprintf('    dip search +/- %g deg\n', p.dip_max);
    fprintf('  running rolling radon ...\n');
end

% --- solve --------------------------------------------------------------
% The grid is already isotropic, so regrid() inside RollingRadon takes its
% skip-interpolation branch and radon_ndh does no aspect correction. Note
% max_frequency is deliberately NOT passed: supplying it would send regrid
% down the resampling path this driver has already done properly.
t0 = tic;
[~, ~, ~, opt_x, opt_y, opt_angle, status_flag] = RollingRadon( ...
    G.x, G.z, G.img, window_samples, [p.dip_max p.dip_accept], ...
    p.plotter, surface_bottom, 0, [], [], p.solver_params);
runtime = toc(t0);

% --- assemble -----------------------------------------------------------
keep = opt_x ~= 0;
keep(1) = true;                       % the first centre legitimately can be 0
slope_x = opt_x(keep);
slopes = opt_angle(:, keep(1:size(opt_angle,2)));
status = status_flag(:, keep(1:size(status_flag,2)));

R = struct();
R.slope_x = slope_x;
R.slope_z = opt_y;
R.slopes = slopes;
R.status = status;        % 0 solved, 1 outside ice, 2 low SNR, 3 slope rejected
R.bed_z = interp1(G.x, G.bed_z, slope_x, 'linear', NaN);
R.lat = interp1(D.dist, D.lat, slope_x, 'linear', NaN);
R.lon = interp1(D.dist, D.lon, slope_x, 'linear', NaN);
R.x = interp1(D.dist, D.x, slope_x, 'linear', NaN);
R.y = interp1(D.dist, D.y, slope_x, 'linear', NaN);
R.grid = struct('x', G.x, 'z', G.z, 'bed_z', G.bed_z, ...
    'grid_spacing', G.grid_spacing, 'c_ice', G.c_ice);
R.param = p;
R.param.data_file = data_file;
R.param.window_samples = window_samples;
R.param.bed_source = D.bed_source;
R.param.surface_source = D.surface_source;
R.param.runtime_s = runtime;
R.param.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

nfin = sum(isfinite(R.slopes(:)));
R.param.n_solved = nfin;
R.param.n_windows = numel(R.slopes);

if p.verbose
    fprintf('  done in %.1f min: %d/%d windows solved (%.1f%%)\n', ...
        runtime/60, nfin, numel(R.slopes), 100*nfin/max(1,numel(R.slopes)));
    if nfin > 0
        v = R.slopes(isfinite(R.slopes));
        fprintf('    dip: median %.2f deg, IQR %.2f to %.2f, range %.2f to %.2f\n', ...
            median(v), prctile(v,25), prctile(v,75), min(v), max(v));
    end
    rej = R.status(~isfinite(R.slopes));
    n = numel(R.slopes);
    fprintf('    rejected: %d outside ice (%.0f%%), %d low SNR (%.0f%%), %d slope gate (%.0f%%)\n', ...
        sum(rej==1), 100*sum(rej==1)/n, sum(rej==2), 100*sum(rej==2)/n, ...
        sum(rej==3), 100*sum(rej==3)/n);
end

% --- save ---------------------------------------------------------------
% Always to a separate file. Writing results back into the echogram is the
% bug that contaminated the shared product in the first place.
if ~isempty(p.out_file)
    out_dir = fileparts(p.out_file);
    if ~isempty(out_dir) && exist(out_dir,'dir') ~= 7
        mkdir(out_dir);
    end
    save(p.out_file, '-struct', 'R', '-v7.3');
    if p.verbose
        fprintf('  wrote %s\n', p.out_file);
    end
end
end

% ------------------------------------------------------------------------
function p = local_options(varargin)
p.grid_spacing = 2;
p.window = [];          % legacy square window (m); sets both if given
p.window_x = 120;
p.window_z = 30;
p.dip_max = 20;
p.dip_accept = [];
p.z_pad_surface = 30;
p.z_pad_bed = 25;
p.z_max = [];
p.detrend_len = 15;
p.smooth_len = 1.5;
p.agc_len = 0;
p.trace_balance = true;
p.layer_file = '';
p.out_file = '';
p.plotter = 0;
p.full_window_in_ice = true;
p.range_resolution = 0.53;   % m in ice, accum radar (~159 MHz bandwidth)
p.solver_params = struct();
p.verbose = true;

if mod(numel(varargin),2) ~= 0
    error('RollingRadon_OPR:badArgs','Options must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    if ~isfield(p, name)
        error('RollingRadon_OPR:unknownOption', ...
            'Unknown option "%s". Valid: %s.', name, strjoin(fieldnames(p)',', '));
    end
    p.(name) = varargin{i+1};
end

% Solver defaults that differ from RollingRadon's built-ins, for reasons
% that matter on real data:
%
%   vr = 1      RollingRadon's continuity filter REPLACES a window that
%               deviates from the one above it with the previous value,
%               and then carries that value forward. On a real line this
%               paints long constant-dip columns through the slope field -
%               fabricated numbers that look like signal. vr = 1 makes the
%               filter accept each measurement instead of substituting.
%               Raise it to re-enable Holschuh's original smoothing.
%   snr_thresh  2 dB passes almost anything once the image is depth
%               detrended; 4 dB actually discriminates.
defaults = struct('vr', 1, 'snr_thresh', 4);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(p.solver_params, fn{i})
        p.solver_params.(fn{i}) = defaults.(fn{i});
    end
end

if ~isempty(p.window)
    p.window_x = p.window;
    p.window_z = p.window;
end
p.window = [p.window_x p.window_z];

if isempty(p.dip_accept)
    p.dip_accept = p.dip_max - 5;
end
if p.dip_accept > p.dip_max
    error('RollingRadon_OPR:badDip', ...
        'dip_accept (%g) cannot exceed dip_max (%g).', p.dip_accept, p.dip_max);
end
if p.dip_max <= 0 || p.dip_max >= 90
    error('RollingRadon_OPR:badDip','dip_max must be in (0,90).');
end
if p.grid_spacing <= 0 || any(p.window <= 0)
    error('RollingRadon_OPR:badSize','grid_spacing and window must be positive.');
end
end
