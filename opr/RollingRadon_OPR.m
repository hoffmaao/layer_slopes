function R = RollingRadon_OPR(data_file, varargin)
% ROLLINGRADON_OPR  Englacial layer slopes from an OPR / CReSIS echogram.
%
%   R = ROLLINGRADON_OPR(data_file)
%   R = ROLLINGRADON_OPR(data_file, 'name', value, ...)
%
% Reads an echogram, conditions it, rolls a window over it, and fits the
% dominant layer slope in each window with a Radon transform, after
% Holschuh et al. (2017). The input file is opened read-only.
%
% SIGN CONVENTION: positive slope = the layer RISES (gets shallower) with
% increasing along-track distance, i.e. d(elevation)/dx.
%
% Options (name/value)
%   grid_spacing   vertical grid spacing (m), default 0.25
%   vert_exag      along-track sampling = grid_spacing*vert_exag, presented
%                  to the Radon as isotropic. Multiplies the apparent slope
%                  by exactly vert_exag, undone on output. Default 20.
%   window_x       along-track window length (m), default 1000
%   window_z       vertical window height (m), default 20
%   overlap_x      step between centres as a fraction of the window,
%   overlap_z      default 0.25 (a quarter window)
%   dip_max        search +/- this, true degrees, default 1
%   dip_accept     discard results beyond this, default dip_max
%   dip_step       search step, true degrees, default 0.005
%   q_thresh       minimum Radon peak/median criterion ratio, default 1.5.
%                  Incoherent noise scores about 1.
%   z_pad_surface  ignore this far below the surface (m), default 30
%   z_pad_bed      stop this far above the bed (m), default 25
%   z_max          hard depth cap (m), [] = from the bed pick
%   bed_default    assumed ice thickness (m) when no bed pick, default 1500
%   exclude_z      n x 2 depth bands (m) to blank, e.g. a merged pulse return
%   smooth_x       along-track low-pass (m), default 60
%   smooth_len     depth low-pass (m), default 1.5
%   detrend_len    depth high-pass (m), default 30. Removes the
%                  power-vs-depth gradient, which is horizontal and
%                  otherwise dominates the window.
%   trace_balance  equalise traces against each other, default true
%   full_window_in_ice  require the whole window inside the gates, default true
%   range_resolution    vertical resolution of the system (m), default 0.53
%   layer_file     explicit CSARP_layer path, '' = derive it
%   out_file       save the result here, '' = do not save
%   verbose        default true
%
% Returns a struct R with
%   .slope_x, .slope_z   window centres (m)
%   .slopes              [n x m] slope in degrees (see convention above)
%   .q                   [n x m] Radon criterion peak/median
%   .status              0 solved, 1 outside ice, 2 low quality, 3 slope gate
%   .lat,.lon,.x,.y      geolocation of each slope column
%   .bed_z               bed depth at each column (m)
%   .param               everything needed to reproduce the run
%
% See also LS_ROLLING_RADON, LS_RADON_DIP, OPR_FLATTEN_GRID, PLOT_SLOPE_FIELD

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

p = local_options(varargin{:});

if p.verbose
    fprintf('RollingRadon_OPR\n');
end

% --- read ---------------------------------------------------------------
D = opr_load_echogram(data_file, struct( ...
    'layer_file', p.layer_file, 'verbose', p.verbose));

if ~any(isfinite(D.bed_twtt)) && isempty(p.z_max) && p.verbose
    fprintf('    no bed pick; assuming %.0f m of ice (bed_default)\n', p.bed_default);
end

% --- grid ---------------------------------------------------------------
G = opr_flatten_grid(D, struct( ...
    'grid_spacing', p.grid_spacing, 'vert_exag', p.vert_exag, ...
    'z_pad_bed', p.z_pad_bed, 'z_max', p.z_max, ...
    'bed_default', p.bed_default, 'exclude_z', p.exclude_z, ...
    'detrend_len', p.detrend_len, 'smooth_len', p.smooth_len, ...
    'smooth_x', p.smooth_x, 'trace_balance', p.trace_balance, ...
    'verbose', p.verbose));

% --- gates --------------------------------------------------------------
% Where the bed is unpicked, gate on the bottom of the gridded column
% rather than NaN: a NaN bound makes the ice-column test false everywhere
% and yields an all-NaN field with no indication of why.
bed_gate = G.bed_z;
nanbed = ~isfinite(bed_gate);
if any(nanbed)
    bed_gate(nanbed) = G.z(end);
end
surf_gate = repmat(p.z_pad_surface, 1, numel(G.x));

if p.window_z < 10*p.range_resolution
    warning('RollingRadon_OPR:thinWindow', ...
        ['A %g m vertical window is only %.0f range-resolution cells. ' ...
         'The slope will be poorly constrained.'], ...
        p.window_z, p.window_z/p.range_resolution);
end

% --- solve --------------------------------------------------------------
S = ls_rolling_radon(G, struct( ...
    'window_x', p.window_x, 'window_z', p.window_z, ...
    'overlap_x', p.overlap_x, 'overlap_z', p.overlap_z, ...
    'slope_max', p.dip_max, 'slope_step', p.dip_step, ...
    'slope_accept', p.dip_accept, 'q_thresh', p.q_thresh, ...
    'surface_z', surf_gate, 'bed_z', bed_gate, ...
    'whole_window_in_ice', p.full_window_in_ice, 'verbose', p.verbose));

% --- assemble -----------------------------------------------------------
R = struct();
R.slope_x = S.slope_x;
R.slope_z = S.slope_z;
R.slopes = S.slopes;
R.q = S.q;
R.status = S.status;
R.bed_z = interp1(G.x, G.bed_z, S.slope_x, 'linear', NaN);
R.lat = interp1(D.dist, D.lat, S.slope_x, 'linear', NaN);
R.lon = interp1(D.dist, D.lon, S.slope_x, 'linear', NaN);
R.x = interp1(D.dist, D.x, S.slope_x, 'linear', NaN);
R.y = interp1(D.dist, D.y, S.slope_x, 'linear', NaN);
R.grid = struct('x', G.x, 'z', G.z, 'bed_z', G.bed_z, ...
    'grid_spacing', G.grid_spacing, 'vert_exag', G.vert_exag, 'c_ice', G.c_ice);

dh = diff(R.slopes, 1, 2);
R.continuity = median(abs(dh(isfinite(dh))));

R.param = p;
R.param.data_file = data_file;
R.param.window_samples = S.window_samples;
R.param.bed_source = D.bed_source;
R.param.surface_source = D.surface_source;
R.param.runtime_s = S.runtime_s;
R.param.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

nfin = sum(isfinite(R.slopes(:)));
R.param.n_solved = nfin;
R.param.n_windows = numel(R.slopes);

if p.verbose
    fprintf('  done in %.1f min: %d/%d windows solved (%.1f%%)\n', ...
        S.runtime_s/60, nfin, numel(R.slopes), ...
        100*nfin/max(1,numel(R.slopes)));
    if nfin > 0
        v = R.slopes(isfinite(R.slopes));
        fprintf('    slope: median %+.3f deg, IQR %+.3f to %+.3f, range %+.3f to %+.3f\n', ...
            median(v), prctile(v,25), prctile(v,75), min(v), max(v));
        fprintf('    sign consistency %.2f, continuity %.3f deg\n', ...
            mean(sign(v)==sign(median(v))), R.continuity);
    end
    st = R.status(~isfinite(R.slopes));
    n = numel(R.slopes);
    fprintf('    rejected: %d outside ice (%.0f%%), %d low quality (%.0f%%), %d slope gate (%.0f%%)\n', ...
        sum(st==1), 100*sum(st==1)/n, sum(st==2), 100*sum(st==2)/n, ...
        sum(st==3), 100*sum(st==3)/n);
end

% --- save ---------------------------------------------------------------
% Always to a separate file; the echogram is never written to.
if ~isempty(p.out_file)
    d = fileparts(p.out_file);
    if ~isempty(d) && exist(d,'dir') ~= 7, mkdir(d); end
    save(p.out_file, '-struct', 'R', '-v7.3');
    if p.verbose, fprintf('  wrote %s\n', p.out_file); end
end
end

% ------------------------------------------------------------------------
function p = local_options(varargin)
p.grid_spacing = 0.25;
p.vert_exag = 20;
p.window_x = 1000;
p.window_z = 20;
p.overlap_x = 0.25;
p.overlap_z = 0.25;
p.dip_max = 1;
p.dip_accept = [];
p.dip_step = 0.005;
p.q_thresh = 1.5;
p.z_pad_surface = 30;
p.z_pad_bed = 25;
p.z_max = [];
p.bed_default = 1500;
p.exclude_z = [];
p.smooth_x = 60;
p.smooth_len = 1.5;
p.detrend_len = 30;
p.trace_balance = true;
p.full_window_in_ice = true;
p.range_resolution = 0.53;
p.layer_file = '';
p.out_file = '';
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

if isempty(p.dip_accept), p.dip_accept = p.dip_max; end
if p.dip_accept > p.dip_max
    error('RollingRadon_OPR:badDip', ...
        'dip_accept (%g) cannot exceed dip_max (%g).', p.dip_accept, p.dip_max);
end
if p.dip_max <= 0 || p.dip_max >= 90
    error('RollingRadon_OPR:badDip','dip_max must be in (0,90).');
end
if p.grid_spacing <= 0 || p.window_x <= 0 || p.window_z <= 0
    error('RollingRadon_OPR:badSize', ...
        'grid_spacing and window sizes must be positive.');
end
end
