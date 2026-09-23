function T = validate_horizon(frame, seed_x, seed_z, varargin)
% VALIDATE_HORIZON  Check the slope field against a reflector tracked by eye.
%
%   T = VALIDATE_HORIZON(frame, seed_x, seed_z)
%   T = VALIDATE_HORIZON(..., 'window_x', [500 1000 2000], 'out_dir', dir)
%
%   frame     OPR frame, e.g. '20250112_01_008' (settings from ls_config)
%   seed_x    along-track distance of a point on a bright reflector (m)
%   seed_z    its depth below the surface (m)
%
% Independent ground truth, with no Radon involved: seed on a reflector,
% follow it trace by trace under a continuity constraint, smooth, and
% differentiate. The pick is drawn over the echogram and saved, so it can be
% checked by eye BEFORE it is trusted - a tracker that jumps between
% reflectors produces confident nonsense.
%
% The solver then runs with the settings in LS_CONFIG at each window_x, and
% each window on the horizon is compared with the horizon's chord dip
% across that same window, which is what the window actually measures.
%
% Returns a table with, per window_x: windows on the horizon, the fraction
% of them solved, median solver and tracked dip, their median difference,
% the regression gain of solver on truth, and the correlation.
%
% Examples (on the CReSIS servers)
%   validate_horizon('20250112_01_008', 10000, 130)
%   validate_horizon('20250108_02_005', 10000, 135)

here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(fullfile(repo,'src'), fullfile(repo,'opr'), fullfile(repo,'examples'));

o.window_x = [500 1000 2000];
o.out_dir = '';
o.track_step = 6;          % m, the most the pick may move between traces
for i = 1:2:numel(varargin)
    if ~isfield(o, varargin{i})
        error('validate_horizon:unknownOption','Unknown option "%s".', varargin{i});
    end
    o.(varargin{i}) = varargin{i+1};
end

cfg = ls_config('frame', frame);
if isempty(o.out_dir), o.out_dir = fullfile(cfg.out_dir, 'diag'); end
if exist(o.out_dir,'dir') ~= 7, mkdir(o.out_dir); end

% --- track the reflector ------------------------------------------------
D = opr_load_echogram(cfg.data_file, struct('verbose', false));
G = opr_flatten_grid(D, struct('grid_spacing',0.25,'vert_exag',20, ...
    'z_max',seed_z+60,'smooth_len',0,'detrend_len',0,'trace_balance',false, ...
    'verbose',false));
dz = G.z(2) - G.z(1);
dx = G.x(2) - G.x(1);

% Condition for TRACKING only: remove the depth power envelope so a
% reflector is a local maximum rather than merely shallower than its
% neighbours, and smooth along track so the pick follows the layer, not
% speckle.
A = G.raw_db;
A = A - movmean(A, round(30/dz), 1, 'omitnan');
A = movmean(A, max(1, round(200/dx)), 2, 'omitnan');
A(~isfinite(A)) = -inf;

[~, c0] = min(abs(G.x - seed_x));
near = find(abs(G.z - seed_z) <= 8);
[~, k] = max(A(near, c0));
ridx = nan(1, numel(G.x));
ridx(c0) = near(k);
step = round(o.track_step/dz);
for c = [c0+1:numel(G.x), c0-1:-1:1]
    prev = ridx(c - sign(c - c0));
    lo = max(1, prev - step); hi = min(numel(G.z), prev + step);
    [~, k] = max(A(lo:hi, c));
    ridx(c) = lo + k - 1;
end
zpk = movmean(G.z(ridx), max(1, round(500/dx)));   % the layer is smooth

f = ls_figure([60 60 1700 700]);
ax = axes(f);
imagesc(ax, G.x/1000, G.z, G.raw_db); colormap(ax, gray);
clim(ax, prctile(G.raw_db(isfinite(G.raw_db)), [8 99.5]));
set(ax,'YDir','reverse'); hold(ax,'on');
plot(ax, G.x/1000, zpk, 'r-', 'LineWidth', 1.4);
plot(ax, seed_x/1000, seed_z, 'co', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(ax,'distance (km)'); ylabel(ax,'depth (m)');
pick_png = fullfile(o.out_dir, sprintf('horizon_%s.png', frame));
exportgraphics(f, pick_png, 'Resolution', 110);
close(f);
fprintf('%s: horizon %.1f m at %.1f km -> %.1f m at %.1f km (check %s)\n', ...
    frame, zpk(1), G.x(1)/1000, zpk(end), G.x(end)/1000, pick_png);

% --- compare the solver at each window length -----------------------------
rows = cell(numel(o.window_x), 8);
for w = 1:numel(o.window_x)
    wx = o.window_x(w);
    opts = cfg.opts;
    opts{find(strcmp(opts, 'window_x')) + 1} = wx;
    R = RollingRadon_OPR(cfg.data_file, opts{:}, 'z_pad_surface', cfg.z_top, ...
        'z_max', cfg.z_bot, 'exclude_z', cfg.exclude_z, 'verbose', false);
    sv = []; tv = []; n_on = 0;
    for c = 1:numel(R.slope_x)
        xa = R.slope_x(c) - wx/2; xb = R.slope_x(c) + wx/2;
        if xa < G.x(1) || xb > G.x(end), continue; end
        za = interp1(G.x, zpk, xa); zb = interp1(G.x, zpk, xb);
        zt = interp1(G.x, zpk, R.slope_x(c));
        [gap, r] = min(abs(R.slope_z - zt));
        if gap > R.param.window_z/4 || R.status(r,c) == 1, continue; end
        n_on = n_on + 1;
        if isfinite(R.slopes(r,c))
            sv(end+1) = R.slopes(r,c);            %#ok<AGROW>
            tv(end+1) = -atand((zb - za)/wx);     %#ok<AGROW> + = rises
        end
    end
    if numel(sv) >= 5
        p = polyfit(tv, sv, 1); rr = corr(sv(:), tv(:));
    else
        p = [NaN NaN]; rr = NaN;
    end
    rows(w,:) = {wx, n_on, numel(sv)/max(1,n_on), median(sv), median(tv), ...
        median(sv - tv), p(1), rr};
end
T = cell2table(rows, 'VariableNames', {'window_x','on_horizon','solved', ...
    'solver_med','truth_med','bias_med','gain','r'});
disp(T);
end
