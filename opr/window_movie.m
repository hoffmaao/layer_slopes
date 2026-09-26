function out_file = window_movie(data_file, out_file, varargin)
% WINDOW_MOVIE  Animate exactly what the Radon transform is being fed.
%
%   window_movie(data_file, out_file, 'name', value, ...)
%
% Steps the solver's own window across the echogram in raster order - left
% to right along a row, then down to the next row - and writes one animated
% GIF frame per window. Each frame shows:
%
%   1. where the window sits in the whole echogram;
%   2. the conditioned window EXACTLY as LS_RADON_DIP receives it, with the
%      fitted slope drawn over it;
%   3. the criterion against candidate slope, with the peak marked.
%
% This is the tool for answering "why did it return that?". A window whose
% criterion curve is flat had nothing to lock onto and should be abstaining;
% one whose fitted line does not lie along the layering is being pulled by
% something else - striping, the power envelope, or a window straddling two
% different dips.
%
% Options
%   z_range     [z0 z1] depth band to scan (m). Default z_pad_surface..z_max.
%   x_range     [x0 x1] along-track span to scan (m). Default: all.
%   max_frames  stop after this many windows. Default 120.
%   stride_x    take every Nth window along a row. Default 1.
%   stride_z    take every Nth row. Default 1.
%   delay       seconds per frame. Default 0.2.
%   semb_thresh semblance gate to mark each window against: a scalar, or
%               the saved result of ROLLINGRADON_OPR (struct or .mat path),
%               whose per-depth calibrated gate is then used. Default []:
%               show the semblance without a verdict.
%   dip_accept  mark windows whose slope exceeds this (true degrees) as
%               rejected, as ROLLINGRADON_OPR does. Default []: no limit.
%   (plus any conditioning/geometry option of ROLLINGRADON_OPR)
%
% See also ROLLINGRADON_OPR, LS_RADON_DIP, LS_ROLLING_RADON

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

o.z_range = [];
o.x_range = [];
o.max_frames = 120;
o.stride_x = 1;
o.stride_z = 1;
o.delay = 0.2;
o.grid_spacing = 0.25;  o.vert_exag = 20;
o.window_x = 1000;      o.window_z = 50;
o.overlap_x = 0.25;     o.overlap_z = 0.25;
o.smooth_x = 60;        o.smooth_len = 1.5;   o.detrend_len = 30;
o.trace_balance = true; o.exclude_z = [];
o.z_pad_surface = 30;   o.z_max = 1500;
o.dip_max = 1;          o.dip_step = 0.005;   o.semb_thresh = [];
o.dip_accept = [];
o.max_excluded = 0.15;

for i = 1:2:numel(varargin)
    if ~isfield(o, varargin{i})
        error('window_movie:unknownOption','Unknown option "%s".', varargin{i});
    end
    o.(varargin{i}) = varargin{i+1};
end

% --- the grid the solver uses ------------------------------------------
D = opr_load_echogram(data_file, struct('verbose',false));
G = opr_flatten_grid(D, struct( ...
    'grid_spacing', o.grid_spacing, 'vert_exag', o.vert_exag, ...
    'smooth_x', o.smooth_x, 'smooth_len', o.smooth_len, ...
    'detrend_len', o.detrend_len, 'trace_balance', o.trace_balance, ...
    'exclude_z', o.exclude_z, 'z_max', o.z_max, 'verbose', false));

[nz, nx] = size(G.img);
wx = max(5, round(o.window_x/G.dx));
wz = max(5, round(o.window_z/G.dz));
sx = max(1, round(wx*o.overlap_x));
sz = max(1, round(wz*o.overlap_z));

c0 = 1:sx:(nx-wx+1);
r0 = 1:sz:(nz-wz+1);

xc_all = G.x(c0 + floor(wx/2));
zc_all = G.z(r0 + floor(wz/2));
if ~isempty(o.x_range)
    keep = xc_all >= o.x_range(1) & xc_all <= o.x_range(2);
    c0 = c0(keep);
end
if isempty(o.z_range)
    o.z_range = [o.z_pad_surface o.z_max];
end
keep = zc_all >= o.z_range(1) & zc_all <= o.z_range(2);
r0 = r0(keep);

c0 = c0(1:o.stride_x:end);
r0 = r0(1:o.stride_z:end);

% The gate to judge each window by, if one was given.
gate = @(z) NaN;
if ischar(o.semb_thresh) || isstring(o.semb_thresh)
    mf = char(o.semb_thresh);
    vars = intersect({'slope_z','semb_thresh'}, who('-file', mf));
    o.semb_thresh = load(mf, vars{:});
end
if isstruct(o.semb_thresh)
    Rg = o.semb_thresh;
    if ~isfield(Rg, 'semb_thresh') || isempty(Rg.semb_thresh)
        % Saved before the gate was stored with the result.
        warning('window_movie:noGate', ...
            'The result holds no semblance gate; showing semblance without a verdict.');
    elseif isscalar(Rg.semb_thresh)
        gate = @(z) Rg.semb_thresh;
    else
        gate = @(z) interp1(Rg.slope_z(:), Rg.semb_thresh(:), z, 'nearest', Inf);
    end
elseif ~isempty(o.semb_thresh)
    gate = @(z) o.semb_thresh;
end

slope_max_app = min(89, atand(o.vert_exag*tand(o.dip_max)));
slope_step_app = atand(o.vert_exag*tand(o.dip_step));
db_lim = prctile(G.raw_db(isfinite(G.raw_db)), [8 99.5]);

fprintf('window_movie: %d x %d windows (%g x %g m), up to %d frames\n', ...
    numel(c0), numel(r0), o.window_x, o.window_z, o.max_frames);

nfr = 0;
for j = 1:numel(r0)                      % row by row ...
    ri = r0(j):(r0(j)+wz-1);
    for i = 1:numel(c0)                  % ... left to right within a row
        if nfr >= o.max_frames, break; end
        ci = c0(i):(c0(i)+wx-1);
        win = G.img(ri, ci);

        [sl_app, ~, crit, slope_axis, sb] = ls_radon_dip(win, G.dz, G.dz, ...
            slope_max_app, slope_step_app, true, 1 - o.max_excluded);
        sl = atand(tand(sl_app)/o.vert_exag);
        thr = gate(mean(G.z(ri)));

        f = ls_figure([50 50 1450 880]);
        tl = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');

        % 1. locator
        ax = nexttile(tl);
        imagesc(ax, G.x/1000, G.z, G.raw_db); colormap(ax,gray); clim(ax,db_lim);
        set(ax,'YDir','reverse'); hold(ax,'on');
        plot(ax, [G.x(ci(1)) G.x(ci(1)) G.x(ci(end)) G.x(ci(end)) G.x(ci(1))]/1000, ...
                 [G.z(ri(1)) G.z(ri(end)) G.z(ri(end)) G.z(ri(1)) G.z(ri(1))], ...
                 'c-', 'LineWidth', 2);
        xlabel(ax,'distance (km)'); ylabel(ax,'depth (m)');

        % 2. the window as the estimator sees it
        ax = nexttile(tl);
        xm = (0:numel(ci)-1)*G.dx;
        imagesc(ax, xm/1000, G.z(ri), win);
        colormap(ax,gray); set(ax,'YDir','reverse'); hold(ax,'on');
        cw = prctile(win(isfinite(win)),[3 97]);
        if numel(cw)==2 && isfinite(cw(1)) && cw(2) > cw(1), clim(ax, cw); end
        if isfinite(sl)
            % A positive slope rises with x, so depth decreases with x.
            plot(ax, xm/1000, mean(G.z(ri)) - (xm-mean(xm))*tand(sl), ...
                'c-', 'LineWidth', 2);
        end
        xlabel(ax,'distance in window (km)'); ylabel(ax,'depth (m)');
        % The same verdict, in the same order, as LS_ROLLING_RADON.
        if isnan(thr)
            title(ax, sprintf('slope %+.4f deg   semblance %.2f', sl, sb));
        elseif ~isfinite(sl) || ~isfinite(sb) || sb < thr
            title(ax, sprintf('slope %+.4f deg   semblance %.2f   rejected (< %.2f)', ...
                sl, sb, thr));
        elseif sl_app <= slope_axis(1) + slope_step_app/2 || ...
                sl_app >= slope_axis(end) - slope_step_app/2
            title(ax, sprintf(['slope %+.4f deg   semblance %.2f   ' ...
                'rejected (best slope at the search edge)'], sl, sb));
        elseif ~isempty(o.dip_accept) && abs(sl) > o.dip_accept
            title(ax, sprintf('slope %+.4f deg   semblance %.2f   rejected (|slope| > %g)', ...
                sl, sb, o.dip_accept));
        else
            title(ax, sprintf('slope %+.4f deg   semblance %.2f   ACCEPTED', sl, sb));
        end

        % 3. the criterion
        ax = nexttile(tl);
        st = atand(tand(slope_axis)/o.vert_exag);    % back to true degrees
        plot(ax, st, crit, 'k-', 'LineWidth', 1.2); hold(ax,'on');
        if isfinite(sl)
            plot(ax, sl, max(crit), 'ro', 'MarkerFaceColor','r');
        end
        grid(ax,'on');
        xlabel(ax,'slope (deg)'); ylabel(ax,'projection variance');
        xlim(ax, [-o.dip_max o.dip_max]);
        title(ax, sprintf('x = %.2f km, z = %.0f m', ...
            mean(G.x(ci))/1000, mean(G.z(ri))));

        drawnow;
        fr = getframe(f);
        [im, cm] = rgb2ind(fr.cdata, 256);
        if nfr == 0
            imwrite(im, cm, out_file, 'gif', 'LoopCount', inf, 'DelayTime', o.delay);
        else
            imwrite(im, cm, out_file, 'gif', 'WriteMode','append', 'DelayTime', o.delay);
        end
        close(f);
        nfr = nfr + 1;
    end
    if nfr >= o.max_frames, break; end
end

fprintf('  wrote %s (%d frames)\n', out_file, nfr);
end
