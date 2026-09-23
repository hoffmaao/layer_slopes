function fig_file = plot_slope_field(R, out_png, varargin)
% PLOT_SLOPE_FIELD  Standard figure: layer slope field over the power image.
%
%   plot_slope_field(R, out_png)
%   plot_slope_field(slope_file, out_png)
%   plot_slope_field(..., 'name', value)
%
% R is the struct returned by ROLLINGRADON_OPR, or the path to a saved one.
% The background is the echogram power in dB (10*log10 of the OPR power
% product), surface-flattened onto the same grid the solver worked on, so
% every slope cell sits exactly over the layering it was measured from.
%
% Options
%   clim_dip     dip colour limits (deg), default symmetric round the p98
%   alpha        slope overlay opacity, default 0.6
%   interp       render the field as a continuous raster rather than the
%                raw window cells, default true
%   interp_smooth  cells of smoothing applied to that raster, default 3
%   despeckle    draw each solved cell as the median of itself and its
%                solved neighbours, default true. A single window that
%                disagrees in sign with all its neighbours otherwise comes
%                out of the interpolation as a bullseye ringed in white.
%                Display only: the saved slopes are never changed.
%   segments     draw a dip tick in each solved cell; default is to draw
%                them only when they would be visibly tilted
%   seg_len      half-length of those ticks (m), default window_x/6
%   dpi          output resolution, default 150
%   title_str    optional figure title; empty (the default) draws none
%
% Returns the path written.
%
% See also ROLLINGRADON_OPR, OPR_FLATTEN_GRID

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

o.clim_dip = [];
o.alpha = 0.6;
o.interp = true;        % render the field as a continuous raster
o.interp_smooth = [];   % smoothing of that raster; [] scales to the
                        % window spacing, which is what the banding is
o.despeckle = true;
o.segments = [];   % default: only when the ticks would be visible
o.seg_len = [];     % default: scale to the window
o.dpi = 150;
o.title_str = '';
for i = 1:2:numel(varargin)
    if ~isfield(o, varargin{i})
        error('plot_slope_field:unknownOption','Unknown option "%s".', varargin{i});
    end
    o.(varargin{i}) = varargin{i+1};
end

if ischar(R) || isstring(R)
    R = load(char(R));
end

% SLOPE_MULTISCALE returns the merged field on .x/.z; ROLLINGRADON_OPR uses
% .slope_x/.slope_z. Accept either.
if ~isfield(R,'slope_x') && isfield(R,'x') && isfield(R,'slopes')
    R.slope_x = R.x;
    R.slope_z = R.z;
end

if isempty(o.seg_len)
    if isfield(R.param,'window_x')
        o.seg_len = R.param.window_x/6;
    else
        o.seg_len = 20;
    end
end

% Rebuild the exact grid the solver saw. Copy across only the fields that
% are actually present, so a result saved by an older or newer version
% still plots instead of erroring on one renamed option.
D = opr_load_echogram(R.param.data_file, struct('verbose', false));
gopt = struct('verbose', false);
for f = {'grid_spacing','z_pad_bed','z_max','bed_default','detrend_len', ...
         'smooth_len','trace_balance','vert_exag','smooth_x','exclude_z'}
    if isfield(R.param, f{1})
        gopt.(f{1}) = R.param.(f{1});
    end
end
G = opr_flatten_grid(D, gopt);

% Lines run to tens of km; metres force an exponent onto the axis.
xkm = G.x/1000;
sxkm = R.slope_x/1000;

dB = G.raw_db;                       % 10*log10(power), surface-flattened
finite_db = dB(isfinite(dB));
db_lim = prctile(finite_db, [8 99.5]);

v = R.slopes(isfinite(R.slopes));
if isempty(o.clim_dip)
    if isempty(v)
        o.clim_dip = [-5 5];
    else
        % No fixed floor: interior layers dip a tenth of a degree, and a
        % hard minimum of +/-1 deg flattens the whole field to one pale
        % colour. Scale to the data, robustly enough to ignore outliers.
        m = prctile(abs(v), 95);
        if ~isfinite(m) || m <= 0
            m = max(abs(v));
        end
        if ~isfinite(m) || m <= 0
            m = 1;
        end
        % Round up to 1, 2, 2.5 or 5 times a power of ten, so the scale
        % ends on a labelled value.
        p10 = 10^floor(log10(m));
        c = [1 2 2.5 5 10]*p10;
        m = c(find(c >= m*(1 - 1e-9), 1));
        o.clim_dip = [-m m];
    end
end

fig = ls_figure([80 80 1750 1080]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

% ---- panel 1: the power image ------------------------------------------
ax1 = nexttile(tl);
imagesc(ax1, xkm, G.z, dB);
colormap(ax1, gray); clim(ax1, db_lim);
set(ax1,'YDir','reverse','Layer','top','TickDir','out','Box','on');
hold(ax1,'on');
if any(isfinite(G.bed_z))
    plot(ax1, xkm, G.bed_z, 'r:', 'LineWidth', 1.6);
end
ylabel(ax1,'depth (m)');
set(ax1,'XTickLabel',[]);
cb1 = colorbar(ax1); cb1.Label.String = 'power (dB)';

% ---- panel 2: slope field over the power image -------------------------
ax2 = nexttile(tl);
imagesc(ax2, xkm, G.z, dB);
colormap(ax2, gray); clim(ax2, db_lim);
set(ax2,'YDir','reverse','Layer','top','TickDir','out','Box','on');
hold(ax2,'on');

% Overlay the dips in their own axes so the two colormaps coexist.
ax3 = axes('Position', ax2.Position, 'Color','none');

% Rolling windows overlap, so plotting the raw cell grid draws a staircase
% of rectangles whose edges are an artefact of the window spacing rather
% than anything in the ice. Interpolating onto a fine grid gives the
% continuous slope raster of Holschuh et al. (2017, fig. 3), where the
% gradient itself is the thing being read. Nick's full RollingRadon has
% interp_method options for this; the public release ships with it off.
if o.interp && nnz(isfinite(R.slopes)) >= 1 && numel(R.slope_x) >= 2 ...
        && numel(R.slope_z) >= 2
    nxq = min(1600, 8*numel(R.slope_x));
    nzq = min(600, 8*numel(R.slope_z));
    xq = linspace(min(sxkm), max(sxkm), nxq);
    zq = linspace(min(R.slope_z), max(R.slope_z), nzq);
    % Adjacent windows overlap heavily, so their estimates differ slightly
    % and leave banding at exactly the window spacing. Smooth over one
    % cell spacing in each direction to remove it without touching the
    % gradient itself.
    if isempty(o.interp_smooth)
        ns_x = max(3, round(nxq/max(1,numel(R.slope_x))));
        ns_z = max(3, round(nzq/max(1,numel(R.slope_z))));
    else
        ns_x = round(o.interp_smooth); ns_z = ns_x;
    end
    % Interpolate ONLY between solved neighbours. Value and support are
    % interpolated separately (normalised convolution) and the raster is
    % drawn only where there is support, fading out as support falls, so a
    % gap between solved windows stays a gap. A scattered interpolant over
    % the convex hull would instead stretch a handful of isolated cells
    % into streaks kilometres long, which reads as structure that is not in
    % the data.
    ok = isfinite(R.slopes);
    V = double(R.slopes);
    if o.despeckle
        V = local_despeckle(V);
    end
    V(~ok) = 0;
    [QX, QZ] = meshgrid(xq, zq);
    Vq = interp2(sxkm, R.slope_z, V, QX, QZ, 'linear', 0);
    Wq = interp2(sxkm, R.slope_z, double(ok), QX, QZ, 'linear', 0);
    if ns_x > 1 || ns_z > 1
        Vq = movmean(movmean(Vq, ns_z, 1), ns_x, 2);
        Wq = movmean(movmean(Wq, ns_z, 1), ns_x, 2);
    end
    S = Vq./Wq;
    fade = min(1, max(0, (Wq - 0.3)/0.4));     % 0 below 30% support, 1 above 70%
    S(fade <= 0) = NaN;
    h = imagesc(ax3, xq, zq, S);
    set(h, 'AlphaData', fade*o.alpha);
else
    h = imagesc(ax3, sxkm, R.slope_z, R.slopes);
    set(h, 'AlphaData', isfinite(R.slopes)*o.alpha);
end
colormap(ax3, ls_slope_colormap());
clim(ax3, o.clim_dip);
set(ax3,'YDir','reverse','Color','none','XTick',[],'YTick',[],'Box','off');
hold(ax3,'on');

% A dip tick is only informative if it is actually visible. At interior
% dips of ~0.1 deg over a window-scaled tick the displacement is a few tens
% of centimetres against a panel hundreds of metres deep, so the ticks
% collapse to flat lines and only add clutter. Drop them automatically.
if isempty(o.segments)
    if isempty(v)
        o.segments = false;
    else
        rise = o.seg_len*abs(tand(prctile(abs(v),90)));
        o.segments = rise > 0.01*(max(G.z)-min(G.z));
    end
end

if o.segments && ~isempty(v)
    for i = 1:numel(R.slope_x)
        for j = 1:numel(R.slope_z)
            d = R.slopes(j,i);
            if ~isfinite(d), continue; end
            plot(ax3, sxkm(i) + [-1 1]*o.seg_len/1000, ...
                      R.slope_z(j) + [-1 1]*o.seg_len*tand(d), ...
                 '-', 'Color', [0 0 0 0.7], 'LineWidth', 0.7);
        end
    end
end
if any(isfinite(G.bed_z))
    plot(ax3, xkm, G.bed_z, 'r:', 'LineWidth', 1.6);
end

linkaxes([ax1 ax2 ax3],'xy');
xlim(ax1,[min(xkm) max(xkm)]); ylim(ax1,[0 max(G.z)]);
ax3.Position = ax2.Position;

xlabel(ax2,'distance (km)');
ylabel(ax2,'depth (m)');
cb2 = colorbar(ax3); cb2.Label.String = 'layer slope (deg)';
cb2.Ticks = local_ticks(o.clim_dip);
drawnow;
cb2.Position([1 3]) = cb1.Position([1 3]);
cb2.Position([2 4]) = [ax2.Position(2) ax2.Position(4)];

% No titles: the figure is meant to be dropped straight into a document,
% where the caption carries the provenance. Pass title_str to add one back.
if ~isempty(o.title_str)
    title(tl, o.title_str, 'FontWeight','bold', 'Interpreter','tex');
end

out_dir = fileparts(out_png);
if ~isempty(out_dir) && exist(out_dir,'dir') ~= 7
    mkdir(out_dir);
end
exportgraphics(fig, out_png, 'Resolution', o.dpi, 'BackgroundColor','w');
close(fig);

fig_file = out_png;
fprintf('  wrote %s\n', out_png);
end

% ------------------------------------------------------------------------
function t = local_ticks(lim)
% Round-number ticks, reaching both ends of the scale whenever a round step
% allows it. Left to itself MATLAB builds 0.3 as 0.30000000000000004, finds
% it past a limit of 0.3 and silently drops the top label.
raw = (lim(2) - lim(1))/6;
p10 = 10^floor(log10(raw));
nice = [1 2 2.5 5 10 20]*p10;
nice = nice(nice >= raw*(1 - 1e-9));
step = nice(1);
for k = 1:numel(nice)
    q = lim/nice(k);
    if all(abs(q - round(q)) < 1e-9)
        step = nice(k);
        break
    end
end
t = round(step*(ceil(lim(1)/step - 1e-9):floor(lim(2)/step + 1e-9)), 12);
end

function M = local_despeckle(A)
% Median of each finite cell and its finite 3 x 3 neighbours. Gaps stay
% gaps: nothing is filled, only isolated outliers are pulled in.
[n, m] = size(A);
P = nan(n+2, m+2);
P(2:end-1, 2:end-1) = A;
stack = nan(n, m, 9);
k = 0;
for di = 0:2
    for dj = 0:2
        k = k + 1;
        stack(:,:,k) = P(1+di:n+di, 1+dj:m+dj);
    end
end
M = median(stack, 3, 'omitnan');
M(~isfinite(A)) = NaN;
end
