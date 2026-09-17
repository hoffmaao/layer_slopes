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
%   alpha        slope overlay opacity, default 0.55
%   segments     also draw a dip tick in each solved cell, default true
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
o.alpha = 0.55;
o.segments = true;
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

if isempty(o.seg_len)
    if isfield(R.param,'window_x')
        o.seg_len = R.param.window_x/6;
    else
        o.seg_len = 20;
    end
end

% Rebuild the exact grid the solver saw.
D = opr_load_echogram(R.param.data_file, struct('verbose', false));
G = opr_flatten_grid(D, struct( ...
    'grid_spacing', R.param.grid_spacing, ...
    'z_pad_bed',    R.param.z_pad_bed, ...
    'z_max',        R.param.z_max, ...
    'detrend_len',  R.param.detrend_len, ...
    'smooth_len',   R.param.smooth_len, ...
    'agc_len',      R.param.agc_len, ...
    'trace_balance',R.param.trace_balance, ...
    'verbose',      false));

dB = G.raw_db;                       % 10*log10(power), surface-flattened
finite_db = dB(isfinite(dB));
db_lim = prctile(finite_db, [8 99.5]);

v = R.slopes(isfinite(R.slopes));
if isempty(o.clim_dip)
    if isempty(v)
        o.clim_dip = [-5 5];
    else
        m = max(1, prctile(abs(v), 98));
        o.clim_dip = [-m m];
    end
end

fig = figure('Visible','off','Color','w','Position',[80 80 1750 1080]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

% ---- panel 1: the power image ------------------------------------------
ax1 = nexttile(tl);
imagesc(ax1, G.x, G.z, dB);
colormap(ax1, gray); clim(ax1, db_lim);
set(ax1,'YDir','reverse','Layer','top','TickDir','out','Box','on');
hold(ax1,'on');
if any(isfinite(G.bed_z))
    plot(ax1, G.x, G.bed_z, 'r:', 'LineWidth', 1.6);
end
ylabel(ax1,'depth (m)');
set(ax1,'XTickLabel',[]);
cb1 = colorbar(ax1); cb1.Label.String = 'power (dB)';

% ---- panel 2: slope field over the power image -------------------------
ax2 = nexttile(tl);
imagesc(ax2, G.x, G.z, dB);
colormap(ax2, gray); clim(ax2, db_lim);
set(ax2,'YDir','reverse','Layer','top','TickDir','out','Box','on');
hold(ax2,'on');

% Overlay the dips in their own axes so the two colormaps coexist.
ax3 = axes('Position', ax2.Position, 'Color','none');
h = imagesc(ax3, R.slope_x, R.slope_z, R.slopes);
set(h, 'AlphaData', isfinite(R.slopes)*o.alpha);
colormap(ax3, b2r2(o.clim_dip(1), o.clim_dip(2)));
clim(ax3, o.clim_dip);
set(ax3,'YDir','reverse','Color','none','XTick',[],'YTick',[],'Box','off');
hold(ax3,'on');

if o.segments && ~isempty(v)
    for i = 1:numel(R.slope_x)
        for j = 1:numel(R.slope_z)
            d = R.slopes(j,i);
            if ~isfinite(d), continue; end
            plot(ax3, R.slope_x(i) + [-1 1]*o.seg_len, ...
                      R.slope_z(j) + [-1 1]*o.seg_len*tand(d), ...
                 '-', 'Color', [0 0 0 0.7], 'LineWidth', 0.7);
        end
    end
end
if any(isfinite(G.bed_z))
    plot(ax3, G.x, G.bed_z, 'r:', 'LineWidth', 1.6);
end

linkaxes([ax1 ax2 ax3],'xy');
xlim(ax1,[min(G.x) max(G.x)]); ylim(ax1,[0 max(G.z)]);
ax3.Position = ax2.Position;

xlabel(ax2,'distance (m)');
ylabel(ax2,'depth (m)');
cb2 = colorbar(ax3); cb2.Label.String = 'layer dip (deg)';
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
