% VALIDATE_HORIZON  Independent ground truth tracked from the image itself.
%
% Andrew's objection: between 6 and 12 km the layering visibly slants upward
% with distance, so the dip should be clearly negative there, and the solver
% reports something much smaller. This checks it directly.
%
% A slope measurement is taken straight off the echogram with no Radon
% involved: seed on a reflector, follow it trace by trace under a continuity
% constraint, smooth, and differentiate. The tracked pick is written out as
% a figure so it can be confirmed by eye BEFORE it is used as truth - a
% tracker that jumps between reflectors produces confident nonsense.
%
% Then the solver is run in several configurations and compared at the
% horizon's own depth and distance.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));
out_dir = '/kucresis/scratch/hoffmana_sta/layer_slopes/products/diag';
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end

data_file = ['/kucresis/scratch/dataproducts/opr_data/accum/' ...
    '2024_Antarctica_Ground2/CSARP_post/CSARP_standard/20250108_02/' ...
    'Data_20250108_02_005.mat'];

if ~exist('seed_x','var'), seed_x = 10000; end   % m along track
if ~exist('seed_z','var'), seed_z = 135;   end   % m depth, on the horizon

D = opr_load_echogram(data_file, struct('verbose',false));
G = opr_flatten_grid(D, struct('grid_spacing',0.25,'z_max',220, ...
    'smooth_len',0,'detrend_len',0,'agc_len',0,'trace_balance',false, ...
    'vert_exag',1,'verbose',false));

dz = G.z(2)-G.z(1);
dx = G.x(2)-G.x(1);

% --- condition for TRACKING only (not for the solver) -------------------
% Remove the depth power envelope so a reflector is a local maximum rather
% than merely shallower than its neighbours, and smooth along track so the
% pick follows the layer instead of speckle.
A = G.raw_db;
A = A - movmean(A, round(30/dz), 1, 'omitnan');
A = movmean(A, round(200/dx), 2, 'omitnan');
A(~isfinite(A)) = -inf;

% --- seed, then propagate under a continuity constraint -----------------
[~, c0] = min(abs(G.x - seed_x));
zwin = 6;                                   % m, max step between traces
[~, r0] = max(A(abs(G.z - seed_z) <= 8, c0));
rlo = find(abs(G.z - seed_z) <= 8, 1);
r0 = rlo + r0 - 1;

nx = numel(G.x);
ridx = nan(1, nx);
ridx(c0) = r0;
step = round(zwin/dz);
for c = c0+1:nx                             % forwards
    lo = max(1, ridx(c-1)-step); hi = min(numel(G.z), ridx(c-1)+step);
    [~, k] = max(A(lo:hi, c)); ridx(c) = lo + k - 1;
end
for c = c0-1:-1:1                           % backwards
    lo = max(1, ridx(c+1)-step); hi = min(numel(G.z), ridx(c+1)+step);
    [~, k] = max(A(lo:hi, c)); ridx(c) = lo + k - 1;
end

zpk = G.z(ridx);
zpk_s = movmean(zpk, round(500/dx));        % the layer is smooth; the pick is not
dip_true = -atand(gradient(zpk_s, dx));     % + = RISES with +x

% --- look at it before believing it -------------------------------------
f = figure('Visible','off','Color','w','Position',[60 60 1700 800]);
ax = axes(f);
imagesc(ax, G.x/1000, G.z, G.raw_db); colormap(ax, gray);
clim(ax, prctile(G.raw_db(isfinite(G.raw_db)), [8 99.5]));
set(ax,'YDir','reverse'); hold(ax,'on');
plot(ax, G.x/1000, zpk_s, 'r-', 'LineWidth', 1.6);
plot(ax, seed_x/1000, seed_z, 'co', 'MarkerSize', 10, 'LineWidth', 2);
xlabel(ax,'distance (km)'); ylabel(ax,'depth (m)');
exportgraphics(f, fullfile(out_dir,'tracked_horizon.png'), 'Resolution',110);
close(f);

fprintf('tracked horizon: %.0f m at 0 km -> %.0f m at %.0f km\n', ...
    zpk_s(1), zpk_s(end), G.x(end)/1000);
seg = G.x >= 6000 & G.x <= 12000;
fprintf('  over 6-12 km: %+.1f m change, mean dip %+.3f deg (%.3f to %.3f)\n', ...
    zpk_s(find(seg,1,'last')) - zpk_s(find(seg,1)), ...
    mean(dip_true(seg)), min(dip_true(seg)), max(dip_true(seg)));
fprintf('  whole line  : mean dip %+.3f deg, range %.3f to %.3f\n', ...
    mean(dip_true), min(dip_true), max(dip_true));
fprintf('  wrote %s - CHECK THE PICK BEFORE TRUSTING THE TABLE\n\n', ...
    fullfile(out_dir,'tracked_horizon.png'));

% --- solver configurations ----------------------------------------------
base = {'grid_spacing',0.25,'vert_exag',20,'window_x',2000,'window_z',20, ...
    'dip_max',1,'dip_accept',0.9,'dip_step',0.005, ...
    'z_pad_surface',30,'z_max',200,'exclude_z',[70 88], ...
    'smooth_len',1.5, ...
    'solver_params',struct('o_f_horizontal',4,'o_f_vertical',4, ...
                           'snr_thresh',2,'snr_fac',1,'vr',1), ...
    'verbose',false};

cases = { ...
 'sx250  hp0   tb on  ',      250,       0,        true
 'sx250  hp20  tb on  ',      250,      20,        true
 'sx60   hp20  tb on  ',       60,      20,        true
 'sx0    hp20  tb on  ',        0,      20,        true
 'sx0    hp20  tb off ',        0,      20,        false
 'sx0    hp0   tb off ',        0,       0,        false };

fprintf('%-22s %-7s %-9s %-9s %-9s %-7s %s\n', ...
    'configuration','n','solver','truth','bias','slope','r');
for i = 1:size(cases,1)
    evalc(['R = RollingRadon_OPR(data_file, base{:}, ' ...
      '''smooth_x'',cases{i,2}, ''detrend_len'',cases{i,3}, ' ...
      '''trace_balance'',cases{i,4});']);
    sv = []; tv = [];
    for c = 1:numel(R.slope_x)
        zt = interp1(G.x, zpk_s,   R.slope_x(c), 'linear', NaN);
        dt = interp1(G.x, dip_true, R.slope_x(c), 'linear', NaN);
        if ~isfinite(zt) || ~isfinite(dt), continue; end
        [~, r] = min(abs(R.slope_z - zt));
        if abs(R.slope_z(r) - zt) > 15, continue; end
        if isfinite(R.slopes(r,c))
            sv(end+1) = R.slopes(r,c); tv(end+1) = dt; %#ok<AGROW>
        end
    end
    if numel(sv) < 5
        fprintf('%-22s %-7d  (too few matches)\n', cases{i,1}, numel(sv));
        continue
    end
    p = polyfit(tv(:), sv(:), 1);           % gain of solver against truth
    fprintf('%-22s %-7d %-9.3f %-9.3f %-9.3f %-7.2f %.2f\n', ...
        cases{i,1}, numel(sv), median(sv), median(tv), ...
        median(sv)-median(tv), p(1), corr(sv(:), tv(:)));
end
