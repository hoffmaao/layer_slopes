function n_fail = test_excluded_band()
% TEST_EXCLUDED_BAND  Tall windows span excluded bands without bias.
%
% A 400 m square window cannot avoid a band such as an accumulation
% radar's pulse-merge return, so a window is allowed to span one: the
% blanked samples are left out of the projections, their normalisation and
% the semblance, and the window abstains only when more than a set fraction
% of it is blanked. This checks that a band leaves the slope unchanged,
% that the fraction rule holds, and that blanking happens before the grid
% is filtered, so a bright band leaves no flat edge in the rows beside it.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));
n_fail = 0;

% ---- the estimator: an 80 x 1600 window, as the solver sees 400 m x 400 m
% Leaving a band out should cost about what losing that much data costs,
% so compare against the same window cut short by the same number of rows.
% One realisation is not enough to tell: an 80-sample-wide window has a
% broad angular peak, and any 10% of the data moves it by its own scatter.
% Measured: RMS 0.365 vs 0.287 deg in the image. The band costs a little
% more, because lines that cross it carry fewer samples, but at 20x
% exaggeration that is 0.018 vs 0.014 deg of true dip.
nx = 80; nz = 1600;
[X, Z] = meshgrid(0:nx-1, 0:nz-1);
e_band = []; e_short = [];
for app = [-3 0 3]                           % apparent slope (deg), rises
    for seed = 1:12
        rng(seed);
        w = sin(2*pi*(Z + tand(app)*X)/32) + 0.8*randn(nz, nx);
        wb = w; wb(700:859, :) = NaN;        % 10% of the window
        e_band(end+1) = ls_radon_dip(wb, 1, 1, 20, 0.1, true, 0.85) - app;           %#ok<AGROW>
        e_short(end+1) = ls_radon_dip(w(1:1440,:), 1, 1, 20, 0.1, true, 0.85) - app; %#ok<AGROW>
    end
end
rms = @(e) sqrt(mean(e.^2));
if rms(e_band) > 1.5*rms(e_short) || abs(mean(e_band)) > 0.25
    fprintf('  FAIL a 10%% band costs far more than 10%% of the data: RMS %.3f vs %.3f, bias %+.3f deg\n', ...
        rms(e_band), rms(e_short), mean(e_band));
    n_fail = n_fail + 1;
else
    fprintf('  ok   a 10%% band costs about what 10%% of the data costs: RMS %.3f vs %.3f deg\n', ...
        rms(e_band), rms(e_short));
end

rng(5);
win = sin(2*pi*(Z + tand(3)*X)/32) + 0.8*randn(nz, nx);
band = win; band(700:859, :) = NaN;
s2 = ls_radon_dip(band, 1, 1, 20, 0.1, true, 1);
wide = win; wide(600:999, :) = NaN;          % 25% of the window
s3 = ls_radon_dip(wide, 1, 1, 20, 0.1, true, 0.85);
if isfinite(s2) || isfinite(s3)
    fprintf('  FAIL windows over the excluded limit did not abstain\n');
    n_fail = n_fail + 1;
else
    fprintf('  ok   abstains when the band exceeds max_excluded\n');
end

% ---- the grid: a bright band must not leak into its neighbours ---------
c = ls_cice();
ntrace = 300; nt = 3000; dt = 1.6667e-9;
D.twtt = (0:nt-1)'*dt;
D.dist = (0:ntrace-1)*5;
D.surface_twtt = zeros(1, ntrace);
D.bed_twtt = nan(1, ntrace);
depth = D.twtt*c/2;
db = -60 - 0.05*depth + 3*randn(nt, ntrace);          % speckle, no layers
db(depth >= 70 & depth <= 88, :) = db(depth >= 70 & depth <= 88, :) + 25;
D.power = 10.^(db/10);
G = opr_flatten_grid(D, struct('grid_spacing',0.25,'vert_exag',20,'z_max',200, ...
    'detrend_len',30,'smooth_len',1.5,'smooth_x',60,'trace_balance',true, ...
    'exclude_z',[70 88],'verbose',false));
rowmean = abs(mean(G.img, 2, 'omitnan'));             % a flat feature = a large row mean
edge = (G.z > 55 & G.z < 70) | (G.z > 88 & G.z < 103);
far = G.z > 120 & G.z < 180;
ratio = max(rowmean(edge))/max(rowmean(far));
if ratio > 2
    fprintf('  FAIL rows beside the band carry its edge: %.1fx the far field\n', ratio);
    n_fail = n_fail + 1;
else
    fprintf('  ok   rows beside a bright band are clean (%.1fx the far field)\n', ratio);
end

if n_fail == 0
    fprintf('test_excluded_band: PASS\n');
else
    fprintf('test_excluded_band: %d FAILURE(S)\n', n_fail);
end
end
