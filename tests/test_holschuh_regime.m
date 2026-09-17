function n_fail = test_holschuh_regime()
% TEST_HOLSCHUH_REGIME  Reproduce the method in the regime it was built for.
%
% Holschuh et al. (2017) applied this to RDS / impulse-radar data: ~195 MHz
% centre frequency, narrow bandwidth, layers tens of metres apart, imaged
% between roughly 400 and 1700 m, with reflector slopes of +/-4 to +/-10 deg
% (their fig. 3). That is a completely different problem from a 600-900 MHz
% accumulation radar looking at 0.1 deg dips in the top 200 m.
%
% This runs a synthetic in Nick's regime through the solver using HIS
% published defaults - no vertical exaggeration, no along-track smoothing,
% no per-trace balancing, o_f 2/6, snr_thresh 2, vr 3, radon method 0 - to
% confirm that none of the fixes here have broken the method as published.
% If this passes, the implementation is faithful and any difficulty on
% accumulation-radar data is about the problem, not the code.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));
n_fail = 0;
cice_import
c = cice;

tmp = tempname; mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp,'s'));

% --- RDS-like acquisition ----------------------------------------------
B = 30e6;                       % narrow bandwidth
range_res = c/(2*B);            % ~2.8 m in ice
dt = 1/(4*B);                   % oversampled
nt = 6000;
twtt = (0:nt-1)'*dt;
depth = twtt*c/2;

ntrace = 1400;
lat = -79.0 - (0:ntrace-1)*1.35e-4;      % ~15 m trace spacing
lon = -110*ones(1,ntrace);
[px,py] = polarstereo_fwd(lat,lon,0);
dist = [0 cumsum(hypot(diff(px),diff(py)))];

lambda = 45;                    % layer spacing, metres - RDS scale
bed = 1700;

% --- a folded layer package, dips sweeping +/-6 deg ---------------------
[Dist, Depth] = meshgrid(dist, depth);
fold = 250*sin(2*pi*Dist/6000);              % gentle large-amplitude folds
true_dip_fn = atand(250*2*pi/6000*cos(2*pi*dist/6000));
db = 14*sin(2*pi*(Depth - fold)/lambda) - 55 - 0.022*Depth;
db(Depth < 150) = -40;                        % surface package
db(Depth > bed) = -125;                       % below the bed
Data = single(10.^(db/10));

Time = twtt; Latitude = lat; Longitude = lon;
Surface = zeros(1,ntrace);
Bottom = (2*bed/c)*ones(1,ntrace);
Elevation = 1500*ones(1,ntrace); GPS_time = (0:ntrace-1)*0.5;
f = fullfile(tmp,'Data_rds_01_001.mat');
save(f,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');

fprintf('  regime: %.1f m range resolution, %.0f m layer spacing, %.0f m traces\n', ...
    range_res, lambda, median(diff(dist)));
fprintf('  true dip sweeps %.1f to %.1f deg\n', min(true_dip_fn), max(true_dip_fn));

% --- Nick's published settings -----------------------------------------
R = RollingRadon_OPR(f, ...
    'grid_spacing', 5, ...        % isotropic, ~2 samples per range cell
    'vert_exag',    1, ...        % OFF - degree-scale dips need no help
    'window_x',     400, ...
    'window_z',     400, ...      % square window, as published
    'dip_max',      20, ...
    'dip_accept',   15, ...
    'z_pad_surface',200, ...
    'z_pad_bed',    50, ...
    'smooth_x',     0, ...        % OFF - would smear a 5 deg layer
    'smooth_len',   0, ...        % OFF
    'detrend_len',  0, ...        % OFF
    'trace_balance',false, ...    % OFF
    'solver_params', struct('o_f_horizontal',2,'o_f_vertical',6, ...
                            'snr_thresh',2,'snr_fac',1,'vr',3, ...
                            'radon_method',0), ...
    'verbose', false);

v = R.slopes(isfinite(R.slopes));
frac = numel(v)/numel(R.slopes);
fprintf('  solved %d/%d windows (%.0f%%), dip range %.1f to %.1f deg\n', ...
    numel(v), numel(R.slopes), 100*frac, min(v), max(v));

if frac < 0.30
    fprintf('  FAIL coverage %.0f%% - the method should work well here\n', 100*frac);
    n_fail = n_fail + 1;
else
    fprintf('  ok   coverage %.0f%%\n', 100*frac);
end

% Compare each solved cell against the true dip at its along-track position.
truth = interp1(dist, true_dip_fn, R.slope_x, 'linear', NaN);
err = [];
for i = 1:numel(R.slope_x)
    col = R.slopes(:,i);
    col = col(isfinite(col));
    if ~isempty(col) && isfinite(truth(i))
        err(end+1) = median(col) - truth(i); %#ok<AGROW>
    end
end
rms = sqrt(mean(err.^2));
fprintf('  dip error vs truth: median %.2f deg, RMS %.2f deg (n=%d)\n', ...
    median(err), rms, numel(err));

if rms > 1.5
    fprintf('  FAIL RMS dip error %.2f deg is too large\n', rms);
    n_fail = n_fail + 1;
else
    fprintf('  ok   RMS dip error %.2f deg\n', rms);
end

% The dip must actually track the fold, not just average to zero.
r = corr(err(:)*0 + truth(~isnan(truth)).', ...
         arrayfun(@(i) median(R.slopes(isfinite(R.slopes(:,i)),i)), ...
                  find(any(isfinite(R.slopes),1))).');
if ~isfinite(r) || r < 0.9
    fprintf('  FAIL recovered dip correlates r=%.2f with truth\n', r);
    n_fail = n_fail + 1;
else
    fprintf('  ok   recovered dip tracks the fold, r = %.3f\n', r);
end

if n_fail == 0
    fprintf('test_holschuh_regime: PASS\n');
else
    fprintf('test_holschuh_regime: %d FAILURE(S)\n', n_fail);
end
end
