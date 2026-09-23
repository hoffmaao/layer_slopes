function n_fail = test_null_gate()
% TEST_NULL_GATE  The default quality gate admits layering and rejects noise.
%
% The gate is semblance along the fitted slope, thresholded per depth
% against the same echogram with each trace jittered in depth (layers no
% longer line up, everything else kept). This checks the guarantee that
% buys: on a profile whose upper part is layered and lower part is pure
% speckle, the speckle is admitted at about the requested false-alarm
% rate, the layering mostly passes, and the dip it returns is right. It is
% checked at two window lengths: the Radon peak/median ratio q that this
% gate replaced admitted almost no layering at window_x 500.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));
n_fail = 0;
c = ls_cice();

tmp = tempname; mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp,'s'));

% ---- an accumulation-radar-like profile --------------------------------
rng(7);
ntrace = 2000;                               % ~10 km at ~5 m
lat = -86.70 - (0:ntrace-1)*4.5e-5;
lon = 68.6*ones(1,ntrace);
[px,py] = ls_polarstereo_fwd(lat,lon,0);
dist = [0 cumsum(hypot(diff(px),diff(py)))];

dt = 1.6667e-9; nt = 4200;                   % ~355 m
twtt = (0:nt-1)'*dt;
depth = twtt*c/2;

true_dip = 0.10;                             % deg, + = RISES with +x
z_layer = 180;                               % layering above, speckle below

% Broadband stratigraphy: a random reflectivity series in depth, smoothed
% to the 0.53 m range resolution, displaced by the dip along track.
zr = (-50:0.05:450)';
refl = conv(randn(size(zr)), gausswin(21)/sum(gausswin(21)), 'same');
refl = refl/std(refl);
[Dist, Depth] = meshgrid(dist, depth);
lay = interp1(zr, refl, Depth + tand(true_dip)*Dist, 'linear', 0);
amp = 10.^((8*lay - 0.05*Depth)/20);
amp(Depth > z_layer) = 10.^(-0.05*Depth(Depth > z_layer)/20);
speckle = abs(randn(size(amp)) + 1i*randn(size(amp))).^2/2;
Data = single(amp.^2 .* speckle);

Time = twtt; Latitude = lat; Longitude = lon;
Surface = zeros(1,ntrace); Bottom = nan(1,ntrace);
Elevation = 3000*ones(1,ntrace); GPS_time = (0:ntrace-1)*0.5;
f = fullfile(tmp,'Data_gate_01_001.mat');
save(f,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');

for wx = [500 1000]
    R = RollingRadon_OPR(f, 'grid_spacing',0.25, 'vert_exag',20, ...
        'window_x',wx, 'window_z',20, 'z_pad_surface',30, 'z_max',340, ...
        'smooth_x',60, 'detrend_len',30, 'dip_max',1, 'verbose',false);
    hz = R.param.window_z/2;
    top = R.slope_z + hz <= z_layer;
    % Speckle windows start one detrend length below the boundary. The
    % boundary itself is a sharp, perfectly flat amplitude step, which is a
    % genuinely coherent horizontal feature: windows touching it are
    % admitted, at its own dip of zero, and are not false alarms.
    bot = R.slope_z - hz >= z_layer + 30;
    inice = R.status ~= 1;

    s_bot = R.slopes(bot,:); in_bot = inice(bot,:);
    fa = nnz(isfinite(s_bot))/nnz(in_bot);
    s_top = R.slopes(top,:); in_top = inice(top,:);
    cov = nnz(isfinite(s_top))/nnz(in_top);
    v = s_top(isfinite(s_top));

    fprintf('  window_x %4d: median gate %.2f, speckle admitted %.1f%%, layering solved %.0f%%\n', ...
        wx, median(R.param.semb_thresh(isfinite(R.param.semb_thresh))), 100*fa, 100*cov);
    if fa > 0.03
        fprintf('  FAIL speckle admitted at %.1f%% (want <= 3%%)\n', 100*fa);
        n_fail = n_fail + 1;
    else
        fprintf('  ok   speckle held to the false-alarm rate\n');
    end
    if cov < 0.5
        fprintf('  FAIL only %.0f%% of layered windows solved\n', 100*cov);
        n_fail = n_fail + 1;
    else
        fprintf('  ok   layering admitted\n');
    end
    if isempty(v) || abs(median(v) - true_dip) > 0.02
        fprintf('  FAIL median dip %.3f, want %.3f\n', median(v), true_dip);
        n_fail = n_fail + 1;
    else
        fprintf('  ok   median dip %.3f (want %.3f)\n', median(v), true_dip);
    end
end

if n_fail == 0
    fprintf('test_null_gate: PASS\n');
else
    fprintf('test_null_gate: %d FAILURE(S)\n', n_fail);
end
end
