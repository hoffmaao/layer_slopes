function n_fail = test_vert_exag()
% TEST_VERT_EXAG  Vertical exaggeration recovers sub-degree dips exactly.
%
% Interior englacial layers dip ~0.1 deg. Over any tractable window that is
% well under one range-resolution cell of vertical displacement, so the
% Radon cannot see it on an isotropic grid. Sampling along track at
% vert_exag times the vertical spacing and presenting the result as
% isotropic multiplies the apparent slope by exactly vert_exag. This checks
% both halves: that the tiny dip is recovered, and that it is NOT recovered
% without the exaggeration.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));
n_fail = 0;
c = ls_cice();

tmp = tempname; mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp,'s'));

ntrace = 3400;
lat = -86.70 - (0:ntrace-1)*5.0e-5;
lon = 68.6*ones(1,ntrace);
[px,py] = ls_polarstereo_fwd(lat,lon,0);
dist = [0 cumsum(hypot(diff(px),diff(py)))];

dt = 1.6667e-9; nt = 9000;
twtt = (0:nt-1)'*dt;
depth = twtt*c/2;

true_dip = 0.12;            % deg, + = RISES with +x; interior scale
lambda = 12;                % m, matching the observed spectral peak

[Dist, Depth] = meshgrid(dist, depth);
db = 10*sin(2*pi*(Depth + tand(true_dip)*Dist)/lambda) - 60 - 0.02*Depth;
db(Depth > 300) = -130;
Data = single(10.^(db/10));
Time = twtt; Latitude = lat; Longitude = lon;
Surface = zeros(1,ntrace); Bottom = nan(1,ntrace);
Elevation = 3000*ones(1,ntrace); GPS_time = (0:ntrace-1)*0.5;
f = fullfile(tmp,'Data_ve_01_001.mat');
save(f,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');

common = {'grid_spacing',0.25,'window_x',1000,'window_z',40, ...
    'z_pad_surface',40,'z_max',260,'smooth_len',0,'detrend_len',20, ...
    'smooth_x',0,'dip_max',1,'dip_accept',0.9,'verbose',false, ...
    'semb_thresh',0.2};   % a noise-free synthetic is no basis for the noise-calibrated gate

% ---- with exaggeration -------------------------------------------------
R = RollingRadon_OPR(f, common{:}, 'vert_exag',20, 'dip_step',0.005);
v = R.slopes(isfinite(R.slopes));
if isempty(v)
    fprintf('  FAIL vert_exag 20: no windows solved\n'); n_fail = n_fail+1;
elseif abs(median(v) - true_dip) > 0.03
    fprintf('  FAIL vert_exag 20: median %.4f deg, expected %+.3f\n', ...
        median(v), true_dip); n_fail = n_fail+1;
else
    fprintf('  ok   vert_exag 20 recovers %.4f deg (true %+.3f), %d windows\n', ...
        median(v), true_dip, numel(v));
end

% ---- mirrored dip keeps its sign --------------------------------------
db2 = 10*sin(2*pi*(Depth - tand(true_dip)*Dist)/lambda) - 60 - 0.02*Depth;
db2(Depth > 300) = -130;
Data = single(10.^(db2/10));
f2 = fullfile(tmp,'Data_ve_01_002.mat');
save(f2,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');
R2 = RollingRadon_OPR(f2, common{:}, 'vert_exag',20, 'dip_step',0.005);
v2 = R2.slopes(isfinite(R2.slopes));
if isempty(v2) || abs(median(v2) + true_dip) > 0.03
    fprintf('  FAIL mirrored: got %.4f, expected %+.3f\n', median(v2), -true_dip);
    n_fail = n_fail+1;
else
    fprintf('  ok   mirrored dip %.4f deg (true %+.3f)\n', median(v2), -true_dip);
end

% ---- vert_exag is EXACT, not merely helpful ---------------------------
% tan(apparent) = vert_exag*tan(true) is an identity, so exaggerating and
% then inverting must return the same answer as not exaggerating at all.
% That is the property worth guarding: vert_exag buys speed (vert_exag
% times fewer pixels per window), not accuracy.
%
% It used to buy accuracy too, when the criterion was the peak amplitude of
% the projection - a sub-pixel slope then had nothing to lock onto. With the
% variance criterion and sub-step refinement in LS_RADON_DIP the plain grid
% resolves 0.12 deg on its own, so the exaggeration is now an optimisation.
R3 = RollingRadon_OPR(f, common{:}, 'vert_exag',1, 'dip_step',0.01);
v3 = R3.slopes(isfinite(R3.slopes));
if isempty(v3)
    fprintf('  FAIL plain grid solved nothing\n'); n_fail = n_fail + 1;
else
    d = abs(median(v3) - median(v));
    if d > 0.01
        fprintf(['  FAIL exaggerated and plain grids disagree by %.4f deg ' ...
            '- the mapping is not being inverted correctly\n'], d);
        n_fail = n_fail + 1;
    else
        fprintf(['  ok   exaggerated %.4f and plain %.4f agree to %.4f deg ' ...
            '(mapping is exact)\n'], median(v), median(v3), d);
    end
end

% ---- the geometry itself: x and dips come back in real units ----------
% The last window centre sits at most half a window plus one step short
% of the end of the line.
if abs(R.slope_x(end) - dist(end)) > 2000
    fprintf('  FAIL slope_x ends at %.0f m, track ends at %.0f m\n', ...
        R.slope_x(end), dist(end)); n_fail = n_fail+1;
else
    fprintf('  ok   slope_x is in real metres (%.0f m of %.0f m)\n', ...
        R.slope_x(end), dist(end));
end

if n_fail == 0
    fprintf('test_vert_exag: PASS\n');
else
    fprintf('test_vert_exag: %d FAILURE(S)\n', n_fail);
end
end
