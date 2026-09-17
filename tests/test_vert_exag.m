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
cice_import
c = cice;

tmp = tempname; mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp,'s'));

ntrace = 3400;
lat = -86.70 - (0:ntrace-1)*5.0e-5;
lon = 68.6*ones(1,ntrace);
[px,py] = polarstereo_fwd(lat,lon,0);
dist = [0 cumsum(hypot(diff(px),diff(py)))];

dt = 1.6667e-9; nt = 9000;
twtt = (0:nt-1)'*dt;
depth = twtt*c/2;

true_dip = 0.12;            % deg - the real scale of interior layer dip
lambda = 12;                % m, matching the observed spectral peak

[Dist, Depth] = meshgrid(dist, depth);
db = 10*sin(2*pi*(Depth - tand(true_dip)*Dist)/lambda) - 60 - 0.02*Depth;
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
    'dip_max',1,'dip_accept',0.9,'verbose',false};

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
db2 = 10*sin(2*pi*(Depth + tand(true_dip)*Dist)/lambda) - 60 - 0.02*Depth;
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

% ---- control: without exaggeration the answer is quantised -------------
% The Radon searches on a fixed angular grid. On an isotropic grid that
% step IS the measurement for a 0.12 deg layer, so every window returns a
% multiple of it; the exaggerated grid resolves inside one step. Compare
% the quantisation, not the median - a coarse grid can still land near the
% right answer by luck, which is exactly what makes it misleading.
R3 = RollingRadon_OPR(f, common{:}, 'vert_exag',1);
v3 = R3.slopes(isfinite(R3.slopes));
q_plain = numel(unique(round(v3,4)));
q_exag  = numel(unique(round(v,4)));
if isempty(v3)
    fprintf('  ok   control: vert_exag 1 solves nothing\n');
elseif q_exag <= q_plain
    fprintf(['  FAIL control: exaggerated grid gives %d distinct dips, ' ...
        'plain grid %d - no gain in resolution\n'], q_exag, q_plain);
    n_fail = n_fail+1;
else
    fprintf(['  ok   control: plain grid returns %d distinct dip value(s) ' ...
        '(quantised at %.3f deg); exaggerated returns %d\n'], ...
        q_plain, 0.1, q_exag);
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
