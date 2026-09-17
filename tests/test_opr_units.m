function n_fail = test_opr_units()
% TEST_OPR_UNITS  End-to-end check of the OPR path on a synthetic echogram.
%
% Writes an OPR-shaped .mat with layers at a known dip, runs the whole
% driver over it, and checks that the dip comes back. This is the test that
% would have caught every unit bug in the original chain at once: the
% 20*log10 power conversion, the travel-time/distance mix-up in regrid, the
% all-NaN Bottom, and the inverted sign.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));
addpath(fullfile(here,'..','opr'));

n_fail = 0;
c = ls_cice();

tmp = tempname; mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp,'s'));

% ---- geometry ----------------------------------------------------------
ntrace = 500;
lat = -77.78 - (0:ntrace-1)*2.0e-5;        % a short straight line
lon = 158.75*ones(1,ntrace);
[px, py] = ls_polarstereo_fwd(lat, lon, 0);
dist = [0 cumsum(hypot(diff(px), diff(py)))];

dt = 1.6667e-9;
nt = 9000;                                  % ~760 m of ice
twtt = (0:nt-1)'*dt;
depth = twtt*c/2;

true_dip = 4;                               % deg, + = RISES with +x
lambda = 45;                                % m
bed_depth = 700;

[Dist, Depth] = meshgrid(dist, depth);
db = 12*sin(2*pi*(Depth + tand(true_dip)*Dist)/lambda) - 60;
db = db - 0.02*Depth;                        % a realistic depth trend
db(Depth > bed_depth) = -130;                % noise below the bed
Data = single(10.^(db/10));                  % OPR stores detected POWER

Time = twtt;
Latitude = lat; Longitude = lon;
Surface = zeros(1,ntrace);
Bottom = nan(1,ntrace);                      % as the real product ships it
Elevation = 2000*ones(1,ntrace);
GPS_time = (0:ntrace-1)*0.5;

data_file = fullfile(tmp,'Data_synth_01_001.mat');
save(data_file,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');

% ---- 1. the loader must not modify the input --------------------------
before = dir(data_file);
vars_before = sort(fieldnames(load(data_file)));
D = opr_load_echogram(data_file, struct('verbose',false));
after = dir(data_file);
vars_after = sort(fieldnames(load(data_file)));
if ~isequal(vars_before, vars_after) || before.bytes ~= after.bytes
    fprintf('  FAIL loader modified the input file\n'); n_fail = n_fail+1;
else
    fprintf('  ok   loader leaves the input file untouched\n');
end

% ---- 2. power -> dB is 10*log10 ---------------------------------------
G = opr_flatten_grid(D, struct('grid_spacing',2,'vert_exag',1, ...
    'z_max',bed_depth-25,'detrend_len',0,'smooth_len',0,'verbose',false));
mid = round(size(G.raw_db,2)/2);
expect = interp1(depth, db(:,round(ntrace/2)), G.z(:), 'linear');
got = G.raw_db(:,mid);
ok = isfinite(expect) & isfinite(got);
err = max(abs(expect(ok)-got(ok)));
if err > 2.5
    fprintf('  FAIL dB conversion off by up to %.2f dB\n', err); n_fail = n_fail+1;
else
    fprintf('  ok   dB conversion matches 10*log10 (max err %.2f dB)\n', err);
end

% ---- 3. the working grid really is isotropic --------------------------
if (G.x(2)-G.x(1)) ~= (G.z(2)-G.z(1))
    fprintf('  FAIL grid is not isotropic: dx %g, dz %g\n', ...
        G.x(2)-G.x(1), G.z(2)-G.z(1));
    n_fail = n_fail+1;
else
    fprintf('  ok   working grid is exactly isotropic (%g m)\n', G.x(2)-G.x(1));
end

% ---- 4. end-to-end dip recovery ---------------------------------------
R = RollingRadon_OPR(data_file, 'grid_spacing',2, 'vert_exag',1, ...
    'window_x',200, 'window_z',200, 'dip_max',20, 'dip_step',0.1, ...
    'z_max',bed_depth-25, 'z_pad_surface',40, 'detrend_len',40, ...
    'smooth_x',0, 'verbose',false);

v = R.slopes(isfinite(R.slopes));
if isempty(v)
    fprintf('  FAIL no windows solved\n'); n_fail = n_fail+1;
else
    med = median(v);
    frac = numel(v)/numel(R.slopes);
    if abs(med - true_dip) > 1.0
        fprintf('  FAIL median dip %.2f, expected %+g\n', med, true_dip);
        n_fail = n_fail+1;
    else
        fprintf('  ok   median dip %.2f (expected %+g), %.0f%% of windows solved\n', ...
            med, true_dip, 100*frac);
    end
end

% ---- 5. the opposite dip comes back with the opposite sign ------------
db2 = 12*sin(2*pi*(Depth - tand(true_dip)*Dist)/lambda) - 60 - 0.02*Depth;
db2(Depth > bed_depth) = -130;
Data = single(10.^(db2/10));
data_file2 = fullfile(tmp,'Data_synth_01_002.mat');
save(data_file2,'Data','Time','Latitude','Longitude','Surface','Bottom', ...
    'Elevation','GPS_time','-v7');
R2 = RollingRadon_OPR(data_file2, 'grid_spacing',2, 'vert_exag',1, ...
    'window_x',200, 'window_z',200, 'dip_max',20, 'dip_step',0.1, ...
    'z_max',bed_depth-25, 'z_pad_surface',40, 'detrend_len',40, ...
    'smooth_x',0, 'verbose',false);
v2 = R2.slopes(isfinite(R2.slopes));
if isempty(v2) || abs(median(v2) + true_dip) > 1.0
    fprintf('  FAIL mirrored dip: got %.2f, expected %+g\n', ...
        median(v2), -true_dip); n_fail = n_fail+1;
else
    fprintf('  ok   mirrored dip %.2f (expected %+g)\n', median(v2), -true_dip);
end

if n_fail == 0
    fprintf('test_opr_units: PASS\n');
else
    fprintf('test_opr_units: %d FAILURE(S)\n', n_fail);
end
end
