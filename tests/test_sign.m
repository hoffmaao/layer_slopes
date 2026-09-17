function n_fail = test_sign()
% TEST_SIGN  radon_ndh recovers the dip of a synthetic layered image.
%
% Convention under test:
%     POSITIVE dip = the layer gets DEEPER with increasing x.
%
% The public SlopeExtraction_Radar release returned the negative of this
% (Nick's full NDH_MatlabTools version corrects it with `slopegrid*-1` at
% the end of RollingRadon; the public release dropped that line). These
% cases pin the convention down so it cannot regress silently.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

n_fail = 0;
tol = 0.35;                      % deg; the search grid is 0.1 deg
lambda = 40;                     % synthetic layer spacing (m)

% ---- 1. isotropic grid, depth increasing downward ----------------------
gs = 2;
x = (0:299)*gs;
z = (0:199)*gs;
[X, Z] = meshgrid(x, z);
for true_dip = [-8 -3 0 3 8]
    img = sin(2*pi*(Z - tand(true_dip)*X)/lambda);
    got = radon_ndh(x, z, img, 20, 0, 0);
    n_fail = n_fail + check(sprintf('isotropic dip %+g', true_dip), ...
        got, true_dip, tol);
end

% ---- 2. anisotropic grid (radon_ndh must correct the aspect) -----------
dx = 6; dz = 2;
x = (0:99)*dx;
z = (0:199)*dz;
[X, Z] = meshgrid(x, z);
for true_dip = [-6 4]
    img = sin(2*pi*(Z - tand(true_dip)*X)/lambda);
    got = radon_ndh(x, z, img, 20, 0, 0);
    n_fail = n_fail + check(sprintf('anisotropic dip %+g', true_dip), ...
        got, true_dip, 1.0);
end

% ---- 3. y axis decreasing (the time_increases_downward == 0 branch) ----
gs = 2;
x = (0:299)*gs;
z = (199:-1:0)*gs;              % decreasing
[X, Z] = meshgrid(x, z);
for true_dip = [-5 5]
    % Z still means depth, so the same expression describes the same
    % physical layer; only the axis ordering changed.
    img = sin(2*pi*(Z - tand(true_dip)*X)/lambda);
    got = radon_ndh(x, z, img, 20, 0, 0);
    n_fail = n_fail + check(sprintf('y-decreasing dip %+g', true_dip), ...
        got, true_dip, tol);
end

% ---- 4. a scalar is returned even when angles tie ----------------------
img = zeros(100,100);
got = radon_ndh((0:99)*2, (0:99)*2, img, 20, 0, 4);
if ~isscalar(got)
    fprintf('  FAIL tie handling: expected a scalar, got %s\n', mat2str(size(got)));
    n_fail = n_fail + 1;
else
    fprintf('  ok   tie handling returns a scalar\n');
end

if n_fail == 0
    fprintf('test_sign: PASS\n');
else
    fprintf('test_sign: %d FAILURE(S)\n', n_fail);
end
end

function bad = check(name, got, want, tol)
if isscalar(got) && abs(got-want) <= tol
    fprintf('  ok   %-26s got %+7.2f (want %+g)\n', name, got, want);
    bad = 0;
else
    fprintf('  FAIL %-26s got %s (want %+g +/- %g)\n', ...
        name, mat2str(got), want, tol);
    bad = 1;
end
end
