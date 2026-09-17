function n_fail = test_sign()
% TEST_SIGN  LS_RADON_DIP recovers the slope of a synthetic layered image.
%
% Convention under test:
%     POSITIVE slope = the layer RISES (gets shallower) with increasing x.
%
% This is d(elevation)/dx, the standard glaciological sense. Holschuh's two
% codebases disagree with each other here - the public SlopeExtraction_Radar
% release produces this sense, while NDH_MatlabTools applies `slopegrid*-1`
% at the end and produces the opposite - so it is worth pinning down
% explicitly rather than inheriting by accident.
%
% Each synthetic is sin(2*pi*(Z + tand(slope)*X)/lambda). Constant phase
% means Z = -tand(slope)*X + c, so depth DECREASES with x for a positive
% slope: the layer rises.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

n_fail = 0;
lambda = 40;                     % synthetic layer spacing (m)
gs = 2;
x = (0:299)*gs;
z = (0:199)*gs;
[X, Z] = meshgrid(x, z);

% ---- 1. isotropic grid -------------------------------------------------
for true_slope = [-8 -3 0 3 8]
    img = sin(2*pi*(Z + tand(true_slope)*X)/lambda);
    got = ls_radon_dip(img, gs, gs, 20, 0.1);
    n_fail = n_fail + check(sprintf('isotropic %+g', true_slope), ...
        got, true_slope, 0.35);
end

% ---- 2. anisotropic grid: the estimator must correct the aspect --------
dxa = 6; dza = 2;
xa = (0:99)*dxa;
za = (0:199)*dza;
[Xa, Za] = meshgrid(xa, za);
for true_slope = [-6 4]
    img = sin(2*pi*(Za + tand(true_slope)*Xa)/lambda);
    got = ls_radon_dip(img, dxa, dza, 20, 0.1);
    n_fail = n_fail + check(sprintf('anisotropic %+g', true_slope), ...
        got, true_slope, 1.0);
end

% ---- 3. sub-step refinement beats the search grid ----------------------
% The search step here is 0.5 deg but the answer should still land near the
% true 2.3 deg, because the criterion peak is interpolated.
img = sin(2*pi*(Z + tand(2.3)*X)/lambda);
got = ls_radon_dip(img, gs, gs, 20, 0.5);
n_fail = n_fail + check('refined (0.5 deg step)', got, 2.3, 0.3);

% ---- 4. a featureless window abstains rather than inventing an angle ---
got = ls_radon_dip(zeros(100,100), gs, gs, 20, 0.1);
if ~isscalar(got) || isfinite(got)
    fprintf('  FAIL flat window: expected NaN, got %s\n', mat2str(got));
    n_fail = n_fail + 1;
else
    fprintf('  ok   flat window abstains (NaN)\n');
end

% ---- 5. quality separates layering from noise --------------------------
% The criterion is the variance of the normalised projection. A max over
% the projection axis does NOT work here: over a few hundred candidate
% angles the largest random max sits several standard deviations up, so
% pure noise outscores real layering and the quality measure inverts.
rng(7);
[~, q_noise] = ls_radon_dip(randn(200,300), gs, gs, 20, 0.1);
[~, q_layer] = ls_radon_dip(sin(2*pi*(Z + tand(4)*X)/lambda), gs, gs, 20, 0.1);
if ~(q_layer > 5*q_noise)
    fprintf('  FAIL quality does not separate layering (%.1f) from noise (%.1f)\n', ...
        q_layer, q_noise);
    n_fail = n_fail + 1;
else
    fprintf('  ok   q separates layering (%.0f) from noise (%.1f)\n', ...
        q_layer, q_noise);
end

if n_fail == 0
    fprintf('test_sign: PASS\n');
else
    fprintf('test_sign: %d FAILURE(S)\n', n_fail);
end
end

function bad = check(name, got, want, tol)
if isscalar(got) && isfinite(got) && abs(got-want) <= tol
    fprintf('  ok   %-24s got %+7.3f (want %+g)\n', name, got, want);
    bad = 0;
else
    fprintf('  FAIL %-24s got %s (want %+g +/- %g)\n', ...
        name, mat2str(got), want, tol);
    bad = 1;
end
end
