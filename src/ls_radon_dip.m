function [slope, q, crit, slope_axis, semb] = ls_radon_dip(win, dx, dz, slope_max, slope_step, refine, min_valid)
% LS_RADON_DIP  Dominant layer slope in one window, by Radon transform.
%
%   [slope, q, crit, slope_axis, semb] = LS_RADON_DIP(win, dx, dz, ...
%                              slope_max, slope_step, refine, min_valid)
%
%   win         [nz x nx] window, rows increasing in DEPTH
%   dx, dz      sample spacing along track and in depth, same units
%   slope_max   search +/- this many degrees
%   slope_step  angular step of the search (degrees)
%   refine      parabolic refinement of the peak (default true)
%   min_valid   fraction of the window that must be defined (not NaN) for
%               a result, default 1. Undefined samples - an excluded depth
%               band, the edge of the record - are left out of the
%               projections, their normalisation and the semblance.
%
%   slope       degrees. POSITIVE = the layer RISES (gets shallower) with
%               increasing x, i.e. d(elevation)/dx - the standard
%               glaciological sense.
%   q           peak criterion divided by the median criterion, kept as a
%               diagnostic. It compares the best slope with the other
%               CANDIDATES, so when the window cannot resolve slopes much
%               finer than the search range it stays small even over clean
%               layering. It is not a reliable quality measure; semb is.
%   crit        the criterion at every candidate slope
%   slope_axis  those candidate slopes
%   semb        semblance along the fitted slope: the fraction of the
%               window's energy that stacks coherently when every column is
%               shifted along that slope. 1 = perfectly continuous layering,
%               about 1/(independent columns) for noise. Unlike q it does not
%               depend on how wide the slope search is.
%
% METHOD. Projecting the window along a family of directions and taking the
% direction whose projection has the largest peak is the rolling Radon
% method of Holschuh et al. (2017, Geophys. Res. Lett. 44, 5561-5570). When
% the projection direction lies along the layering, every sample of a layer
% falls in the same bin and the projection has sharp maxima; off the layer
% orientation the energy smears out. This is an independent implementation
% of that idea.
%
% Two details matter for correctness:
%
%   * The Radon transform assumes square pixels, so a window whose samples
%     are not square is resampled first. Passing an already-isotropic grid
%     makes that a no-op, which is much cheaper.
%   * Longer chords through the window accumulate more samples, so the raw
%     projection favours whichever direction crosses the most pixels. The
%     criterion is normalised by the chord length so that orientation, not
%     window geometry, decides.
%
% See also LS_ROLLING_RADON

if nargin < 6 || isempty(refine), refine = true; end
if nargin < 7 || isempty(min_valid), min_valid = 1; end

slope_axis = -slope_max:slope_step:slope_max;
slope = NaN; q = NaN; semb = NaN;
crit = nan(size(slope_axis));

if ~ismatrix(win) || any(size(win) < 4) || ~any(isfinite(win(:)))
    return
end

% --- square pixels ------------------------------------------------------
% radon() works in pixel space, so anisotropic sampling would bias the
% angle. Resample to the finer of the two spacings when they differ.
if abs(dx - dz) > 1e-12*max(dx,dz)
    ds = min(dx, dz);
    nz = size(win,1); nx = size(win,2);
    zq = (0:ds:(nz-1)*dz);
    xq = (0:ds:(nx-1)*dx);
    [XQ, ZQ] = meshgrid(xq, zq);
    win = interp2((0:nx-1)*dx, ((0:nz-1)*dz).', double(win), XQ, ZQ, 'linear');
end

win = double(win);
valid = isfinite(win);
if mean(valid(:)) < min_valid - 1e-12
    % Too much of the window is undefined - an excluded band or the edge of
    % the record - for what is left to speak for it. Abstain.
    return
end

% --- normalise ----------------------------------------------------------
% Remove the mean so the transform responds to structure rather than to the
% window's overall brightness, then scale to unit range so the criterion is
% comparable between windows. Undefined samples are set to the mean, which
% is zero after this, so they add nothing to any projection.
win = win - mean(win(valid));
win(~valid) = 0;
s = max(abs(win(:)));
if s <= 0
    return
end
win = win/s;
win0 = win;
win0(~valid) = NaN;       % the semblance skips undefined samples outright

% --- taper to a disc ----------------------------------------------------
% Without this the window corners dominate the longest chords and the
% criterion picks up the shape of the rectangle.
[nz, nx] = size(win);
[XX, ZZ] = meshgrid(linspace(-1,1,nx), linspace(-1,1,nz));
rr = hypot(XX, ZZ);
taper = 0.5*(1 + cos(pi*min(1, max(0, (rr-0.7)/0.3))));   % flat, then cosine
taper(rr >= 1) = 0;
full = all(valid(:));
if ~full
    taper = taper.*valid;     % only defined samples count, in both below
end
win = win.*taper;

% --- transform ----------------------------------------------------------
% MATLAB's radon() measures theta from the x axis; a layer whose slope is
% `slope` degrees in the rising-positive sense projects coherently at
% theta = 90 + slope. tests/test_sign.m pins this down against synthetics.
theta = 90 + slope_axis;
[R, ~] = radon(win, theta);

% Chord-length normalisation: the projection of a uniform disc gives the
% number of samples contributing to each bin, so dividing by it removes the
% geometric bias toward directions that cross more of the window.
if full
    Rn = local_chord(taper, theta);
else
    Rn = radon(taper, theta);
end
Rn(Rn < 0.05*max(Rn(:))) = NaN;
Rc = R./Rn;

% Criterion: the VARIANCE of the normalised projection.
%
% The obvious choice is the projection's peak amplitude, but a max over the
% projection axis is itself a noisy statistic - over a few hundred candidate
% angles the largest random max is several standard deviations up, so pure
% noise scores HIGHER than real layering and the quality measure inverts.
% Variance averages over the whole projection axis instead: an aligned
% projection has sharp layer peaks and large variance, a misaligned one is
% smeared and flat. It is stable, and it discriminates the right way round.
crit = var(Rc, 0, 1, 'omitnan');
crit = crit(:).';

if all(~isfinite(crit))
    crit = nan(size(slope_axis));
    return
end

[peak, ia] = max(crit);
slope = slope_axis(ia);

% Quality: how much better the winning orientation scores than a typical
% one. Dimensionless, so it is comparable between windows and data sets.
% Incoherent noise has no preferred orientation and lands near 1.
good = isfinite(crit) & crit > 0;
m = median(crit(good));
if m > 0
    q = peak/m;
else
    q = 0;
end

% --- sub-step refinement ------------------------------------------------
% The search step is a hard quantum otherwise, and on interior layers the
% whole signal can be a few steps wide.
if refine && ia > 1 && ia < numel(crit)
    y1 = crit(ia-1); y2 = crit(ia); y3 = crit(ia+1);
    if all(isfinite([y1 y2 y3]))
        den = y1 - 2*y2 + y3;
        if den ~= 0
            d = 0.5*(y1 - y3)/den;
            if abs(d) <= 1
                slope = slope_axis(ia) + d*slope_step;
            end
        end
    end
end

if nargout > 4 && isfinite(slope)
    semb = local_semblance(win0, slope);
end
end

function S = local_semblance(win, slope)
% Shift every column along the slope, stack, and compare the energy of the
% stack with the energy of the columns. Rows are depth, so a layer that
% rises with +x sits at a smaller row index further right.
% Remove each column's mean first: a gain stripe is constant down a column,
% stacks coherently along ANY slope, and would otherwise count as layering.
win = win - mean(win, 1, 'omitnan');
[nz, nx] = size(win);
jc = (nx+1)/2;
r = (1:nz)';
acc = zeros(nz,1); en = zeros(nz,1); cnt = zeros(nz,1);
for j = 1:nx
    col = interp1(r, win(:,j), r - (j-jc)*tand(slope), 'linear', NaN);
    ok = isfinite(col);
    acc(ok) = acc(ok) + col(ok);
    en(ok) = en(ok) + col(ok).^2;
    cnt(ok) = cnt(ok) + 1;
end
use = cnt >= 0.8*nx;
if nnz(use) < 3 || sum(en(use)) <= 0
    S = NaN;
else
    S = sum(acc(use).^2 ./ cnt(use)) / sum(en(use));
end
end

function Rn = local_chord(taper, theta)
% The chord-length normalisation of a fully defined window depends only on
% the window size and the slope search, which are the same for every
% window of a run, so compute it once. It is one of the two transforms per
% window, and on a tall window the more expensive one to repeat.
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
key = sprintf('%d_%d_%.10g_%.10g_%d', size(taper,1), size(taper,2), ...
    theta(1), theta(end), numel(theta));
if isKey(cache, key)
    Rn = cache(key);
else
    Rn = radon(taper, theta);
    cache(key) = Rn;
end
end
