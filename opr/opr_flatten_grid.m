function G = opr_flatten_grid(D, opt)
% OPR_FLATTEN_GRID  Surface-flatten an echogram onto an isotropic depth grid.
%
%   G = OPR_FLATTEN_GRID(D, opt)   D from OPR_LOAD_ECHOGRAM
%
% opt fields
%   .grid_spacing  isotropic grid spacing (m), default 2
%   .z_pad_bed     stop this far above the bed (m), default 25
%   .z_max         hard depth cap (m), [] = derive from the bed pick
%   .bed_default   assumed ice thickness (m) when no bed pick exists,
%                  default 1500
%   .smooth_len    low-pass length along depth (m), 0 = off, default 1.5
%   .detrend_len   high-pass length along depth (m): a local quadratic fit
%                  over this length is subtracted. 0 = off, default 15
%   .trace_balance equalise traces against each other, default true
%   .smooth_x      along-track low-pass length (m), 0 = off. Suppresses the
%                  vertical striping left by trace-to-trace gain changes.
%   .exclude_z     N x 2 array of [z_start z_end] depth bands (m) to blank,
%                  e.g. where the first and second pulse returns merge.
%   .vert_exag     along-track spacing = grid_spacing*vert_exag, presented
%                  to the solver as isotropic. Default 1.
%   .agc_len       running-RMS normalisation length (m), 0 = off, default 0
%   .c_ice         wave speed (m/s), default LS_CICE()
%   .verbose
%
% Returns
%   .x       [1 x nx] along-track distance (m), spacing grid_spacing
%   .z       [1 x nz] depth below the ice surface (m), spacing grid_spacing
%   .img     [nz x nx] log power (dB), depth-detrended
%   .raw_db  [nz x nx] log power (dB) before detrending
%   .bed_z   [1 x nx] bed depth below surface (m)
%
% Why this step exists at all:
%
%   * RADON NEEDS SQUARE PIXELS. LS_RADON_DIP resamples to a 1:1 aspect
%     internally, and on a native OPR grid (dz 0.14 m, dx 5.9 m) that means
%     upsampling every window 42x along track. Gridding once, isotropically
%     and up front, makes the aspect correction a no-op and is the
%     difference between a tractable run and an intractable one.
%   * CROP TO THE ICE COLUMN. The example frame carries 32066 samples
%     spanning ~4.5 km of range for ~1 km of ice. Everything below the bed
%     is noise, and gridding it is what made the original regrid call ask
%     for an impossible allocation.
%   * REMOVE THE DEPTH TREND. Radar power falls off steeply with depth and
%     that trend is horizontal. Left in place it is the strongest linear
%     feature in the window and biases every dip estimate toward zero.
%
% See also ROLLINGRADON_OPR, OPR_LOAD_ECHOGRAM

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

if nargin < 2, opt = struct(); end
if ~isfield(opt,'grid_spacing'), opt.grid_spacing = 2;   end
if ~isfield(opt,'z_pad_bed'),    opt.z_pad_bed = 25;     end
if ~isfield(opt,'z_max'),        opt.z_max = [];         end
if ~isfield(opt,'bed_default'),  opt.bed_default = 1500; end
if ~isfield(opt,'detrend_len'),  opt.detrend_len = 15;   end
if ~isfield(opt,'smooth_len'),   opt.smooth_len = 1.5;   end
if ~isfield(opt,'agc_len'),      opt.agc_len = 0;        end
if ~isfield(opt,'trace_balance'),opt.trace_balance = true; end
if ~isfield(opt,'vert_exag'),    opt.vert_exag = 1;      end
if ~isfield(opt,'smooth_x'),     opt.smooth_x = 0;       end
if ~isfield(opt,'exclude_z'),    opt.exclude_z = [];     end
if ~isfield(opt,'verbose'),      opt.verbose = true;     end
if ~isfield(opt,'c_ice') || isempty(opt.c_ice)
    opt.c_ice = ls_cice();
end

c = opt.c_ice;
gs = opt.grid_spacing;
dz = gs;                      % vertical sample spacing (m)
dx = gs*opt.vert_exag;        % along-track sample spacing (m)

surf = D.surface_twtt;
surf(~isfinite(surf)) = 0;

bed_z = (D.bed_twtt - surf)*c/2;          % bed depth below surface (m)

% --- depth extent -------------------------------------------------------
if ~isempty(opt.z_max)
    z_bot = opt.z_max;
elseif any(isfinite(bed_z))
    z_bot = max(bed_z(isfinite(bed_z))) - opt.z_pad_bed;
else
    % No bed pick anywhere. Assume a nominal ice thickness rather than
    % gridding the whole record - for this radar the record runs to ~4.5 km
    % of range against ice that is nowhere near that thick, and the surplus
    % is pure noise. bed_default is the knob for this.
    z_bot = opt.bed_default - opt.z_pad_bed;
    z_rec = (D.twtt(end) - median(surf))*c/2;
    if z_bot > z_rec
        z_bot = z_rec;
    end
    warning('opr_flatten_grid:noBed', ...
        ['No bed pick available; assuming %.0f m of ice and gridding to ' ...
         '%.0f m. Set bed_default (or z_max) if that is wrong for this ' ...
         'site.'], opt.bed_default, z_bot);
end
if z_bot <= dz
    error('opr_flatten_grid:emptyColumn', ...
        'Usable ice column is empty (bottom at %.1f m).', z_bot);
end

% VERTICAL EXAGGERATION. Interior englacial layers dip by ~0.1 deg, and
% over any window short enough to be useful that is far less than one
% range-resolution cell of vertical displacement - unmeasurable. Sampling
% along track at vert_exag times the vertical spacing and then presenting
% the grid to the solver AS IF it were isotropic multiplies the apparent
% slope by vert_exag, which moves a 0.1 deg dip into the range the Radon
% transform can actually resolve, and costs vert_exag times fewer pixels
% per window. The mapping is exact:
%
%       tan(apparent) = vert_exag * tan(true)
%
% so nothing is approximated; RollingRadon_OPR inverts it on the way out.
% vert_exag = 1 reproduces the plain isotropic grid.
nz = floor(z_bot/dz) + 1;
nx = floor(D.dist(end)/dx) + 1;
z = (0:nz-1)*dz;
x = (0:nx-1)*dx;              % real along-track distance (m)
x_solver = (0:nx-1)*dz;       % same columns, presented as isotropic

% --- power -> dB --------------------------------------------------------
pw = D.power;
pw(~isfinite(pw) | pw <= 0) = NaN;
if ~any(isfinite(pw(:)))
    error('opr_flatten_grid:noPower','Data contains no positive samples.');
end
pw(~isfinite(pw)) = min(pw(isfinite(pw)));
db = 10*log10(pw);          % OPR standard/qlook Data is detected POWER

% --- anti-alias before vertical decimation ------------------------------
dz_native = median(diff(D.twtt))*c/2;
decim = dz/dz_native;
if decim > 1.5
    db = movmean(db, max(2, round(decim)), 1);
end

% --- resample each trace onto depth-below-surface ------------------------
nx_in = numel(D.dist);
db_z = zeros(nz, nx_in);
for n = 1:nx_in
    db_z(:,n) = interp1(D.twtt, db(:,n), surf(n) + 2*z(:)/c, 'linear', NaN);
end

% --- resample along track -----------------------------------------------
img = interp1(D.dist(:), db_z.', x(:), 'linear', NaN).';

% Carry the nearest valid trace into any all-NaN edge column so a single
% bad column does not poison every window that overlaps it.
bad = all(~isfinite(img),1);
if any(bad) && ~all(bad)
    good = find(~bad);
    for k = find(bad)
        [~, j] = min(abs(good-k));
        img(:,k) = img(:,good(j));
    end
end

G = struct();
G.x = x;
G.x_solver = x_solver;
G.z = z;
G.dx = dx;
G.dz = dz;
G.vert_exag = opt.vert_exag;
G.grid_spacing = gs;
G.dz_native = dz_native;
G.raw_db = img;
G.c_ice = c;

% --- band-pass along depth, then AGC ------------------------------------
% Two things in a radargram are much stronger linear features than the
% layering, and the Radon transform will happily lock onto either:
%
%   * the power-vs-depth envelope, which is horizontal, and biases every
%     dip toward zero;
%   * trace-to-trace gain and coupling variation, which shows up as
%     vertical striping and biases dips toward steep values.
%
% The band-pass keeps only the depth wavelengths the layering actually
% occupies, and the AGC divides out the running amplitude so a bright
% trace does not outvote its neighbours. Without these the solver measures
% the gain structure instead of the stratigraphy.
work = img;

% EXCLUDED DEPTH BANDS. Some depth bands carry strong horizontal features
% that are not stratigraphy: the surface ringdown of a ground-based system,
% and the band where an accumulation radar's first and second pulse returns
% merge. They are blanked to NaN BEFORE any filtering, and every filter
% below skips NaN, so a band cannot leak into its neighbours: filtered
% across it, the depth high-pass subtracts the band's brightness from the
% rows beside it and leaves a flat edge that a tall window would fit.
% Windows count only defined samples, and abstain when more than
% max_excluded of them are blanked (see ROLLINGRADON_OPR).
G.exclude_z = opt.exclude_z;
blank = false(nz, 1);
if ~isempty(opt.exclude_z)
    ez = opt.exclude_z;
    if size(ez,2) ~= 2
        error('opr_flatten_grid:badExclude', ...
            'exclude_z must be an N x 2 array of [z_start z_end] in metres.');
    end
    for k = 1:size(ez,1)
        blank = blank | (z(:) >= min(ez(k,:)) & z(:) <= max(ez(k,:)));
    end
    work(blank,:) = NaN;
    nmask = nnz(blank);
    if opt.verbose
        fprintf('    excluded %d of %d depth samples (%s m)\n', ...
            nmask, nz, strjoin(arrayfun(@(k) sprintf('%g-%g', ...
            ez(k,1), ez(k,2)), 1:size(ez,1), 'UniformOutput', false), ', '));
    end
end

% low-pass: drop sub-resolution noise
if opt.smooth_len > 0
    nlo = max(1, round(opt.smooth_len/dz));
    if nlo > 1
        work = movmean(work, nlo, 1, 'omitnan');
        work(blank,:) = NaN;      % omitnan fills a band's edges; keep it exact
    end
    G.smooth_samples = nlo;
else
    G.smooth_samples = 0;
end

% high-pass: remove the long-wavelength power envelope
if opt.detrend_len > 0
    nsm = max(3, round(opt.detrend_len/dz));
    if mod(nsm,2) == 0, nsm = nsm+1; end
    work = work - local_quadratic_trend(work, nsm);
    G.detrend_samples = nsm;
else
    G.detrend_samples = 0;
end

% Trace balance: one scalar per trace, so profile-to-profile gain and
% coupling swings stop dominating the Radon. This is deliberately NOT a
% running AGC - a running window also flattens the amplitude variation
% DOWN each trace, and RollingRadon's SNR gate is computed from exactly
% that variation, so a full AGC silently rejects every window. Balancing
% traces against each other fixes the striping and leaves the gate intact.
if opt.trace_balance
    env = median(abs(work), 1, 'omitnan');           % [1 x nx]
    ref = median(env(isfinite(env) & env > 0));
    if ~isfinite(ref) || ref <= 0, ref = 1; end
    env(~isfinite(env) | env < 1e-3*ref) = 1e-3*ref;
    work = work .* (ref./env);
end

% ALONG-TRACK SMOOTHING. Trace-to-trace gain and coupling changes appear
% as vertical stripes, and a stripe is a strong linear feature that the
% Radon will happily fit. Smoothing along track removes them, and at these
% dips it is close to free: a layer dipping 0.1 deg moves 0.05 m across
% 30 m of track and 0.17 m across 100 m, both well inside one 0.53 m range
% cell, so the layer geometry the solver is trying to measure is untouched
% while the stripes are averaged away.
if opt.smooth_x > 0
    nxs = max(1, round(opt.smooth_x/dx));
    if mod(nxs,2) == 0, nxs = nxs+1; end
    if nxs > 1
        work = movmean(work, nxs, 2, 'omitnan');
    end
    G.smooth_x_samples = nxs;
else
    G.smooth_x_samples = 0;
end

% Optional running AGC, off by default. If you switch it on, drop
% solver_params.snr_thresh to match - the gate is amplitude-based.
if opt.agc_len > 0
    nagc = max(3, round(opt.agc_len/dz));
    if mod(nagc,2) == 0, nagc = nagc+1; end
    env = movmean(abs(work), nagc, 1, 'omitnan');
    ref = median(env(isfinite(env) & env > 0));
    if ~isfinite(ref) || ref <= 0, ref = 1; end
    env(~isfinite(env) | env < 1e-6*ref) = 1e-6*ref;
    work = work.*(ref./env);
    G.agc_samples = nagc;
else
    G.agc_samples = 0;
end

work(blank,:) = NaN;
G.img = work;

G.bed_z = interp1(D.dist(:), bed_z(:), x(:), 'linear', NaN).';
G.surface_z = zeros(1,nx);

if opt.verbose
    fprintf('    gridded %d x %d (dz %.2f m, dx %.2f m, vert exag %gx), depth 0-%.0f m, %.0f m along track\n', ...
        nz, nx, dz, dx, opt.vert_exag, z(end), x(end));
end
end

% ------------------------------------------------------------------------
function T = local_quadratic_trend(A, n)
% Least-squares quadratic through each sample's n-sample neighbourhood down
% the trace, evaluated at that sample; undefined samples are left out, and
% at the top and bottom of the grid or beside an excluded band the fit is
% one-sided. A moving mean is biased wherever the neighbourhood is
% one-sided, and a moving line cannot follow the curvature of the power
% envelope - the bright firn falls away by ~25 dB within 50 m near 150-200 m
% on this radar - and both leave a flat residual. A tall window reads that
% residual as layering, and so does the noise copy the gate is set from,
% which then puts the gate above the real data.
%
% Moments are taken about each sample (d = offset in samples), so the
% 3 x 3 normal equations stay well conditioned however deep the grid.
h = floor(n/2);
d = (-h:h).';
ok = isfinite(A);
w = double(ok);
y = A;  y(~ok) = 0;
mom = @(v, k) conv2(v, flipud(d.^k), 'same');     % sum over the window of v.*d.^k
m0 = mom(w, 0); m1 = mom(w, 1); m2 = mom(w, 2); m3 = mom(w, 3); m4 = mom(w, 4);
b0 = mom(y, 0); b1 = mom(y, 1); b2 = mom(y, 2);
det3 = @(a,b,c, e,f,g, i,j,k) a.*(f.*k - g.*j) - b.*(e.*k - g.*i) + c.*(e.*j - f.*i);
Dm = det3(m0,m1,m2, m1,m2,m3, m2,m3,m4);
T = det3(b0,m1,m2, b1,m2,m3, b2,m3,m4)./Dm;         % the fit at d = 0
flat = m0 < 5 | ~(abs(Dm) > 1e-9*abs(m0.*m2.*m4));  % too few points for a curve
T(flat) = b0(flat)./max(m0(flat), 1);
T(~ok) = NaN;
end
