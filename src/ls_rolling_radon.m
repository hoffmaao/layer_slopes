function S = ls_rolling_radon(G, p)
% LS_ROLLING_RADON  Roll a window over a gridded echogram and fit layer slope.
%
%   S = LS_ROLLING_RADON(G, p)
%
%   G   grid from OPR_FLATTEN_GRID: .x (m), .x_solver, .z (m), .img,
%       .dx, .dz, .vert_exag
%   p   parameters:
%       .window_x, .window_z   window size in metres
%       .overlap_x, .overlap_z fraction of a window between centres, in
%                              (0,1]; 0.25 steps by a quarter window
%       .slope_max             search +/- this, in TRUE degrees
%       .slope_step            search step, in TRUE degrees
%       .slope_accept          discard results beyond this (true degrees)
%       .q_thresh              minimum peak/median criterion ratio
%       .surface_z, .bed_z     [1 x nx] gates in depth (m)
%       .whole_window_in_ice   require the full window inside the gates
%       .verbose
%
% Returns
%   .slope_x [1 x m]  window centres along track (m)
%   .slope_z [1 x n]  window centres in depth (m)
%   .slopes  [n x m]  slope, degrees, POSITIVE = layer RISES with +x
%   .q       [n x m]  criterion peak/median for each window
%   .status  [n x m]  0 solved, 1 outside ice, 2 low quality,
%                     3 rejected by the slope gate
%
% The window is rectangular on purpose. On ice-penetrating radar the two
% dimensions do different jobs: the along-track extent sets the smallest
% slope that is measurable at all, because a slope displaces a layer by
% window_x*tan(slope) and that has to exceed a useful part of one range
% cell; the vertical extent only has to hold enough layer cycles to define
% an orientation, and wants to stay short enough that the power envelope is
% not the strongest feature in the window.
%
% See also LS_RADON_DIP, OPR_FLATTEN_GRID

if ~isfield(p,'overlap_x') || isempty(p.overlap_x), p.overlap_x = 0.25; end
if ~isfield(p,'overlap_z') || isempty(p.overlap_z), p.overlap_z = 0.25; end
if ~isfield(p,'q_thresh')  || isempty(p.q_thresh),  p.q_thresh  = 1.5;  end
if ~isfield(p,'slope_accept') || isempty(p.slope_accept)
    p.slope_accept = p.slope_max;
end
if ~isfield(p,'whole_window_in_ice'), p.whole_window_in_ice = true; end
if ~isfield(p,'verbose'), p.verbose = true; end

[nz, nx] = size(G.img);

wx = max(5, round(p.window_x/G.dx));
wz = max(5, round(p.window_z/G.dz));
if wx > nx || wz > nz
    error('ls_rolling_radon:windowTooLarge', ...
        'Window %d x %d samples does not fit the %d x %d grid.', wx, wz, nx, nz);
end

sx = max(1, round(wx*p.overlap_x));
sz = max(1, round(wz*p.overlap_z));

c0 = 1:sx:(nx-wx+1);
r0 = 1:sz:(nz-wz+1);

% The grid is presented to the Radon as isotropic; vert_exag multiplies the
% apparent slope by an exactly known factor, undone on the way out.
ve = G.vert_exag;
slope_max_app = min(89, atand(ve*tand(p.slope_max)));
slope_step_app = atand(ve*tand(p.slope_step));

slopes = nan(numel(r0), numel(c0));
qq = nan(size(slopes));
status = zeros(size(slopes));

xc = G.x(c0 + floor(wx/2));
zc = G.z(r0 + floor(wz/2));

if p.verbose
    fprintf('  rolling %d x %d windows (%d x %d samples, step %d x %d)\n', ...
        numel(c0), numel(r0), wx, wz, sx, sz);
    fprintf('    one range cell across the window = %.4f deg\n', ...
        atand(0.53/p.window_x));
end

t0 = tic;
for i = 1:numel(c0)
    ci = c0(i):(c0(i)+wx-1);
    for j = 1:numel(r0)
        ri = r0(j):(r0(j)+wz-1);

        % --- ice-column gate -------------------------------------------
        if isfield(p,'bed_z') && ~isempty(p.bed_z)
            sg = interp1(G.x, p.surface_z, xc(i), 'linear', 'extrap');
            bg = interp1(G.x, p.bed_z,     xc(i), 'linear', 'extrap');
            if p.whole_window_in_ice
                lo = G.z(ri(1)); hi = G.z(ri(end));
            else
                lo = zc(j); hi = zc(j);
            end
            if ~(lo > sg && hi < bg)
                status(j,i) = 1;
                continue
            end
        end

        win = G.img(ri, ci);
        [sl_app, q] = ls_radon_dip(win, G.dz, G.dz, ...
            slope_max_app, slope_step_app);

        if ~isfinite(sl_app) || ~isfinite(q) || q < p.q_thresh
            status(j,i) = 2;
            qq(j,i) = q;
            continue
        end

        sl = atand(tand(sl_app)/ve);
        if abs(sl) > p.slope_accept
            status(j,i) = 3;
            qq(j,i) = q;
            continue
        end

        slopes(j,i) = sl;
        qq(j,i) = q;
    end
    if p.verbose && mod(i, max(1,round(numel(c0)/8))) == 0
        fprintf('    %3.0f%%  %.1f min\n', 100*i/numel(c0), toc(t0)/60);
    end
end

S = struct();
S.slope_x = xc(:).';
S.slope_z = zc(:).';
S.slopes = slopes;
S.q = qq;
S.status = status;
S.window_samples = [wx wz];
S.step_samples = [sx sz];
S.runtime_s = toc(t0);
end
