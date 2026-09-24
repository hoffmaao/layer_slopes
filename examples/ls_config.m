function cfg = ls_config(varargin)
% LS_CONFIG  One place where the target frame and solver settings live.
%
%   cfg = LS_CONFIG()                  defaults
%   cfg = LS_CONFIG('frame','...')     override any field
%
% The three example scripts (RUN_SLOPES, MAKE_SLOPE_FIGURE,
% MAKE_WINDOW_MOVIE) all read this, so they cannot drift apart.
%
% THE WINDOW is 1000 m along track by 50 m in depth: square in the image
% the Radon transform sees, 200 x 200 samples at a 0.25 m grid and 20x
% vertical exaggeration. That is Holschuh et al.'s (2017) square window,
% carried into the exaggerated grid this radar needs. Measured against a
% reflector tracked by hand on 20250108_02_005 and 20250112_01_008:
%
%   window (m)   samples      solved   gain   r       adjacent-window change
%   1000 x 20    200 x 80     90-91%   1.03   0.99    0.031-0.034 deg
%   1000 x 50    200 x 200    94%      1.01   0.99    0.023-0.027 deg
%   400 x 400    80 x 1600    99-100%  0.71   0.88    0.017-0.021 deg
%   200 x 200    40 x 800     96%      0.75   0.67    0.034-0.043 deg
%
% The along-track length sets how precisely a dip is measured; the height
% stacks more layers onto one estimate but averages the dip over that much
% depth. 400 m squares in metres read the reflector 25-35% low for that
% reason, and 200 m is too short along track for a tenth of a degree.
% 1000 x 50 keeps the precision of the long window and adds enough height
% for a smooth, complete field.
%
% The other settings, and why:
%
%   detrend_len  A window spans part of the power-vs-depth decay, and that
%                gradient is horizontal. The 30 m high-pass - a local
%                quadratic fit, subtracted - removes it, including the
%                curvature where the firn's brightness falls away.
%   smooth_x     Trace-to-trace gain changes draw vertical stripes, which
%                the Radon will happily fit. 60 m removes them; much more
%                than that starts smearing the layers themselves.
%   exclude_z    Bands that are strong and horizontal but are not ice, left
%                out of every window that spans them (see max_excluded in
%                ROLLINGRADON_OPR).

cfg.season    = '2024_Antarctica_Ground2';
cfg.product   = 'CSARP_post/CSARP_standard_HH';   % SAR-focused: migrated
cfg.frame     = '20250108_02_005';
cfg.data_root = '/kucresis/scratch/dataproducts/opr_data/accum';
cfg.out_root  = '/kucresis/scratch/hoffmana_sta/layer_slopes/products';

% Depth range: from the surface to 1500 m. Two bands are not ice and are
% excluded: the surface ringdown of this ground-based system (0-28 m) and
% the merge of the first and second pulse returns near 75 m.
cfg.z_top     = 0;
cfg.z_bot     = 1500;
cfg.exclude_z = [0 28; 70 88];

% One colour scale for every figure, so frames and segments can be read
% against each other. Auto-scaling per frame ranged from +/-0.18 to
% +/-1.0 deg across one season's products, and the same colour meant a
% five-fold different dip from one figure to the next.
cfg.clim_dip  = [-0.3 0.3];      % deg

cfg.opts = { ...
    'grid_spacing',  0.25, ...   % m; the system resolves 0.53 m in ice
    'vert_exag',     20, ...     % exact, and 20x fewer pixels per window
    'window_x',      1000, ...   % m; 200 samples in the solver's image
    'window_z',      50, ...     % m; 200 samples
    'dip_max',       1, ...      % deg, true
    'dip_step',      0.005, ...  % deg, true
    'smooth_x',      60, ...     % m, along-track destriping
    'smooth_len',    1.5, ...    % m, depth low-pass
    'detrend_len',   30};        % m, depth high-pass

for i = 1:2:numel(varargin)
    v = varargin{i+1};
    if isstring(v) && isscalar(v)
        v = char(v);     % "..." and '...' both work; the paths below index chars
    end
    cfg.(char(varargin{i})) = v;
end

day_seg = cfg.frame(1:11);
cfg.data_file = fullfile(cfg.data_root, cfg.season, cfg.product, day_seg, ...
    ['Data_' cfg.frame '.mat']);
cfg.out_dir = fullfile(cfg.out_root, cfg.season);
cfg.mat_file = fullfile(cfg.out_dir, sprintf('slopes_%s.mat', cfg.frame));
cfg.fig_file = fullfile(cfg.out_dir, sprintf('slopes_%s.png', cfg.frame));
cfg.gif_file = fullfile(cfg.out_dir, sprintf('windows_%s.gif', cfg.frame));

% Everything the solver needs, assembled once.
cfg.solver_args = [cfg.opts, { ...
    'z_pad_surface', cfg.z_top, 'z_max', cfg.z_bot, ...
    'exclude_z', cfg.exclude_z}];
end
