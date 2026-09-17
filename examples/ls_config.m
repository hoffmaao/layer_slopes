function cfg = ls_config(varargin)
% LS_CONFIG  One place where the target frame and solver settings live.
%
%   cfg = LS_CONFIG()                  defaults
%   cfg = LS_CONFIG('frame','...')     override any field
%
% The three example scripts (RUN_SLOPES, MAKE_SLOPE_FIGURE,
% MAKE_WINDOW_MOVIE) all read this, so they cannot drift apart.
%
% SETTINGS ARE CALIBRATED, not guessed. tests/validate_horizon.m tracks the
% bright reflector trace by trace with no Radon involved, differentiates it,
% and compares. Against that truth these give a regression gain of 0.96 and
% r = 0.98 - the magnitude is right, not just the sign. What set them:
%
%   detrend_len  A window spans part of the power-vs-depth decay, and that
%                gradient is horizontal. Without the high-pass it dominates
%                the window: gain 0.70, r 0.78. With it: gain 0.86, r 0.97.
%   window_x     A long window averages over a range of true dips and
%                regresses toward their mean. Gain runs 0.99 at 750 m, 0.96
%                at 1000 m, 0.86 at 2000 m, 0.71 at 3000 m, while the
%                smallest measurable dip goes the other way: 0.040, 0.030,
%                0.015, 0.010 deg. 1000 m balances the two.
%   window_z     Sampling MORE layers does not help. Gain falls from 0.86
%                at 20 m to 0.74 at 100 m and coverage collapses, because
%                dip varies with depth and a tall window averages across it.
%   smooth_x     Trace-to-trace gain changes draw vertical stripes, which
%                the Radon will happily fit. 60 m removes them; much more
%                than that starts smearing the layers themselves.

cfg.season    = '2024_Antarctica_Ground2';
cfg.product   = 'CSARP_post/CSARP_standard_HH';   % SAR-focused: migrated
cfg.frame     = '20250108_02_005';
cfg.data_root = '/kucresis/scratch/dataproducts/opr_data/accum';
cfg.out_root  = '/kucresis/scratch/hoffmana_sta/layer_slopes/products';

% Depth window. The first and second pulse returns merge near 75 m on this
% system; that band is blanked so windows overlapping it abstain, while the
% layering above and below is still solved.
cfg.z_top     = 28;
cfg.z_bot     = 200;
cfg.exclude_z = [70 88];

cfg.opts = { ...
    'grid_spacing',  0.25, ...   % m; the system resolves 0.53 m in ice
    'vert_exag',     20, ...     % exact, and 20x fewer pixels per window
    'window_x',      1000, ...   % m
    'window_z',      20, ...     % m
    'dip_max',       1, ...      % deg, true
    'dip_step',      0.005, ...  % deg, true
    'smooth_x',      60, ...     % m, along-track destriping
    'smooth_len',    1.5, ...    % m, depth low-pass
    'detrend_len',   30};        % m, depth high-pass

for i = 1:2:numel(varargin)
    cfg.(varargin{i}) = varargin{i+1};
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
