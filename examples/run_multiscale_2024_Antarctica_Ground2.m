% RUN_MULTISCALE_2024_ANTARCTICA_GROUND2
%
% Multi-scale layer slope field for a mega-dune profile, resolving small
% features and small changes in dip.
%
%   matlab -batch "run_multiscale_2024_Antarctica_Ground2"
%   matlab -batch "frame='20250108_02_001'; run_multiscale_2024_Antarctica_Ground2"
%
% WHY MULTI-SCALE. On this radar the two window dimensions have independent
% limits, measured on this frame in tests/sweep_scales.m:
%
%   window_z is cheap. The system resolves 0.53 m and layers sit ~8 m apart,
%   so a 10 m window still holds several layer cycles. Fine depth detail is
%   what the image resolution buys you.
%
%   window_x is expensive. A dip of theta displaces a layer by
%   window_x*tan(theta), and that has to exceed a useful part of one range
%   cell. At 500 m the floor is 0.061 deg, the size of the signal itself,
%   and the field falls apart (sign consistency 0.61). At 2000 m the floor
%   is 0.015 deg and it holds (0.82).
%
% So no single window does both. SLOPE_MULTISCALE runs several and keeps the
% finest one whose answer the next coarser scale corroborates. On this frame
% that takes coverage from 25% (finest scale alone, sign consistency 0.67)
% to 58% at sign consistency 0.90, with 72% of cells coming from the finest
% scale - the small features survive where they are real.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'));
addpath(fullfile(repo,'opr'));

if ~exist('season','var'),  season  = '2024_Antarctica_Ground2';   end
if ~exist('product','var'), product = 'CSARP_post/CSARP_standard'; end
if ~exist('frame','var'),   frame   = '20250108_02_005';           end
if ~exist('data_root','var')
    data_root = '/kucresis/scratch/dataproducts/opr_data/accum';
end
if ~exist('out_dir','var')
    out_dir = fullfile('/kucresis/scratch/hoffmana_sta/layer_slopes/products', season);
end
if ~exist('scales','var')
    scales = [4000 30; 2000 20; 1000 10];   % [window_x window_z], coarse first
end

day_seg = frame(1:11);
data_file = fullfile(data_root, season, product, day_seg, ['Data_' frame '.mat']);
out_file = fullfile(out_dir, sprintf('multiscale_%s.mat', frame));
fig_file = fullfile(out_dir, sprintf('multiscale_%s.png', frame));

fprintf('=== multi-scale layer slopes: %s %s ===\n', season, frame);

M = slope_multiscale(data_file, ...
    'scales',        scales, ...
    'grid_spacing',  0.25, ...   % m; the system resolves 0.53 m
    'vert_exag',     20, ...     % makes 0.1 deg present to the Radon as 2 deg
    'dip_max',       1, ...      % deg true
    'dip_accept',    0.9, ...
    'dip_step',      0.005, ...
    'z_pad_surface', 30, ...
    'z_max',         200, ...    % below this the frame is at the noise floor
    'exclude_z',     [70 88], ...% merged first/second pulse return
    'smooth_x',      250, ...    % along-track destriping
    'smooth_len',    1.5, ...
    'detrend_len',   0, ...
    'solver_params', struct('o_f_horizontal',4,'o_f_vertical',4, ...
                            'snr_thresh',2,'snr_fac',1,'vr',1), ...
    'out_file',      out_file, ...
    'verbose',       true);

plot_slope_field(M, fig_file);
fprintf('=== complete ===\n');
