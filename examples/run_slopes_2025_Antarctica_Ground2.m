% RUN_SLOPES_2025_ANTARCTICA_GROUND2
%
% Englacial layer slopes for a Ridge A accumulation-radar HH profile from
% the 2025_Antarctica_Ground2 season.
%
% Run headless on the CReSIS machines with
%
%   cd <repo>/examples
%   matlab -batch "run_slopes_2025_Antarctica_Ground2"
%
% Override the frame or any solver setting by defining the variable before
% calling, e.g.
%
%   matlab -batch "frame='20260109_02_001'; run_slopes_2025_Antarctica_Ground2"
%
% Products go to out_dir; the echogram itself is opened read-only and is
% never modified.

%% ---- repo paths -------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'));
addpath(fullfile(repo,'opr'));

%% ---- what to process --------------------------------------------------
if ~exist('season','var'),   season   = '2025_Antarctica_Ground2'; end
if ~exist('product','var'),  product  = 'CSARP_standard_HH';       end
if ~exist('frame','var'),    frame    = '20260109_02_003';         end
if ~exist('data_root','var')
    data_root = '/kucresis/scratch/dataproducts/opr_data/accum';
end
if ~exist('out_dir','var')
    out_dir = fullfile('/kucresis/scratch/hoffmana_sta/layer_slopes/products', season);
end

day_seg = frame(1:11);            % 20260109_02
data_file = fullfile(data_root, season, product, day_seg, ...
    ['Data_' frame '.mat']);
out_file = fullfile(out_dir, sprintf('slopes_%s_%s.mat', product, frame));

%% ---- solver settings --------------------------------------------------
% grid_spacing is the main cost/resolution knob. The Radon transform needs
% square pixels, and the native grid is 0.14 m vertically against 5.9 m
% along track, so the echogram is regridded isotropically once up front
% rather than 42x per window inside radon_ndh.
opts = { ...
    'grid_spacing',  0.25, ...  % m. The system resolves 0.53 m in ice, so
                        ...     % the grid has to sample that: at 1-2 m the
                        ...     % anti-alias average wipes the layering out
                        ...     % before the Radon transform ever sees it.
    'window_x',      120, ...   % m along track (~20 traces at 5.9 m spacing)
    'window_z',      30, ...    % m vertical (~57 range cells). Short enough
                        ...     % that system power drift does not become the
                        ...     % strongest linear feature in the window.
    'dip_max',       12, ...    % deg. Internal layers here dip a few degrees.
    'z_pad_surface', 30, ...    % m, skip the firn / feedthrough
    'z_pad_bed',     25, ...    % m, stop above the bed
    'smooth_len',    0, ...     % m, no extra low-pass: it would remove the
                        ...     % 0.53 m layering we are trying to measure.
    'detrend_len',   10, ...    % m, depth high-pass. Only ~10%% of the
                        ...     % variance sits below 4 m wavelength, so the
                        ...     % long envelope has to go or it dominates.
    'out_file',      out_file, ...
    'verbose',       true};

fprintf('=== layer slopes: %s %s %s ===\n', season, product, frame);
fprintf('  in  : %s\n', data_file);
fprintf('  out : %s\n', out_file);

R = RollingRadon_OPR(data_file, opts{:});

%% ---- standard figure --------------------------------------------------
% Slope field over the echogram power image, 10*log10(Data).
fig_file = fullfile(out_dir, sprintf('slopes_%s_%s.png', product, frame));
plot_slope_field(R, fig_file);

fprintf('=== complete ===\n');
