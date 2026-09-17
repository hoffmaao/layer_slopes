% RUN_SLOPES_2024_ANTARCTICA_GROUND2
%
% Englacial layer slopes for the mega-dune profiles of the
% 2024_Antarctica_Ground2 accumulation-radar season (~86.7 S).
%
% Defaults to 20250108_02 frame 005 - the 20 km profile shown in
% imb.picker, whose CSARP_post/CSARP_standard echogram is the same image
% the picks are made on. That is deliberate: the slope field should be
% measured from exactly the product a person looks at.
%
%   matlab -batch "run_slopes_2024_Antarctica_Ground2"
%   matlab -batch "frame='20250108_02_001'; run_slopes_2024_Antarctica_Ground2"
%
% ---------------------------------------------------------------------
% GEOMETRY. In the picker view the bright horizon runs from ~1.85 us at
% 0 km to ~1.45 us at 20 km: about 34 m of relief over 20,000 m, a dip
% near 0.1 deg. Two things follow, and both are settings below.
%
%  * A 0.1 deg dip displaces a layer by 1.7 m over a 1000 m window and by
%    0.2 m over a 120 m one - less than the 0.53 m range resolution. The
%    window has to be ~1 km along track for the dip to exist in the data
%    at all.
%  * radon_ndh searched on a hard-coded 0.1 deg grid, so the whole signal
%    fitted inside one increment. dip_step is now explicit, and vert_exag
%    samples along track at 20x the vertical spacing so a 0.1 deg dip
%    presents to the Radon as ~2 deg. The mapping tan(apparent) =
%    vert_exag*tan(true) is exact and is inverted on output.
%
% No bed pick exists at this latitude, so the depth window is set
% explicitly to the interval that carries coherent layering.
% ---------------------------------------------------------------------

%% ---- repo paths -------------------------------------------------------
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'));
addpath(fullfile(repo,'opr'));

%% ---- what to process --------------------------------------------------
if ~exist('season','var'),  season  = '2024_Antarctica_Ground2';      end
if ~exist('product','var'), product = 'CSARP_post/CSARP_standard';    end
if ~exist('frame','var'),   frame   = '20250108_02_005';              end
if ~exist('data_root','var')
    data_root = '/kucresis/scratch/dataproducts/opr_data/accum';
end
if ~exist('out_dir','var')
    out_dir = fullfile('/kucresis/scratch/hoffmana_sta/layer_slopes/products', season);
end

day_seg = frame(1:11);
data_file = fullfile(data_root, season, product, day_seg, ['Data_' frame '.mat']);
pname = strrep(product,'/','_');
out_file = fullfile(out_dir, sprintf('slopes_%s_%s.mat', pname, frame));
fig_file = fullfile(out_dir, sprintf('slopes_%s_%s.png', pname, frame));

%% ---- depth window -----------------------------------------------------
% 60-300 m spans the horizons visible in the picker view (0.85-3.2 us).
% The first and second pulse returns merge near 75 m on this system, and
% that band is a strong horizontal feature with no stratigraphic meaning.
% Rather than starting below it, blank it and solve on both sides: the
% shallow layering above it is some of the clearest in the profile.
if ~exist('z_top','var'),     z_top = 28;        end
if ~exist('z_bot','var'),     z_bot = 300;       end
if ~exist('exclude_z','var'), exclude_z = [70 88]; end
if ~exist('window_z','var'),  window_z = 40;     end
if ~exist('window_x','var'),  window_x = 2000;   end

%% ---- solver settings --------------------------------------------------
opts = { ...
    'grid_spacing',  0.25, ...  % m vertical; the system resolves 0.53 m
    'vert_exag',     20, ...    % along track at 5 m, near the 5.8 m native
                        ...     % spacing; makes 0.1 deg present as ~2 deg
    'window_x',      window_x, ...% m. A 0.1 deg dip throws a layer 3.5 m over
                        ...     % this, against a 0.53 m range cell. Long
                        ...     % windows are what make the dip measurable;
                        ...     % the overlap below buys the detail back.
    'window_z',      window_z, ...% m. Small enough to fit between the
                        ...     % surface and the excluded pulse-merge band,
                        ...     % so the shallow section is solved too.
    'exclude_z',     exclude_z, ...
    'dip_max',       1, ...     % deg true
    'dip_accept',    0.9, ...   % deg true
    'dip_step',      0.005, ... % deg true; was effectively 0.1 before
    'z_pad_surface', z_top, ...
    'z_max',         z_bot, ...
    'smooth_x',      250, ...   % m along-track low-pass. Trace-to-trace gain
                        ...     % changes show up as vertical stripes, and a
                        ...     % stripe is a strong linear feature the Radon
                        ...     % will fit. Removing them is nearly free: a
                        ...     % 0.1 deg layer moves 0.4 m over 250 m, under
                        ...     % one range cell. This is the single biggest
                        ...     % lever on how coherent the field looks.
    'smooth_len',    1.5, ...   % m, depth low-pass at the layer scale
    'detrend_len',   0, ...     % m, no depth high-pass
    'solver_params', struct( ...
        'vr', 1, ...            % keep each measurement, do not substitute
        'snr_thresh', 4, ...
        'o_f_horizontal', 16, ...  % step the long window finely so the field
        'o_f_vertical', 12), ...   % is sampled at ~125 m, not ~1 km
    'out_file',      out_file, ...
    'verbose',       true};

fprintf('=== layer slopes: %s %s %s ===\n', season, product, frame);
fprintf('  in  : %s\n', data_file);
fprintf('  out : %s\n', out_file);

R = RollingRadon_OPR(data_file, opts{:});

%% ---- standard figure --------------------------------------------------
plot_slope_field(R, fig_file);

fprintf('=== complete ===\n');
