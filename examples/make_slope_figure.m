% MAKE_SLOPE_FIGURE  Slope field over the power image, from a saved result.
%
%   matlab -batch "make_slope_figure"
%
% Run RUN_SLOPES first. The background is the echogram power in dB, and the
% slope field is drawn over it as a continuous raster rather than as the raw
% window cells - overlapping windows plotted as cells make a staircase whose
% edges come from the window spacing, not from the ice.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'), fullfile(repo,'opr'), here);

args = {};
if exist('frame','var'), args = [args {'frame', frame}]; end
cfg = ls_config(args{:});

if exist(cfg.mat_file,'file') ~= 2
    error('make_slope_figure:noResult', ...
        'No result at %s. Run run_slopes first.', cfg.mat_file);
end

plot_slope_field(cfg.mat_file, cfg.fig_file, 'clim_dip', cfg.clim_dip);
