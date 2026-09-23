% MAKE_WINDOW_MOVIE  Animate the windows the Radon transform is fed.
%
%   matlab -batch "make_window_movie"
%
% Scans the solver's own windows in raster order - left to right along a
% row, then down to the next row - writing one GIF frame per window. Each
% frame shows where the window sits, the conditioned window with the fitted
% slope drawn on it with its semblance, and the criterion against candidate
% slope. After RUN_SLOPES, each window is marked against the calibrated gate.
%
% Use it to see why a window returned what it did: a flat criterion means
% there was nothing to lock onto, and a fitted line that does not lie along
% the layering means something else is pulling it.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'), fullfile(repo,'opr'), here);

args = {};
if exist('frame','var'), args = [args {'frame', frame}]; end
cfg = ls_config(args{:});

if ~exist('max_frames','var'), max_frames = 60; end
if ~exist('stride_x','var'),   stride_x = 6;    end
if ~exist('stride_z','var'),   stride_z = 2;    end

% Judge each window against the gate the solver calibrated, if it has run.
gate = {};
if exist(cfg.mat_file, 'file') == 2
    gate = {'semb_thresh', cfg.mat_file};
end

window_movie(cfg.data_file, cfg.gif_file, cfg.solver_args{:}, gate{:}, ...
    'max_frames', max_frames, 'stride_x', stride_x, 'stride_z', stride_z);
