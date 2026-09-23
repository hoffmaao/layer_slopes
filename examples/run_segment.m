% RUN_SEGMENT  Process every frame of one segment and save per-frame results.
%
%   matlab -batch "day_seg='20250111_01'; run_segment"
%   matlab -batch "day_seg='20250111_01'; z_top=200; z_bot=700; run_segment"
%
% A segment is split into frames for storage, but a layer does not stop at a
% frame boundary - so a feature of interest usually needs the whole line.
% This runs the solver over every frame in the segment with identical
% settings and writes one .mat and one .png each, plus a summary table.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'), fullfile(repo,'opr'), here);

if ~exist('day_seg','var'), day_seg = '20250111_01'; end
if ~exist('z_top','var'),   z_top = 200; end
if ~exist('z_bot','var'),   z_bot = 700; end
if ~exist('make_figs','var'), make_figs = true; end

cfg = ls_config('frame', [day_seg '_001']);
seg_dir = fileparts(cfg.data_file);
f = dir(fullfile(seg_dir, ['Data_' day_seg '_*.mat']));
f = f(~contains({f.name}, '_img_'));
if isempty(f)
    error('run_segment:noFrames','No frames found in %s', seg_dir);
end

fprintf('=== segment %s: %d frames, depth %g-%g m ===\n', ...
    day_seg, numel(f), z_top, z_bot);

summary = {};
for k = 1:numel(f)
    frame = regexprep(f(k).name, '^Data_|\.mat$', '');
    c = ls_config('frame', frame);
    fprintf('\n--- %s (%d/%d) ---\n', frame, k, numel(f));
    try
        R = RollingRadon_OPR(c.data_file, c.opts{:}, ...
            'z_pad_surface', z_top, 'z_max', z_bot, ...
            'exclude_z', c.exclude_z, ...
            'out_file', c.mat_file, 'verbose', true);
        if make_figs
            plot_slope_field(R, c.fig_file, 'clim_dip', c.clim_dip);
        end
        v = R.slopes(isfinite(R.slopes));
        summary(end+1,:) = {frame, R.grid.x(end)/1000, ...
            100*numel(v)/numel(R.slopes), median(v), ...
            mean(sign(v)==sign(median(v))), R.continuity}; %#ok<SAGROW>
    catch ME
        fprintf('  SKIPPED: %s\n', ME.message);
        summary(end+1,:) = {frame, NaN, 0, NaN, NaN, NaN}; %#ok<SAGROW>
    end
end

fprintf('\n=== summary: %s, %g-%g m ===\n', day_seg, z_top, z_bot);
fprintf('%-22s %-8s %-8s %-9s %-7s %s\n', ...
    'frame','km','solved','median','signed','continuity');
for k = 1:size(summary,1)
    fprintf('%-22s %-8.1f %-8.0f%% %-9.3f %-7.2f %.3f\n', summary{k,:});
end
