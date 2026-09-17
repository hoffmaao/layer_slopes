% SWEEP_COVERAGE  How small must the window be to fill the domain?
%
% Edge loss is window_x/2 at each end and window_z/2 top and bottom, so a
% 2 km window cannot reach within 1 km of either end of the line. Shrinking
% the window costs dip sensitivity, because a 0.1 deg dip has to displace a
% layer by more than one 0.53 m range cell across it. vert_exag is raised to
% keep the apparent dip well inside the search range.
%
% "x span" is the fraction of the line covered by window centres; "z max" is
% the deepest solved cell; "signed" is the fraction of cells agreeing with
% the median sign.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));

data_file = ['/kucresis/scratch/dataproducts/opr_data/accum/' ...
    '2024_Antarctica_Ground2/CSARP_post/CSARP_standard/20250108_02/' ...
    'Data_20250108_02_005.mat'];

base = {'grid_spacing',0.25,'dip_max',1,'dip_accept',0.9,'dip_step',0.005, ...
    'z_pad_surface',28,'z_max',1500,'exclude_z',[70 88], ...
    'smooth_x',250,'smooth_len',1.5,'detrend_len',0,'verbose',false};

cases = { ...
  % label        window_x  window_z  vert_exag  of_h  of_v
  'wx2000 wz40 ',   2000,     40,       20,      4,    4
  'wx1000 wz40 ',   1000,     40,       40,      4,    4
  'wx600  wz40 ',    600,     40,       60,      4,    4
  'wx400  wz40 ',    400,     40,       80,      4,    4
  'wx600  wz20 ',    600,     20,       60,      4,    4
  'wx400  wz20 ',    400,     20,       80,      4,    4 };

fprintf('%-13s %-7s %-6s %-7s %-7s %-7s %-7s %s\n', ...
    'case','cells','pct','x span','z max','median','contin','signed');
for i = 1:size(cases,1)
    t = tic;
    evalc(['R = RollingRadon_OPR(data_file, base{:}, ' ...
      '''window_x'',cases{i,2}, ''window_z'',cases{i,3}, ' ...
      '''vert_exag'',cases{i,4}, ''solver_params'', struct(''vr'',1, ' ...
      '''snr_thresh'',4,''o_f_horizontal'',cases{i,5},''o_f_vertical'',cases{i,6}));']);
    v = R.slopes(isfinite(R.slopes));
    if isempty(v)
        fprintf('%-13s %-7d %-6.1f %s\n', cases{i,1}, 0, 0, '(none)'); continue
    end
    [~, zi] = find(isfinite(R.slopes).');
    zmax = max(R.slope_z(unique(zi)));
    xspan = (max(R.slope_x)-min(R.slope_x))/R.grid.x(end);
    fprintf('%-13s %-7d %-6.1f %-7.2f %-7.0f %-7.3f %-7.3f %-5.2f  [%.0f s]\n', ...
        cases{i,1}, numel(v), 100*numel(v)/numel(R.slopes), xspan, zmax, ...
        median(v), R.continuity, mean(sign(v)==sign(median(v))), toc(t));
end
