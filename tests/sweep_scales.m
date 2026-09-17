% SWEEP_SCALES  Which window scales does the accumulation radar support?
%
% The along-track and vertical requirements are independent here:
%   * window_z sets the depth detail. The radar resolves 0.53 m and layers
%     sit ~8 m apart, so a short window still contains several cycles - this
%     is where the fine resolution of the image pays off.
%   * window_x sets the smallest dip that is measurable at all. A dip of
%     theta displaces a layer by window_x*tan(theta); that has to exceed a
%     useful fraction of one 0.53 m range cell.
% The floor printed for each row is the dip that gives exactly one range
% cell of displacement across the window.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));

data_file = ['/kucresis/scratch/dataproducts/opr_data/accum/' ...
    '2024_Antarctica_Ground2/CSARP_post/CSARP_standard/20250108_02/' ...
    'Data_20250108_02_005.mat'];

base = {'grid_spacing',0.25,'dip_max',1,'dip_accept',0.9,'dip_step',0.005, ...
    'z_pad_surface',30,'z_max',200,'exclude_z',[70 88], ...
    'smooth_x',250,'smooth_len',1.5,'detrend_len',0,'trace_balance',true, ...
    'solver_params',struct('o_f_horizontal',4,'o_f_vertical',4, ...
                           'snr_thresh',2,'snr_fac',1,'vr',1), ...
    'verbose',false};

wxs = [500 1000 2000 4000];
wzs = [10 15 20 30];

fprintf('%-7s %-7s %-8s %-7s %-6s %-8s %-7s %-6s\n', ...
    'wx(m)','wz(m)','floor','cells','pct','median','contin','signed');
for wz = wzs
    for wx = wxs
        ve = max(5, round(2000/wx)*20);     % keep the apparent dip sensible
        evalc(['R = RollingRadon_OPR(data_file, base{:}, ' ...
           '''window_x'',wx, ''window_z'',wz, ''vert_exag'',ve);']);
        v = R.slopes(isfinite(R.slopes));
        fl = atand(0.53/wx);
        if isempty(v)
            fprintf('%-7d %-7d %-8.3f %-7d %-6.1f   (none)\n', wx, wz, fl, 0, 0);
        else
            fprintf('%-7d %-7d %-8.3f %-7d %-6.1f %-8.3f %-7.3f %-6.2f\n', ...
                wx, wz, fl, numel(v), 100*numel(v)/numel(R.slopes), ...
                median(v), R.continuity, mean(sign(v)==sign(median(v))));
        end
    end
end
