% SWEEP_FAITHFUL  Are the additions here helping, or suppressing layers?
%
% Turns each addition off in turn on the real accumulation-radar frame and
% measures the effect. "signed" is the fraction of solved cells agreeing
% with the median sign, which is the honest coherence measure because
% smoothing cannot inflate it. "contin" is the median change between
% adjacent cells, which smoothing DOES inflate, so read the two together.
%
% Run against the 30-70 m band, where this radar genuinely has layering
% (2.05 dB band-passed contrast, discrete spectral peaks near 8 m).

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src')); addpath(fullfile(here,'..','opr'));

data_file = ['/kucresis/scratch/dataproducts/opr_data/accum/' ...
    '2024_Antarctica_Ground2/CSARP_post/CSARP_standard/20250108_02/' ...
    'Data_20250108_02_005.mat'];

common = {'grid_spacing',0.25,'window_x',2000,'window_z',30, ...
    'dip_max',1,'dip_accept',0.9,'dip_step',0.005, ...
    'z_pad_surface',30,'z_max',70,'verbose',false};

nick = struct('o_f_horizontal',2,'o_f_vertical',6,'snr_thresh',2, ...
              'snr_fac',1,'vr',3,'radon_method',0);
mine = struct('o_f_horizontal',4,'o_f_vertical',4,'snr_thresh',4, ...
              'vr',1,'radon_method',0);

cases = { ...
% label                 vexag smooth_x smooth_len tbal   params
 'Nick defaults       ',   1,    0,       0,     false,  nick
 'Nick + vert_exag    ',  20,    0,       0,     false,  nick
 've + trace balance  ',  20,    0,       0,     true,   nick
 've + tb + smooth_len',  20,    0,       1.5,   true,   nick
 've + tb + sl + sx250',  20,  250,       1.5,   true,   nick
 'all + my gating     ',  20,  250,       1.5,   true,   mine
 'all, sx = 60        ',  20,   60,       1.5,   true,   mine
 'all, no trace bal   ',  20,  250,       1.5,   false,  mine };

fprintf('%-22s %-7s %-6s %-8s %-7s %-6s\n', ...
    'configuration','cells','pct','median','contin','signed');
for i = 1:size(cases,1)
    evalc(['R = RollingRadon_OPR(data_file, common{:}, ' ...
      '''vert_exag'',cases{i,2}, ''smooth_x'',cases{i,3}, ' ...
      '''smooth_len'',cases{i,4}, ''detrend_len'',0, ' ...
      '''trace_balance'',cases{i,5}, ''solver_params'',cases{i,6});']);
    v = R.slopes(isfinite(R.slopes));
    if isempty(v)
        fprintf('%-22s %-7d %-6.1f   (none solved)\n', cases{i,1}, 0, 0);
        continue
    end
    fprintf('%-22s %-7d %-6.1f %-8.3f %-7.3f %-6.2f\n', ...
        cases{i,1}, numel(v), 100*numel(v)/numel(R.slopes), ...
        median(v), R.continuity, mean(sign(v)==sign(median(v))));
end
