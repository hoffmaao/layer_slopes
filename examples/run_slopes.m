% RUN_SLOPES  Compute the layer slope field and save it.
%
%   matlab -batch "run_slopes"
%   matlab -batch "frame='20250108_02_001'; run_slopes"
%
% Edit examples/ls_config.m to change the target or the settings.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
repo = fileparts(here);
addpath(fullfile(repo,'src'), fullfile(repo,'opr'), here);

args = {};
if exist('frame','var'),   args = [args {'frame', frame}];     end
if exist('product','var'), args = [args {'product', product}]; end
if exist('season','var'),  args = [args {'season', season}];   end
cfg = ls_config(args{:});

fprintf('=== layer slopes: %s ===\n', cfg.frame);
fprintf('  in  : %s\n', cfg.data_file);
fprintf('  out : %s\n', cfg.mat_file);

R = RollingRadon_OPR(cfg.data_file, cfg.solver_args{:}, ...
    'out_file', cfg.mat_file, 'verbose', true);

fprintf('=== complete ===\n');
