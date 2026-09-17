function run_tests()
% RUN_TESTS  Run the layer-slope regression suite.
%
%   matlab -batch "run_tests"
%
% Exits non-zero if anything fails, so it can gate a deploy.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here,'..','src'));
addpath(fullfile(here,'..','opr'));

tests = {@test_sign, @test_regrid, @test_opr_units};
total = 0;

for i = 1:numel(tests)
    name = func2str(tests{i});
    fprintf('\n=== %s ===\n', name);
    try
        total = total + tests{i}();
    catch ME
        fprintf('  ERROR %s: %s\n', name, ME.message);
        for k = 1:numel(ME.stack)
            fprintf('      at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        total = total + 1;
    end
end

fprintf('\n========================================\n');
if total == 0
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('%d FAILURE(S)\n', total);
    exit(1);
end
end
