% run every test in this folder
%
%   >> run tests/runAllTests
%
% or from a shell:
%   matlab -batch "run('tests/runAllTests.m')"

testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
results = runtests(testDir);
disp(table(results))
assert(all([results.Passed]),'%d of %d tests failed',sum(~[results.Passed]),numel(results))
