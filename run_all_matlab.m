% RUN_ALL_MATLAB Run all MATLAB estimation and plotting programs in order.
%
% Each estimation driver determines its paths from its own file location.
% graph_rep.m runs last and produces Figures 1, 3, and 4 after the main
% estimation program has created the required results MAT file.

clear;
clc;
diary off;

runner_path = mfilename('fullpath');
if isempty(runner_path)
    error('Run run_all_matlab.m as a script file.');
end
replication_dir = fileparts(runner_path);
cd(replication_dir);

fprintf('\nRunning preferred specification and counterfactuals...\n');
run(fullfile(replication_dir, 'main_UniformMSRP.m'));

% The analysis scripts begin with CLEAR ALL, so reconstruct runner state.
replication_dir = fileparts(mfilename('fullpath'));
cd(replication_dir);
fprintf('\nRunning alternative BLP-instrument specification...\n');
run(fullfile(replication_dir, 'demand_BLPIV.m'));

replication_dir = fileparts(mfilename('fullpath'));
cd(replication_dir);
fprintf('\nRunning alternative price-over-income specification...\n');
run(fullfile(replication_dir, 'demand_priceoverincome', ...
    'demand_priceoverincome.m'));

replication_dir = fileparts(mfilename('fullpath'));
cd(replication_dir);
fprintf('\nRunning graphing...\n');
run(fullfile(replication_dir, 'graph_rep.m'));

fprintf('\nAll MATLAB estimation and plotting programs completed.\n');
