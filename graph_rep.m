% Replication plotting driver for Figures 1, 3, and 4.
%
% INPUTS
%   Figure 1: NEVproduction.xls.
%   Figure 3: subsidy&income.xls (sheet1).
%   Figure 4: dated results MAT file produced by main_UniformMSRP.m.
%
% OUTPUTS
%   evproduction.tex, incomesubsidygraph.tex, and incomeandsimulated.tex.
%
% EXECUTION NOTES
%   The script identifies the replication directory, adds the bundled
%   matlab2tikz source directory automatically, and uses the most recently
%   generated results_*.mat file for Figure 4. Open this file in UTF-8 mode.

%% Initialize plotting environment
SCRIPT_FULL_PATH = mfilename('fullpath');
if isempty(SCRIPT_FULL_PATH)
    REPLICATION_DIR = pwd;
else
    REPLICATION_DIR = fileparts(SCRIPT_FULL_PATH);
end
cd(REPLICATION_DIR);

MATLAB2TIKZ_SRC = fullfile(REPLICATION_DIR, ...
    'matlab2tikz-matlab2tikz-806c97d', 'src');
if ~isfolder(MATLAB2TIKZ_SRC)
    error('Bundled matlab2tikz source directory not found: %s', MATLAB2TIKZ_SRC);
end
addpath(MATLAB2TIKZ_SRC);

%% Figure 1: EV production portion
evprod = readtable('NEVproduction.xls','Format','auto') ;
evprod.Properties.VariableNames={'time','evprod', 'totalprod' };
evprod.t = datetime(evprod.time,'InputFormat','yyyy-MM');

evprod.evp = evprod.evprod ./ evprod.totalprod ;
evprod(isnan(evprod.evp)==1,:)=[];
evprod.evp = evprod.evp *100;

plot(evprod.t,evprod.evp,'--bo','LineWidth',2,'color','blue','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red' )
xlabel('Year')
ylabel('Electric Vehicle Production Percentage (%)')

matlab2tikz('evproduction.tex')

clear evprod

%% Figure 3: Subsidy Distribution Based on Income
clear

%average
data = xlsread('subsidy&income.xls','sheet1');
%sales_weighted average
% data = xlsread('subsidy&income.xls','sheet2');

category = data(:, 3);
num_subsets = max(unique(data(:,3)));
avg = zeros(1, num_subsets);
max_val = zeros(1, num_subsets);
min_val = zeros(1, num_subsets);
low = zeros(1, num_subsets);
high = zeros(1, num_subsets);
range = zeros(num_subsets, 2);

for subset_idx = 1:num_subsets
    subset_data = data(category == subset_idx, :);
    range(subset_idx, :) = ceil([min(subset_data(:, 1)), max(subset_data(:, 1))])*12/1000;
    subsidy = subset_data(:, 2)*10;
    avg(subset_idx) = mean(subsidy);
    max_val(subset_idx) = max(subsidy);
    min_val(subset_idx) = min(subsidy);
    low(subset_idx) = avg(subset_idx) - min_val(subset_idx);
    high(subset_idx) = max_val(subset_idx) - avg(subset_idx);
end

fig = figure(1);
set(fig, 'Position', [300, 200, 800, 500], 'Color', [1 1 1])
errorbar(1:num_subsets, avg, low, high, 'ok', 'MarkerSize', 4, 'CapSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2);  

box on
axis([0.5 num_subsets+0.5 15 25])

set(gca, 'xtick', 1:num_subsets, 'xticklabel', {'<= 160', '160-190', '190-230','230-300', '>=300'}, ...
    'ytick', 15:2.5:25,'yticklabel',15:2.5:25, 'FontSize', 14, 'Linewidth', 0.7, 'FontName', 'Times New Roman', 'XTickLabelRotation', 0)

% title('title', 'FontName', 'Times New Roman', 'fontsize', 16);
xlabel('Annual Income (RMB ,000)');
ylabel('Subsidy (RMB ,000)');

%average
matlab2tikz('incomesubsidygraph.tex')
%sales_weighted average
% matlab2tikz('incomesubsidygraph_salesweighted.tex')


%% Figure 4: buyer income versus simulated population income
% Figure 3 clears the workspace, so reconstruct the replication root.
SCRIPT_FULL_PATH = mfilename('fullpath');
REPLICATION_DIR = fileparts(SCRIPT_FULL_PATH);
result_candidates = dir(fullfile(REPLICATION_DIR, 'diary', 'results_*.mat'));
if isempty(result_candidates)
    error(['No Figure 4 result file was found in %s. Run ' ...
        'main_UniformMSRP.m first.'], fullfile(REPLICATION_DIR, 'diary'));
end
[~, newest_result_index] = max([result_candidates.datenum]);
result_file = fullfile(result_candidates(newest_result_index).folder, ...
    result_candidates(newest_result_index).name);
fprintf('Figure 4 input: %s\n', result_file);

result_data = load(result_file, 's_inc', 'demogr','data');


figure;
h1 = histogram(result_data.s_inc * 12 / 1000);
hold on;
h2 = histogram(result_data.demogr(:),40);

xlabel('Annual Income (RMB ''000)');
ylabel('Frequency');

hLegend = legend([h1, h2], ...
    'Vehicle Buyers'' Incomes', ...
    'Population Incomes', ...
    'Location', 'northeast');

hLegend.ItemTokenSize = [7 7];

matlab2tikz('incomeandsimulated.tex');
hold off
