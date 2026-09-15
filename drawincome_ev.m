function [Income_Mean, EV_Adoption_Mean,EV_Adoption_Std,...
    GroupCount] = drawincome_ev(income_sorted, s_iev_sorted, n_bins)
% DRAWINCOME_EV Summarize simulated EV adoption across income quantile bins.
%
% INPUTS
%   income_sorted - simulated income values, in the income units used
%                   by the calling section of the main script.
%   s_iev_sorted  - aligned individual EV-adoption probabilities.
%   n_bins        - positive integer number of requested quantile groups.
%
% OUTPUTS
%   Income_Mean      - mean income in each nonempty bin.
%   EV_Adoption_Mean - mean EV-adoption probability in each bin.
%   EV_Adoption_Std  - within-bin adoption standard deviation.
%   GroupCount       - counts returned by GRPSTATS for the grouped variables.

% Using quantile for boxing data
% n_bins = 10;
edges = quantile(income_sorted, linspace(0, 1, n_bins+1));
bin_idx = discretize(income_sorted, edges);

% Getting necessary statistics
[mean_vals, std_vals, min_vals, max_vals, count_vals] = ...
    grpstats([income_sorted, s_iev_sorted], bin_idx, ...
    {'mean', 'std', 'min', 'max', 'numel'});

% Creating bin_stats table
Income_Mean = mean_vals(:, 1);
EV_Adoption_Mean = mean_vals(:, 2);
EV_Adoption_Std = std_vals(:, 2);
GroupCount = count_vals;

bin_stats = table();
bin_stats.Bin = (1:size(mean_vals, 1))';
bin_stats.Income_Mean = Income_Mean;
bin_stats.EV_Adoption_Mean = EV_Adoption_Mean;
bin_stats.Income_Std = std_vals(:, 1);
bin_stats.EV_Adoption_Std = EV_Adoption_Std;
bin_stats.Income_Min = min_vals(:, 1);
bin_stats.EV_Adoption_Min = min_vals(:, 2);
bin_stats.Income_Max = max_vals(:, 1);
bin_stats.EV_Adoption_Max = max_vals(:, 2);
bin_stats.GroupCount = GroupCount;

disp(bin_stats);

end
