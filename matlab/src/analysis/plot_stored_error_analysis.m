%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stored CMM Error Analysis Results
% Loads and plots previously computed error analysis data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

% Load stored results
if ~exist('cmm_error_analysis_fresh.mat', 'file')
    error('cmm_error_analysis_fresh.mat not found. Please run cmm_error_analysis.m first to generate the data.');
end

fprintf('Loading stored error analysis results...\n');
load('cmm_error_analysis_fresh.mat');

% Extract data from results structure
N_remap_test = results_struct.N_remap_test;
N_remap_benchmark = results_struct.N_remap_benchmark;
target_times = results_struct.target_times;
errors = results_struct.errors;

fprintf('Loaded results:\n');
fprintf('  N_remap test values: %s\n', mat2str(N_remap_test));
fprintf('  N_remap benchmark: %d\n', N_remap_benchmark);
fprintf('  Target times: %s\n', mat2str(target_times));
fprintf('  Number of methods: %d\n', length(N_remap_test));
fprintf('  Number of time points: %d\n', length(target_times));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results
figure('Position', [100, 100, 1400, 900]);

% Colors, markers, and line styles for different times
n_times = length(target_times);
if n_times == 2
    colors = {'b', 'r'};
    markers = {'o', 's'};
    line_styles = {'-', '--'};
    time_labels = {sprintf('t = %.0f', target_times(1)), sprintf('t = %.0f', target_times(2))};
elseif n_times == 3
    colors = {'b', 'r', 'g'};
    markers = {'o', 's', '^'};
    line_styles = {'-', '--', ':'};
    time_labels = {sprintf('t = %.0f', target_times(1)), sprintf('t = %.0f', target_times(2)), sprintf('t = %.0f', target_times(3))};
else
    colors = lines(n_times);
    markers = {'o', 's', '^', 'd', 'v', '<', '>', 'p', 'h'};
    line_styles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-'};
    time_labels = cell(n_times, 1);
    for i = 1:n_times
        time_labels{i} = sprintf('t = %.0f', target_times(i));
    end
end

% Plot 1: Electric Field Error
subplot(2, 2, 1);
hold on;
for t_idx = 1:n_times
    if iscell(colors)
        color = colors{t_idx};
        marker = markers{t_idx};
        line_style = line_styles{t_idx};
    else
        color = colors(t_idx,:);
        marker = markers{min(t_idx, length(markers))};
        line_style = line_styles{min(t_idx, length(line_styles))};
    end
    loglog(N_remap_test, errors.electric_field(:, t_idx), ...
           [marker line_style], 'Color', color, 'LineWidth', 2.5, 'MarkerSize', 10, ...
           'DisplayName', time_labels{t_idx});
end
xlabel('N_{remap}');
ylabel('Electric Field L2 Relative Error');
title('Electric Field Error vs Remapping Frequency');
legend('Location', 'best');
grid on;

% Plot 2: Distribution Function Error  
subplot(2, 2, 2);
hold on;
for t_idx = 1:n_times
    if iscell(colors)
        color = colors{t_idx};
        marker = markers{t_idx};
        line_style = line_styles{t_idx};
    else
        color = colors(t_idx,:);
        marker = markers{min(t_idx, length(markers))};
        line_style = line_styles{min(t_idx, length(line_styles))};
    end
    loglog(N_remap_test, errors.distribution_function(:, t_idx), ...
           [marker line_style], 'Color', color, 'LineWidth', 2.5, 'MarkerSize', 10, ...
           'DisplayName', time_labels{t_idx});
end
xlabel('N_{remap}');
ylabel('Distribution Function L2 Relative Error');
title('Distribution Function Error vs Remapping Frequency');
legend('Location', 'best');
grid on;

% Plot 3: Potential Energy Error
subplot(2, 2, 3);
hold on;
for t_idx = 1:n_times
    if iscell(colors)
        color = colors{t_idx};
        marker = markers{t_idx};
        line_style = line_styles{t_idx};
    else
        color = colors(t_idx,:);
        marker = markers{min(t_idx, length(markers))};
        line_style = line_styles{min(t_idx, length(line_styles))};
    end
    loglog(N_remap_test, errors.potential_energy(:, t_idx), ...
           [marker line_style], 'Color', color, 'LineWidth', 2.5, 'MarkerSize', 10, ...
           'DisplayName', time_labels{t_idx});
end
xlabel('N_{remap}');
ylabel('Potential Energy Relative Error');
title('Potential Energy Error vs Remapping Frequency');
legend('Location', 'best');
grid on;

% Plot 4: Combined error comparison at last time point
subplot(2, 2, 4);
hold on;
t_idx = n_times; % Last time point
loglog(N_remap_test, errors.electric_field(:, t_idx), 'b o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Electric Field');
loglog(N_remap_test, errors.distribution_function(:, t_idx), 'r s-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Distribution Function');
loglog(N_remap_test, errors.potential_energy(:, t_idx), 'g ^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Potential Energy');
xlabel('N_{remap}');
ylabel('Relative Error');
title(sprintf('All Errors at %s', time_labels{t_idx}));
legend('Location', 'best');
grid on;

sgtitle('CMM Error Analysis: Stored Results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print summary
fprintf('\n=== STORED ERROR ANALYSIS SUMMARY ===\n');
fprintf('Benchmark: N_remap = %d (%.2f seconds)\n', N_remap_benchmark, results_struct.benchmark_time);
fprintf('\nError Summary:\n');

for t_idx = 1:n_times
    fprintf('\nErrors at t = %.1f:\n', target_times(t_idx));
    fprintf('%-8s | %12s | %15s | %15s\n', 'N_remap', 'E-field', 'Dist. Func.', 'Pot. Energy');
    fprintf('%-8s | %12s | %15s | %15s\n', repmat('-', 1, 8), repmat('-', 1, 12), repmat('-', 1, 15), repmat('-', 1, 15));
    
    for method_idx = 1:length(N_remap_test)
        fprintf('%-8d | %10.2e | %13.2e | %13.2e\n', ...
                N_remap_test(method_idx), ...
                errors.electric_field(method_idx, t_idx), ...
                errors.distribution_function(method_idx, t_idx), ...
                errors.potential_energy(method_idx, t_idx));
    end
end

% Print simulation times if available
if isfield(results_struct, 'test_results') && ~isempty(results_struct.test_results)
    fprintf('\nSimulation Times:\n');
    fprintf('%-8s | %12s\n', 'N_remap', 'Time (s)');
    fprintf('%-8s | %12s\n', repmat('-', 1, 8), repmat('-', 1, 12));
    fprintf('%-8d | %10.2f\n', N_remap_benchmark, results_struct.benchmark_time);
    for method_idx = 1:length(N_remap_test)
        fprintf('%-8d | %10.2f\n', N_remap_test(method_idx), results_struct.test_results{method_idx}.simulation_time);
    end
end

fprintf('\nPlots generated from stored data in: cmm_error_analysis_fresh.mat\n');