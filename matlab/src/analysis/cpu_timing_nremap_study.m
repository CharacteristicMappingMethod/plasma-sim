%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPU Timing Analysis for CMM with Different N_remap Values
% 
% This script analyzes CPU time per iteration vs simulation time for
% the CMM method with varying N_remap parameters on Landau damping case.
%
% N_remap values tested: [5, 10, 20, 30, 40, 1000]
% Test case: Landau Damping
% Time limit: 40.0
% Method: CMM 
%
% Output:
%   - Plot of CPU time per iteration vs simulation time
%   - Saved data file with all timing results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%% Configuration
fprintf('=== CPU Timing Analysis for CMM N_remap Study ===\n');
fprintf('Test case: Landau Damping\n');
fprintf('Time limit: 10.0\n');
fprintf('Method: CMM\n\n');

% Define N_remap values to test
N_remap_values = [5,10, 20, 30,40,1000];
target_time = 40.0;
num_tests = length(N_remap_values);

% Define data filename and initialize flag
data_filename = 'cpu_timing_nremap_landau_damping_Tend40.mat';
run_simulations = true;  % Default to running simulations

% Check for existing data
if exist(data_filename, 'file')
    fprintf('Found existing CPU timing data. Loading...\n');
    load(data_filename, 'results');
    fprintf('Loaded data for N_remap values: [%s]\n', join(string(results.N_remap_values), ', '));
    fprintf('Skipping simulations - using existing data.\n\n');
    
    % Verify loaded data has required fields
    if ~isfield(results, 'cpu_times') || ~isfield(results, 'time_arrays')
        fprintf('Warning: Existing data is incomplete. Re-running simulations...\n');
        run_simulations = true;
    else
        run_simulations = false;
    end
else
    fprintf('No existing data found. Will run simulations and save results.\n\n');
    run_simulations = true;
end

if run_simulations
    % Initialize storage for results
    results = struct();
    results.N_remap_values = N_remap_values;
    results.target_time = target_time;
    results.cpu_times = cell(num_tests, 1);
    results.time_arrays = cell(num_tests, 1);
    results.total_cpu_time = zeros(num_tests, 1);
    results.mean_cpu_per_iter = zeros(num_tests, 1);
    results.final_time = zeros(num_tests, 1);

%% Main Simulation Loop
fprintf('Starting simulations...\n');

for i = 1:num_tests
    current_N_remap = N_remap_values(i);
    fprintf('Running simulation %d/%d (N_remap = %d)...\n', i, num_tests, current_N_remap);
    
    % Save current state and clear to avoid parameter conflicts (following measure_remap_error.m pattern)
    save('temp_cpu_timing_config.mat', 'N_remap_values', 'target_time', 'num_tests', 'results', 'i', 'data_filename');
    
    % Clear workspace and reload basics
    clear all
    addpath(genpath('./src/'),genpath('./params/'))
    DEFAULTS
    
    % Reload configuration
    load('temp_cpu_timing_config.mat');
    current_N_remap = N_remap_values(i);
    
    % Load fresh Landau damping parameters
    PARAMS_landau_damping;
    
    % Override specific parameters for this test
    params.method = "CMM";
    params.N_remap = current_N_remap;
    params.Tend = target_time;
    
    % Display current simulation parameters
    fprintf('  Method: %s\n', params.method);
    fprintf('  N_remap: %d\n', params.N_remap);
    fprintf('  Target time: %.1f\n', params.Tend);
    fprintf('  Grid size: %dx%d\n', params.Nx, params.Nv);
    fprintf('  Time step: %.3f\n', params.dt);
    
    % Run simulation and measure total time
    tic_total = tic();
    [params_result, data_result] = Sim(params);
    total_simulation_time = toc(tic_total);
    
    % Extract timing data
    cpu_times = params_result.tcpu;
    time_array = params_result.time_array;
    
    % Filter data to only include time <= target_time
    valid_indices = time_array <= target_time;
    cpu_times_filtered = cpu_times(valid_indices);
    time_array_filtered = time_array(valid_indices);
    
    % Store results
    results.cpu_times{i} = cpu_times_filtered;
    results.time_arrays{i} = time_array_filtered;
    results.total_cpu_time(i) = total_simulation_time;
    results.mean_cpu_per_iter(i) = mean(cpu_times_filtered);
    results.final_time(i) = time_array_filtered(end);
    
    fprintf('  Completed in %.2f seconds\n', total_simulation_time);
    fprintf('  Mean CPU time per iteration: %.4f seconds\n', results.mean_cpu_per_iter(i));
    fprintf('  Final simulation time: %.2f\n', results.final_time(i));
    fprintf('  Total iterations: %d\n\n', length(cpu_times_filtered));
    
    % Update progress display
    fprintf('Progress: %d/%d simulations completed\n', i, num_tests);
    
    % Save partial results in case of interruption
    save('temp_cpu_timing_results.mat', 'results', 'N_remap_values', 'target_time');
end

% Clean up temporary files
if exist('temp_cpu_timing_config.mat', 'file')
    delete('temp_cpu_timing_config.mat');
end
if exist('temp_cpu_timing_results.mat', 'file')
    delete('temp_cpu_timing_results.mat');
end

fprintf('All simulations completed!\n');

    % Save results data for future use
    save(data_filename, 'results');
    fprintf('CPU timing data saved as: %s\n', data_filename);
    fprintf('Next time you run this script, it will load existing data instead of re-running simulations.\n\n');
    
end  % End of if run_simulations

% Ensure variables are available for both simulation and loaded data paths
N_remap_values = results.N_remap_values;
target_time = results.target_time;
num_tests = length(N_remap_values);
data_filename = 'cpu_timing_nremap_landau_damping_Tend40.mat';  % Re-define for consistency

%% Data Processing and Analysis
fprintf('=== Analysis Summary ===\n');
fprintf('N_remap    Mean CPU/iter (s)    Total CPU (s)    Final Time\n');
fprintf('-------    -----------------    -------------    ----------\n');
for i = 1:num_tests
    fprintf('%6d    %17.6f    %13.2f    %10.2f\n', ...
        N_remap_values(i), results.mean_cpu_per_iter(i), ...
        results.total_cpu_time(i), results.final_time(i));
end
fprintf('\n');

%% Plotting and Visualization
fprintf('Creating plots...\n');

% Create main figure
figure('Position', [100, 100, 1200, 800]);

% Define colors for different N_remap values
colors = lines(num_tests);

% Subplot 1: CPU time per iteration vs simulation time
subplot(2, 2, 1);
hold on;
legend_entries = {};
for i = 1:num_tests
    plot(results.time_arrays{i}, results.cpu_times{i}, ...
         'Color', colors(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('N_{remap} = %d', N_remap_values(i)));
    legend_entries{end+1} = sprintf('N_{remap} = %d', N_remap_values(i));
end
hold off;
xlabel('Simulation Time');
ylabel('CPU Time per Iteration (s)');
title('CPU Time per Iteration vs Simulation Time');
legend('Location', 'best', 'FontSize', 10);
grid on;
xlim([0, target_time]);

% Subplot 2: Mean CPU time per iteration for each N_remap
subplot(2, 2, 2);
bar(1:num_tests, results.mean_cpu_per_iter, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('N_{remap} Value');
ylabel('Mean CPU Time per Iteration (s)');
title('Mean CPU Time per Iteration vs N_{remap}');
set(gca, 'XTick', 1:num_tests, 'XTickLabel', N_remap_values);
grid on;
for i = 1:num_tests
    text(i, results.mean_cpu_per_iter(i) + max(results.mean_cpu_per_iter)*0.02, ...
         sprintf('%.4f', results.mean_cpu_per_iter(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

% Subplot 3: Total CPU time for each simulation
subplot(2, 2, 3);
bar(1:num_tests, results.total_cpu_time, 'FaceColor', [0.8, 0.4, 0.2]);
xlabel('N_{remap} Value');
ylabel('Total CPU Time (s)');
title('Total CPU Time vs N_{remap}');
set(gca, 'XTick', 1:num_tests, 'XTickLabel', N_remap_values);
grid on;
for i = 1:num_tests
    text(i, results.total_cpu_time(i) + max(results.total_cpu_time)*0.02, ...
         sprintf('%.1f', results.total_cpu_time(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

% Subplot 4: CPU efficiency (iterations per second)
subplot(2, 2, 4);
iterations_per_second = 1 ./ results.mean_cpu_per_iter;
bar(1:num_tests, iterations_per_second, 'FaceColor', [0.2, 0.8, 0.4]);
xlabel('N_{remap} Value');
ylabel('Iterations per Second');
title('Computational Efficiency vs N_{remap}');
set(gca, 'XTick', 1:num_tests, 'XTickLabel', N_remap_values);
grid on;
for i = 1:num_tests
    text(i, iterations_per_second(i) + max(iterations_per_second)*0.02, ...
         sprintf('%.1f', iterations_per_second(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

% Add main title
sgtitle('CMM N_{remap} CPU Timing Analysis - Landau Damping Case', 'FontSize', 16, 'FontWeight', 'bold');

% Optional: Save the plot (create directory if needed)
plot_filename = 'cpu_timing_nremap_landau_damping';
try
    % Try to save to images directory
    if ~exist('./images', 'dir')
        mkdir('./images');
        fprintf('Created ./images directory\n');
    end
    saveas(gcf, ['./images/' plot_filename '.png']);
    saveas(gcf, ['./images/' plot_filename '.fig']);
    fprintf('Plot saved as: ./images/%s.png and ./images/%s.fig\n', plot_filename, plot_filename);
catch
    % If that fails, save to current directory
    saveas(gcf, [plot_filename '.png']);
    saveas(gcf, [plot_filename '.fig']);
    fprintf('Plot saved to current directory as: %s.png and %s.fig\n', plot_filename, plot_filename);
end

%% Final Summary
fprintf('\n=== Final Summary ===\n');
fprintf('N_remap values analyzed: [%s]\n', join(string(N_remap_values), ', '));
fprintf('Target simulation time: %.1f\n', target_time);
fprintf('Number of configurations: %d\n', num_tests);

% Find optimal N_remap (balance between speed and accuracy)
[~, fastest_idx] = min(results.mean_cpu_per_iter);
[~, slowest_idx] = max(results.mean_cpu_per_iter);

fprintf('\nPerformance Analysis:\n');
fprintf('Fastest configuration: N_remap = %d (%.6f s/iter)\n', ...
    N_remap_values(fastest_idx), results.mean_cpu_per_iter(fastest_idx));
fprintf('Slowest configuration: N_remap = %d (%.6f s/iter)\n', ...
    N_remap_values(slowest_idx), results.mean_cpu_per_iter(slowest_idx));
fprintf('Speed ratio (slowest/fastest): %.2fx\n', ...
    results.mean_cpu_per_iter(slowest_idx) / results.mean_cpu_per_iter(fastest_idx));

fprintf('\nPlots and data have been saved to the analysis directory.\n');
fprintf('Analysis complete!\n');