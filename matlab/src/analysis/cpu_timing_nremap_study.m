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
fprintf('Test case: Two Stream Instability\n');
fprintf('Time limit: 40.0\n');
fprintf('Method: CMM\n\n');

% Define N_remap values to test
N_remap_values = [5,10, 20, 30,40,1000];
target_time = 40.0;
num_tests = length(N_remap_values);

% Define data filename and initialize flag
data_filename = 'cpu_timing_nremap_two_stream_Tend40.mat';
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
    
    % Load fresh two stream instability parameters
    PARAMS_two_stream;
    
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
data_filename = 'cpu_timing_nremap_two_stream_Tend10.mat';  % Re-define for consistency

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

% Set default interpreter to avoid LaTeX warnings
set(0, 'DefaultTextInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');

% Publication-quality plot with professional formatting
fig = figure('Name', 'CMM N_remap CPU Timing Analysis', 'Position', [100, 100, 1200, 500], ...
             'Units', 'pixels', 'PaperPositionMode', 'auto', 'Color', 'white');

% Set paper size for publication
set(fig, 'PaperUnits', 'inches', 'PaperSize', [12, 5]);

% Set publication font properties with larger text
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);

% Professional color scheme
colors = lines(num_tests);

% Subplot 1: CPU time per iteration vs simulation time
subplot(1, 2, 1);
hold on;
for i = 1:num_tests
    plot(results.time_arrays{i}, results.cpu_times{i}, ...
         'Color', colors(i, :), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('$N_{\\mathrm{remap}} = %d$', N_remap_values(i)));
end
hold off;
xlabel('Time $t$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('CPU Time per Iteration (s)', 'FontSize', 18);
title('CPU Time per Iteration vs Simulation Time', 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend('Location', 'best', 'FontSize', 14, 'Interpreter', 'latex', 'Box', 'on');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;
xlim([0, target_time]);

% Subplot 2: Mean CPU time per iteration for each N_remap
subplot(1, 2, 2);
% Create bars with colors matching the first subplot
bar_handles = bar(1:num_tests, results.mean_cpu_per_iter, 'EdgeColor', 'k', 'LineWidth', 1.5);
% Set individual bar colors to match the line plot colors
for i = 1:num_tests
    bar_handles.FaceColor = 'flat';
    bar_handles.CData(i,:) = colors(i,:);
end
xlabel('$N_{\mathrm{remap}}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Mean CPU Time per Iteration (s)', 'FontSize', 18);
title('Mean CPU Time per Iteration vs $N_{\mathrm{remap}}$', 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(gca, 'XTick', 1:num_tests, 'XTickLabel', N_remap_values);
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;

% Main title removed as requested

% Set background to white
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');

% Save plots with professional quality
plot_filename = 'cpu_timing_nremap_two_stream';
images_dir = './images';

% Create images directory if it doesn't exist
if ~exist(images_dir, 'dir')
    mkdir(images_dir);
    fprintf('Created %s directory\n', images_dir);
end

% High-quality export settings
png_path = fullfile(images_dir, [plot_filename '.png']);
fig_path = fullfile(images_dir, [plot_filename '.fig']);
eps_path = fullfile(images_dir, [plot_filename '.eps']);
pdf_path = fullfile(images_dir, [plot_filename '.pdf']);

% Save in multiple formats for publication
print(png_path, '-dpng', '-r300');
print(eps_path, '-depsc2', '-r300');
print(pdf_path, '-dpdf', '-r300');
saveas(gcf, fig_path);

fprintf('Professional plots saved as:\n');
fprintf('  PNG: %s\n', png_path);
fprintf('  EPS: %s\n', eps_path);
fprintf('  PDF: %s\n', pdf_path);
fprintf('  FIG: %s\n', fig_path);

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

% Reset default interpreters
set(0, 'DefaultTextInterpreter', 'tex');
set(0, 'DefaultLegendInterpreter', 'tex');
set(0, 'DefaultAxesTickLabelInterpreter', 'tex');