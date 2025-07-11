%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remapping frequency study for CMM method
% Compares CPU time vs time for different remapping frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup base parameters
PARAMS_two_stream;
params.Tend = 40;  % Set end time to 10 for quick testing
params.method = "CMM";

% Different remapping frequencies to test
N_remap_values = [10, 20, 40, 80, 10000000];  % 10000000 = no remapping (NuFi-like)
method_names = {'CMM (N_{remap}=10)', 'CMM (N_{remap}=20)', 'CMM (N_{remap}=40)', 'CMM (N_{remap}=80)', 'NuFi (no remap)'};

% Storage for results
results = cell(length(N_remap_values), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulations
fprintf('Running remapping frequency study...\n');
for i = 1:length(N_remap_values)
    fprintf('Running simulation %d/%d: N_remap = %d\n', i, length(N_remap_values), N_remap_values(i));
    
    % Set remapping frequency
    params.N_remap = N_remap_values(i);
    
    % Run simulation
    tic_total = tic();
    [params_result, data] = Sim(params);
    total_time = toc(tic_total);
    
    % Store results
    results{i}.params = params_result;
    results{i}.data = data;
    results{i}.total_time = total_time;
    results{i}.N_remap = N_remap_values(i);
    results{i}.method_name = method_names{i};
    
    fprintf('  Completed in %.2f seconds\n', total_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results
figure('Position', [100, 100, 1200, 800]);

% Plot 1: CPU time per iteration vs simulation time
subplot(2, 2, 1);
hold on;
colors = {'b', 'r', 'g', 'm', 'c'};
line_styles = {'-', '-', '-', '-', '--'};
for i = 1:length(results)
    plot(results{i}.params.time_array, results{i}.params.tcpu, ...
         'Color', colors{i}, 'LineStyle', line_styles{i}, 'LineWidth', 2, ...
         'DisplayName', results{i}.method_name);
end
xlabel('Simulation Time');
ylabel('CPU Time per Iteration (s)');
title('CPU Time per Iteration vs Simulation Time');
legend('Location', 'best');
grid on;

% Plot 2: Cumulative CPU time vs simulation time
subplot(2, 2, 2);
hold on;
for i = 1:length(results)
    cumulative_time = cumsum(results{i}.params.tcpu);
    plot(results{i}.params.time_array, cumulative_time, ...
         'Color', colors{i}, 'LineStyle', line_styles{i}, 'LineWidth', 2, ...
         'DisplayName', results{i}.method_name);
end
xlabel('Simulation Time');
ylabel('Cumulative CPU Time (s)');
title('Cumulative CPU Time vs Simulation Time');
legend('Location', 'best');
grid on;

% Plot 3: Total simulation time comparison
subplot(2, 2, 3);
total_times = zeros(length(results), 1);
for i = 1:length(results)
    total_times(i) = results{i}.total_time;
end
color_matrix = [0 0 1; 1 0 0; 0 1 0; 1 0 1; 0 1 1]; % blue, red, green, magenta, cyan
b = bar(total_times, 'FaceColor', 'flat');
b.CData = color_matrix;
set(gca, 'XTickLabel', method_names);
xtickangle(45);
ylabel('Total CPU Time (s)');
title('Total Simulation Time Comparison');
grid on;
% Add value labels on bars
for i = 1:length(total_times)
    text(i, total_times(i) + max(total_times)*0.02, sprintf('%.1fs', total_times(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

% Plot 4: Average CPU time per iteration
subplot(2, 2, 4);
avg_times = zeros(length(results), 1);
for i = 1:length(results)
    avg_times(i) = mean(results{i}.params.tcpu);
end
b = bar(avg_times, 'FaceColor', 'flat');
b.CData = color_matrix;
set(gca, 'XTickLabel', method_names);
xtickangle(45);
ylabel('Average CPU Time per Iteration (s)');
title('Average CPU Time per Iteration');
grid on;
% Add value labels on bars
for i = 1:length(avg_times)
    text(i, avg_times(i) + max(avg_times)*0.02, sprintf('%.3fs', avg_times(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

sgtitle('Remapping Frequency Study: CPU Performance Analysis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save results
save('remapping_study_results.mat', 'results', 'N_remap_values', 'method_names');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print summary
fprintf('\n=== REMAPPING STUDY SUMMARY ===\n');
fprintf('Simulation time: 0 to %.1f\n', params.Tend);
fprintf('%-25s | %12s | %15s\n', 'Method', 'Total Time', 'Avg Time/Iter');
fprintf('%-25s | %12s | %15s\n', repmat('-', 1, 25), repmat('-', 1, 12), repmat('-', 1, 15));
for i = 1:length(results)
    fprintf('%-25s | %10.2f s | %13.4f s\n', ...
            results{i}.method_name, ...
            results{i}.total_time, ...
            mean(results{i}.params.tcpu));
end
fprintf('\n');

% Find the optimal remapping frequency
cmm_times = zeros(length(results)-1, 1);
for i = 1:length(results)-1
    cmm_times(i) = results{i}.total_time;
end
[min_time, min_idx] = min(cmm_times); % Exclude NuFi from optimization
fprintf('Fastest CMM method: %s (%.2f s)\n', results{min_idx}.method_name, min_time);
fprintf('NuFi time: %.2f s\n', results{end}.total_time);
if min_time < results{end}.total_time
    speedup = results{end}.total_time / min_time;
    fprintf('Speedup over NuFi: %.2fx\n', speedup);
else
    slowdown = min_time / results{end}.total_time;
    fprintf('Slowdown vs NuFi: %.2fx\n', slowdown);
end