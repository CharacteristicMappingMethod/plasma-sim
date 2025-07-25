%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CMM Error Analysis - Fresh Simulations
% Runs simulations with different N_remap values and computes errors
% Benchmark: N_remap = 100000 (equivalent to NuFi)
% Error analysis at t = 5 and t = 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup parameters
PARAMS_two_stream;
params.Tend = 20;  % Run until t=20
params.method = "CMM";
params.dt_save = 5;  % Save every 5 time units to capture target times

% Target times for error analysis  
target_times = [5, 15, 20];  % Use times that align with save points

% Remapping frequencies to test
N_remap_test = [10, 20, 40, 80];
N_remap_benchmark = 100000;  % Benchmark (equivalent to NuFi)

fprintf('=== CMM ERROR ANALYSIS ===\n');
fprintf('Running simulations for error analysis...\n');
fprintf('Benchmark: N_remap = %d\n', N_remap_benchmark);
fprintf('Test cases: %s\n', mat2str(N_remap_test));
fprintf('Target times: %s\n', mat2str(target_times));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run benchmark simulation
fprintf('\nRunning benchmark simulation (N_remap = %d)...\n', N_remap_benchmark);
params.N_remap = N_remap_benchmark;
tic_bench = tic();
[params_benchmark, data_benchmark] = Sim(params);
bench_time = toc(tic_bench);
fprintf('Benchmark completed in %.2f seconds\n', bench_time);

% Find time indices for benchmark
benchmark_time_indices = zeros(length(target_times), 1);
for t_idx = 1:length(target_times)
    [~, benchmark_time_indices(t_idx)] = min(abs(data_benchmark.time - target_times(t_idx)));
    actual_time = data_benchmark.time(benchmark_time_indices(t_idx));
    fprintf('Benchmark t=%.1f -> actual t=%.3f (index %d)\n', target_times(t_idx), actual_time, benchmark_time_indices(t_idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run test simulations and compute errors
n_methods = length(N_remap_test);
n_times = length(target_times);

% Storage for errors
errors = struct();
errors.electric_field = zeros(n_methods, n_times);
errors.distribution_function = zeros(n_methods, n_times);
errors.potential_energy = zeros(n_methods, n_times);

% Storage for simulation data
test_results = cell(n_methods, 1);

for method_idx = 1:n_methods
    N_remap_current = N_remap_test(method_idx);
    fprintf('\nRunning test simulation %d/%d (N_remap = %d)...\n', method_idx, n_methods, N_remap_current);
    
    % Set remapping frequency
    params.N_remap = N_remap_current;
    
    % Run simulation
    tic_test = tic();
    [params_test, data_test] = Sim(params);
    test_time = toc(tic_test);
    fprintf('Test simulation completed in %.2f seconds\n', test_time);
    
    % Store results
    test_results{method_idx}.params = params_test;
    test_results{method_idx}.data = data_test;
    test_results{method_idx}.N_remap = N_remap_current;
    test_results{method_idx}.simulation_time = test_time;
    
    % Compute errors at each target time
    for t_idx = 1:n_times
        target_time = target_times(t_idx);
        
        % Find corresponding time index in test data
        [~, test_time_idx] = min(abs(data_test.time - target_time));
        actual_test_time = data_test.time(test_time_idx);
        
        % Get benchmark data at this time
        bench_time_idx = benchmark_time_indices(t_idx);
        benchmark_E = data_benchmark.Efield(:, bench_time_idx);
        benchmark_f = data_benchmark.fs(:, :, bench_time_idx, 1); % First species
        
        % Get test data at this time
        test_E = data_test.Efield(:, test_time_idx);
        test_f = data_test.fs(:, :, test_time_idx, 1);
        
        % Grid spacing
        dx = params_benchmark.grids(1).dx;
        dv = params_benchmark.grids(1).dv;
        
        % Compute benchmark potential energy
        benchmark_Epot = 0.5 * sum(benchmark_E.^2) * dx;
        
        % Compute test potential energy
        test_Epot = 0.5 * sum(test_E.^2) * dx;
        
        % 1. Electric field L2 relative error
        E_error = sqrt(sum((test_E - benchmark_E).^2) * dx) / sqrt(sum(benchmark_E.^2) * dx);
        errors.electric_field(method_idx, t_idx) = E_error;
        
        % 2. Distribution function L2 relative error
        f_error = sqrt(sum(sum((test_f - benchmark_f).^2)) * dx * dv) / ...
                  sqrt(sum(sum(benchmark_f.^2)) * dx * dv);
        errors.distribution_function(method_idx, t_idx) = f_error;
        
        % 3. Potential energy relative error
        Epot_error = abs(test_Epot - benchmark_Epot) / abs(benchmark_Epot);
        errors.potential_energy(method_idx, t_idx) = Epot_error;
        
        fprintf('  t=%.1f (%.3f): E_field=%.2e, f=%.2e, Epot=%.2e\n', ...
                target_time, actual_test_time, E_error, f_error, Epot_error);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results
figure('Position', [100, 100, 1400, 900]);

% Colors and markers for different times
colors = {'b', 'r', 'g'};
markers = {'o', 's', '^'};
time_labels = {'t = 5', 't = 15', 't = 50'};

% Plot 1: Electric Field Error
subplot(2, 2, 1);
hold on;
for t_idx = 1:n_times
    loglog(N_remap_test, errors.electric_field(:, t_idx), ...
           [colors{t_idx} markers{t_idx} '-'], 'LineWidth', 2, 'MarkerSize', 8, ...
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
    loglog(N_remap_test, errors.distribution_function(:, t_idx), ...
           [colors{t_idx} markers{t_idx} '-'], 'LineWidth', 2, 'MarkerSize', 8, ...
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
    loglog(N_remap_test, errors.potential_energy(:, t_idx), ...
           [colors{t_idx} markers{t_idx} '-'], 'LineWidth', 2, 'MarkerSize', 8, ...
           'DisplayName', time_labels{t_idx});
end
xlabel('N_{remap}');
ylabel('Potential Energy Relative Error');
title('Potential Energy Error vs Remapping Frequency');
legend('Location', 'best');
grid on;

% Plot 4: Combined error comparison at t=50
subplot(2, 2, 4);
hold on;
t_idx = 3; % t = 50
loglog(N_remap_test, errors.electric_field(:, t_idx), 'b o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Electric Field');
loglog(N_remap_test, errors.distribution_function(:, t_idx), 'r s-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Distribution Function');
loglog(N_remap_test, errors.potential_energy(:, t_idx), 'g ^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Potential Energy');
xlabel('N_{remap}');
ylabel('Relative Error');
title('All Errors at t = 50');
legend('Location', 'best');
grid on;

sgtitle('CMM Error Analysis: Fresh Simulations');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print summary
fprintf('\n=== ERROR ANALYSIS SUMMARY ===\n');
fprintf('Benchmark: N_remap = %d (%.2f seconds)\n', N_remap_benchmark, bench_time);
fprintf('\nError Summary:\n');

for t_idx = 1:n_times
    fprintf('\nErrors at t = %.1f:\n', target_times(t_idx));
    fprintf('%-8s | %12s | %15s | %15s\n', 'N_remap', 'E-field', 'Dist. Func.', 'Pot. Energy');
    fprintf('%-8s | %12s | %15s | %15s\n', repmat('-', 1, 8), repmat('-', 1, 12), repmat('-', 1, 15), repmat('-', 1, 15));
    
    for method_idx = 1:n_methods
        fprintf('%-8d | %10.2e | %13.2e | %13.2e\n', ...
                N_remap_test(method_idx), ...
                errors.electric_field(method_idx, t_idx), ...
                errors.distribution_function(method_idx, t_idx), ...
                errors.potential_energy(method_idx, t_idx));
    end
end

% Print simulation times
fprintf('\nSimulation Times:\n');
fprintf('%-8s | %12s\n', 'N_remap', 'Time (s)');
fprintf('%-8s | %12s\n', repmat('-', 1, 8), repmat('-', 1, 12));
fprintf('%-8d | %10.2f\n', N_remap_benchmark, bench_time);
for method_idx = 1:n_methods
    fprintf('%-8d | %10.2f\n', N_remap_test(method_idx), test_results{method_idx}.simulation_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save results
results_struct = struct();
results_struct.N_remap_test = N_remap_test;
results_struct.N_remap_benchmark = N_remap_benchmark;
results_struct.target_times = target_times;
results_struct.errors = errors;
results_struct.test_results = test_results;
results_struct.benchmark_data = data_benchmark;
results_struct.benchmark_params = params_benchmark;
results_struct.benchmark_time = bench_time;

save('cmm_error_analysis_fresh.mat', 'results_struct');
fprintf('\nResults saved to: cmm_error_analysis_fresh.mat\n');