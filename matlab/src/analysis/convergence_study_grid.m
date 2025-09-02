%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid Convergence Study Script
% Measures how L∞ error decreases as grid size increases
% Grid sizes: 64, 128, 256, 512, 1024 (finest grid 1024 is benchmark)
% CMM method with N_remap = 5, measured at time t = 10
% Uses same interpolation as compose_maps_numerical (interp2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

% Configuration
grid_sizes = [64, 128, 256, 512, 1024];
target_time = 10.0;
N_remap_value = 5;

% Choose test case
fprintf('Available test cases:\n');
fprintf('1. Two Stream Instability\n');
fprintf('2. Landau Damping\n');
test_case_choice = input('Enter your choice (1 or 2): ');

% Choose benchmark method
fprintf('\nBenchmark method options:\n');
fprintf('1. CMM with 1024x1024 grid (default)\n');
fprintf('2. NuFi with 1024x1024 grid\n');
benchmark_choice = input('Enter benchmark choice (1 or 2, default=1): ');
if isempty(benchmark_choice)
    benchmark_choice = 1;
end

switch test_case_choice
    case 1
        test_case_name = "two_stream";
        params_script = "PARAMS_two_stream";
        fprintf('Selected: Two Stream Instability\n');
    case 2
        test_case_name = "landau_damping";
        params_script = "PARAMS_landau_damping";
        fprintf('Selected: Landau Damping\n');
    otherwise
        error('Invalid choice. Please run the script again and choose 1 or 2.');
end

switch benchmark_choice
    case 1
        benchmark_method = "CMM";
        fprintf('Selected benchmark: CMM with 1024x1024 grid\n');
    case 2
        benchmark_method = "NuFi";
        fprintf('Selected benchmark: NuFi with 1024x1024 grid\n');
    otherwise
        error('Invalid benchmark choice. Please run the script again and choose 1 or 2.');
end

fprintf('\n=== GRID CONVERGENCE STUDY (%s) ===\n', test_case_name);
fprintf('Grid sizes: %s\n', mat2str(grid_sizes));
fprintf('N_remap = %d, Target time = %.1f\n', N_remap_value, target_time);
fprintf('Benchmark: %s with %dx%d grid\n', benchmark_method, grid_sizes(end), grid_sizes(end));

% Initialize storage for results
distribution_functions = cell(length(grid_sizes), 1);
simulation_times = zeros(length(grid_sizes), 1);
grid_params = cell(length(grid_sizes), 1);

%% Run simulations for each grid size
for i = 1:length(grid_sizes)
    grid_size = grid_sizes(i);
    fprintf('Running simulation %d/%d: Grid %dx%d...\n', i, length(grid_sizes), grid_size, grid_size);
    
    % Save current state and clear to avoid parameter conflicts
    temp_filename = sprintf('temp_convergence_%d.mat', i);
    save('temp_convergence_config.mat', 'test_case_name', 'params_script', 'grid_sizes', 'target_time', 'N_remap_value', 'benchmark_method', 'i');
    
    % Clear workspace and reload basics (following measure_remap_error.m pattern)
    clear all
    addpath(genpath('./src/'),genpath('./params/'))
    load('temp_convergence_config.mat');
    
    % Get current grid size
    grid_size = grid_sizes(i);
    
    % Load fresh parameters for current grid
    eval(params_script);
    
    % Determine which method to use for this grid
    if i == length(grid_sizes)  % Last grid (1024x1024) - use benchmark method
        params.method = benchmark_method;
        fprintf('  Using %s method for benchmark grid\n', benchmark_method);
    else  % All other grids - use CMM
        params.method = "CMM";
        fprintf('  Using CMM method for grid %dx%d\n', grid_size, grid_size);
    end
    
    params.N_remap = N_remap_value;
    params.Tend = target_time;
    
    % Set grid resolution
    params.Nx = grid_size;
    params.Nv = grid_size;
    
    % Run simulation
    [sim_params, sim_data] = Sim(params);
    
    % Find time step closest to target time
    [~, time_idx] = min(abs(sim_data.time - target_time));
    actual_time = sim_data.time(time_idx);
    
    % Save results for this grid size
    grid_results = struct();
    grid_results.distribution_function = sim_data.fs(:,:,time_idx);
    grid_results.simulation_time = actual_time;
    grid_results.grid_params = sim_params;
    grid_results.grid_size = grid_size;
    
    temp_filename = sprintf('temp_convergence_%d.mat', i);
    save(temp_filename, 'grid_results');
    
    method_used = params.method;
    fprintf('  Completed: Grid %dx%d (%s), actual time = %.4f\n', grid_size, grid_size, method_used, actual_time);
end

% Load all results back
fprintf('\nLoading all simulation results...\n');
clear all
addpath(genpath('./src/'),genpath('./params/'))
load('temp_convergence_config.mat');

% Initialize storage arrays
distribution_functions = cell(length(grid_sizes), 1);
simulation_times = zeros(length(grid_sizes), 1);
grid_params = cell(length(grid_sizes), 1);

% Load results from temporary files
for i = 1:length(grid_sizes)
    temp_filename = sprintf('temp_convergence_%d.mat', i);
    temp_data = load(temp_filename);
    
    distribution_functions{i} = temp_data.grid_results.distribution_function;
    simulation_times(i) = temp_data.grid_results.simulation_time;
    grid_params{i} = temp_data.grid_results.grid_params;
    
    fprintf('Loaded results for grid %dx%d\n', grid_sizes(i), grid_sizes(i));
end

%% Extract benchmark (finest grid) solution
benchmark_f = distribution_functions{end};  % 1024x1024 grid
benchmark_params = grid_params{end};

% Get benchmark grid coordinates
x_bench = benchmark_params.grids(1).x;
v_bench = benchmark_params.grids(1).v;

fprintf('\nBenchmark grid: %dx%d (%s method)\n', length(x_bench), length(v_bench), benchmark_method);

%% Compute errors using interpolation method 
L_inf_errors = zeros(length(grid_sizes)-1, 1);
L_inf_rel_errors = zeros(length(grid_sizes)-1, 1);

fprintf('\nError Analysis (MATLAB interp2 Method):\n');
fprintf('%-10s %-15s %-15s %-10s %-15s\n', 'Grid', 'L∞ Absolute', 'L∞ Relative', 'Time', 'Method');
fprintf('%-10s %-15s %-15s %-10s %-15s\n', '----', '-----------', '-----------', '----', '------');

 

% Create benchmark grid meshgrid once
[X_bench, V_bench] = meshgrid(x_bench, v_bench);

for i = 1:length(grid_sizes)-1
    grid_size = grid_sizes(i);
    coarse_f = distribution_functions{i};
    coarse_params = grid_params{i};
    
    % Get coarse grid coordinates
    x_coarse = coarse_params.grids(1).x;
    v_coarse = coarse_params.grids(1).v;
    [X_coarse, V_coarse] = meshgrid(x_coarse, v_coarse);
    
    % Interpolate benchmark solution to coarse grid coordinates
    f_bench_interp = interp2(X_bench, V_bench, benchmark_f, X_coarse, V_coarse, 'cubic');
    
    % Calculate errors
    error_abs = abs(coarse_f - f_bench_interp);
    L_inf_abs = max(error_abs(:));
    L_inf_rel = L_inf_abs / max(abs(f_bench_interp(:)));
    
    L_inf_errors(i) = L_inf_abs;
    L_inf_rel_errors(i) = L_inf_rel;
    
    method_str = 'Cubic-interp2';
    
    fprintf('%-10s %-15.6e %-15.6e %-10.4f %-15s\n', sprintf('%dx%d', grid_size, grid_size), ...
            L_inf_abs, L_inf_rel, simulation_times(i), method_str);
end

% If using interpolation gives very small errors, also try subsampling for comparison
if all(L_inf_errors < 1e-12)
    fprintf('\n⚠️  Interpolation errors very small. Comparing with subsampling method:\n');
    fprintf('%-10s %-15s %-15s %-15s\n', 'Grid', 'Interpolation', 'Subsampling', 'Ratio');
    fprintf('%-10s %-15s %-15s %-15s\n', '----', '------------', '-----------', '-----');
    
    for i = 1:length(grid_sizes)-1
        grid_size = grid_sizes(i);
        coarse_f = distribution_functions{i};
        
        % Subsampling method for comparison
        ratio = grid_sizes(end) / grid_size;
        if mod(ratio, 1) == 0
            f_bench_subsampled = benchmark_f(1:ratio:end, 1:ratio:end);
            subsample_error = max(abs(coarse_f(:) - f_bench_subsampled(:)));
            error_ratio = L_inf_errors(i) / max(subsample_error, eps);
            
            fprintf('%-10s %-15.6e %-15.6e %-15.2f\n', sprintf('%dx%d', grid_size, grid_size), ...
                    L_inf_errors(i), subsample_error, error_ratio);
        end
    end
    
    fprintf('\nNote: Large ratio suggests interpolation may be too accurate for convergence study\n');
end

% If all errors are still zero, try alternative approach
if all(L_inf_errors < 1e-14)
    fprintf('\n⚠️  All errors are near machine precision. Trying alternative approaches...\n');
    
    % Alternative: Compare adjacent grids in the sequence
    fprintf('\nAlternative Analysis (Adjacent Grid Comparison):\n');
    fprintf('%-15s %-15s %-15s\n', 'Comparison', 'L∞ Absolute', 'L∞ Relative');
    fprintf('%-15s %-15s %-15s\n', '----------', '-----------', '-----------');
    
    for i = 1:length(grid_sizes)-1
        coarse_f = distribution_functions{i};
        fine_f = distribution_functions{i+1};
        
        % Subsample fine grid to coarse resolution
        ratio = grid_sizes(i+1) / grid_sizes(i);
        fine_subsampled = fine_f(1:ratio:end, 1:ratio:end);
        
        alt_error = abs(coarse_f - fine_subsampled);
        alt_L_inf = max(alt_error(:));
        alt_L_inf_rel = alt_L_inf / max(abs(fine_subsampled(:)));
        
        fprintf('%-15s %-15.6e %-15.6e\n', ...
                sprintf('%dx%d vs %dx%d', grid_sizes(i), grid_sizes(i), grid_sizes(i+1), grid_sizes(i+1)), ...
                alt_L_inf, alt_L_inf_rel);
        
        % Update main arrays if this gives better results
        if i < length(L_inf_errors) && alt_L_inf > L_inf_errors(i)
            L_inf_errors(i) = alt_L_inf;
            L_inf_rel_errors(i) = alt_L_inf_rel;
        end
    end
    
    fprintf('\nSuggestions:\n');
    fprintf('• Try later time (t > 5) when instabilities are more developed\n');
    fprintf('• Use less frequent remapping (N_remap = 1000) to see grid effects\n');
    fprintf('• Consider different test case parameters\n');
end

%% Calculate convergence rates
fprintf('\nConvergence Analysis:\n');
fprintf('%-15s %-15s %-15s %-15s\n', 'Grid Pair', 'L∞ Error', 'h', 'Conv. Rate');
fprintf('%-15s %-15s %-15s %-15s\n', '---------', '--------', '-', '---------');

convergence_rates = zeros(length(L_inf_errors)-1, 1);
for i = 1:length(L_inf_errors)-1
    h_coarse = 1/grid_sizes(i);      % Grid spacing for coarser grid
    h_fine = 1/grid_sizes(i+1);     % Grid spacing for finer grid
    
    error_coarse = L_inf_errors(i);
    error_fine = L_inf_errors(i+1);
    
    % Convergence rate: log(error_coarse/error_fine) / log(h_coarse/h_fine)
    conv_rate = log(error_coarse/error_fine) / log(h_coarse/h_fine);
    convergence_rates(i) = conv_rate;
    
    fprintf('%-15s %-15.6e %-15.6e %-15.2f\n', ...
            sprintf('%d→%d', grid_sizes(i), grid_sizes(i+1)), ...
            error_fine, h_fine, conv_rate);
end

mean_conv_rate = mean(convergence_rates);
fprintf('\nMean convergence rate (pairwise): %.2f\n', mean_conv_rate);

%% Regression-based convergence analysis
fprintf('\nRegression-based Convergence Analysis:\n');

% Filter out any zero or negative errors for log fitting
valid_indices = L_inf_errors > 0;
if sum(valid_indices) >= 2
    h_values = 1./grid_sizes(1:end-1);  % Grid spacing values
    valid_h = h_values(valid_indices);
    valid_errors = L_inf_errors(valid_indices);
    
    % Fit log(error) = log(C) + p*log(h), where p is the convergence rate
    % This gives: error = C * h^p
    log_h = log(valid_h);
    log_error = log(valid_errors);
    
    % Linear regression: log_error = intercept + slope * log_h
    coeffs = polyfit(log_h, log_error, 1);
    regression_conv_rate = coeffs(1);  % slope = convergence rate
    log_C = coeffs(2);                 % intercept = log(C)
    C_constant = exp(log_C);
    
    % Calculate R-squared for goodness of fit
    log_error_fitted = polyval(coeffs, log_h);
    SS_res = sum((log_error - log_error_fitted).^2);
    SS_tot = sum((log_error - mean(log_error)).^2);
    R_squared = 1 - SS_res/SS_tot;
    
    fprintf('Fitted model: Error = %.3e * h^%.2f\n', C_constant, regression_conv_rate);
    fprintf('Regression convergence rate: %.2f\n', regression_conv_rate);
    fprintf('R-squared (goodness of fit): %.4f\n', R_squared);
    
    % Compare methods
    fprintf('\nComparison of Methods:\n');
    fprintf('%-25s %-15s\n', 'Method', 'Conv. Rate');
    fprintf('%-25s %-15s\n', '------', '----------');
    fprintf('%-25s %-15.2f\n', 'Pairwise (mean)', mean_conv_rate);
    fprintf('%-25s %-15.2f\n', 'Regression fit', regression_conv_rate);
    
    if abs(mean_conv_rate - regression_conv_rate) < 0.5
        fprintf('✓ Methods agree well (difference < 0.5)\n');
    else
        fprintf('! Methods differ significantly (%.2f vs %.2f)\n', mean_conv_rate, regression_conv_rate);
        if R_squared > 0.95
            fprintf('  Regression fit is excellent (R² > 0.95) - prefer regression rate\n');
        elseif R_squared < 0.8
            fprintf('  Regression fit is poor (R² < 0.8) - prefer pairwise rate\n');
        else
            fprintf('  Regression fit is reasonable - both rates are valid\n');
        end
    end
    
else
    regression_conv_rate = NaN;
    R_squared = NaN;
    fprintf('⚠️  Not enough valid data points for regression analysis\n');
end

%% Visualizations

% 1. Convergence plot
figure('Name', sprintf('Grid Convergence Study (%s)', test_case_name), 'Position', [100, 100, 1200, 800]);

subplot(2,3,1);
loglog(grid_sizes(1:end-1), L_inf_errors, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
% Add regression fit line if available
if exist('regression_conv_rate', 'var') && ~isnan(regression_conv_rate)
    h_fit = 1./grid_sizes(1:end-1);
    error_fit = C_constant * (h_fit).^regression_conv_rate;
    loglog(grid_sizes(1:end-1), error_fit, '--r', 'LineWidth', 1.5);
    legend('Computed', sprintf('Fit: h^{%.1f}', regression_conv_rate), 'Location', 'best');
else
    legend('Computed', 'Location', 'best');
end
hold off;
xlabel('Grid Size');
ylabel('L∞ Error');
title('Convergence: L∞ Error vs Grid Size');
grid on;

subplot(2,3,2);
plot(1:length(convergence_rates), convergence_rates, 's-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Refinement Level');
ylabel('Convergence Rate');
title('Convergence Rate');
yline(mean_conv_rate, '--r', sprintf('Mean = %.2f', mean_conv_rate));
grid on;

% 2. Distribution function comparison (show 3 representative grids)
subplot(2,3,4);
imagesc(distribution_functions{1});  % 64x64
title(sprintf('64×64 Grid (t=%.1f)', target_time));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,5);
imagesc(distribution_functions{3});  % 256x256
title(sprintf('256×256 Grid (t=%.1f)', target_time));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,6);
imagesc(benchmark_f);  
title(sprintf('1024×1024 Grid - %s Benchmark (t=%.1f)', benchmark_method, target_time));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

% 3. Error visualization
subplot(2,3,3);
% Interpolate benchmark to 64x64 grid for visualization
x_64 = grid_params{1}.grids(1).x;
v_64 = grid_params{1}.grids(1).v;
[X_64, V_64] = meshgrid(x_64, v_64);
f_bench_64 = interp2(X_bench, V_bench, benchmark_f, X_64, V_64, 'cubic');
error_vis = abs(distribution_functions{1} - f_bench_64);
imagesc(error_vis);
title(sprintf('Error: %dx%d vs Benchmark', grid_sizes(1), grid_sizes(1)));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

sgtitle(sprintf('Grid Convergence Study: %s (N_{remap}=%d, t=%.1f, %s benchmark)', ...
                test_case_name, N_remap_value, target_time, benchmark_method));

%% Save results
results = struct();
results.test_case = test_case_name;
results.grid_sizes = grid_sizes;
results.target_time = target_time;
results.N_remap = N_remap_value;
results.benchmark_method = benchmark_method;
results.simulation_times = simulation_times;
results.L_inf_errors = L_inf_errors;
results.L_inf_rel_errors = L_inf_rel_errors;
results.convergence_rates = convergence_rates;
results.mean_convergence_rate = mean_conv_rate;
if exist('regression_conv_rate', 'var')
    results.regression_convergence_rate = regression_conv_rate;
    results.R_squared = R_squared;
end
results.distribution_functions = distribution_functions;
results.grid_params = grid_params;

results_filename = sprintf('convergence_study_results_%s_t%.1f_Nremap%d_%s.mat', ...
                          test_case_name, target_time, N_remap_value, benchmark_method);
save(results_filename, 'results');
fprintf('\nResults saved to %s\n', results_filename);

% Clean up temporary files
delete('temp_convergence_config.mat');
for i = 1:length(grid_sizes)
    temp_filename = sprintf('temp_convergence_%d.mat', i);
    if exist(temp_filename, 'file')
        delete(temp_filename);
    end
end
% Clean up NuFi temporary files
if exist('temp_before_nufi.mat', 'file')
    delete('temp_before_nufi.mat');
end
if exist('temp_nufi_results.mat', 'file')
    delete('temp_nufi_results.mat');
end
fprintf('Temporary files cleaned up.\n');

% Summary
fprintf('\n=== CONVERGENCE STUDY SUMMARY ===\n');
fprintf('Test case: %s\n', test_case_name);
fprintf('Time: t = %.1f, N_remap = %d\n', target_time, N_remap_value);
fprintf('Grids tested: %s\n', mat2str(grid_sizes));
fprintf('Benchmark: %s with %dx%d grid\n', benchmark_method, grid_sizes(end), grid_sizes(end));
fprintf('Mean convergence rate: %.2f\n', mean_conv_rate);
