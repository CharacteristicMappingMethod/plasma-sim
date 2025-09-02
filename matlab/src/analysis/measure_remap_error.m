%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to measure error between CMM N_remap = 10 and N_remap = 1000
% Benchmark: NuFi method, Test cases: CMM with different N_remap
% Test cases: Two Stream or Landau Damping
% Error measured at time = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

% Choose test case
fprintf('Available test cases:\n');
fprintf('1. Two Stream Instability\n');
fprintf('2. Landau Damping\n');
test_case_choice = input('Enter your choice (1 or 2): ');

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

% Check for existing NuFi benchmark data
benchmark_filename = sprintf('benchmark_NuFi_%s_Tend10.mat', test_case_name);
if exist(benchmark_filename, 'file')
    fprintf('Found existing NuFi benchmark data for %s. Loading...\n', test_case_name);
    load(benchmark_filename, 'params_benchmark', 'data_benchmark');
    fprintf('Loaded benchmark data from %s\n', benchmark_filename);
else
    % Run simulation with NuFi method (benchmark)
    fprintf('No existing benchmark data found. Running NuFi simulation (%s)...\n', test_case_name);
    eval(params_script);
    params.method = "NuFi";
    params.Tend = 10;
    [params_benchmark, data_benchmark] = Sim(params);
    
    % Save benchmark data persistently
    save(benchmark_filename, 'params_benchmark', 'data_benchmark');
    fprintf('Benchmark data saved to %s\n', benchmark_filename);
end

% Save temporary config and clear everything
save('temp_config.mat', 'test_case_name', 'params_script', 'benchmark_filename');
clear all
addpath(genpath('./src/'),genpath('./params/'))
load('temp_config.mat');

% Check for existing CMM N_remap=1000 data
cmm1000_filename = sprintf('benchmark_CMM_Nremap1000_%s_Tend10.mat', test_case_name);
if exist(cmm1000_filename, 'file')
    fprintf('Found existing CMM N_remap=1000 data for %s. Loading...\n', test_case_name);
    load(cmm1000_filename, 'params_cmm1000', 'data_cmm1000');
    fprintf('Loaded CMM 1000 data from %s\n', cmm1000_filename);
else
    % Run simulation with CMM N_remap = 1000
    fprintf('No existing CMM 1000 data found. Running CMM simulation with N_remap = 1000 (%s)...\n', test_case_name);
    eval(params_script);
    params.N_remap = 1000;
    params.method = "CMM";
    params.Tend = 10;
    [params_cmm1000, data_cmm1000] = Sim(params);
    
    % Save CMM 1000 data persistently
    save(cmm1000_filename, 'params_cmm1000', 'data_cmm1000');
    fprintf('CMM 1000 data saved to %s\n', cmm1000_filename);
end

% Update temporary config and clear everything
save('temp_config.mat', 'test_case_name', 'params_script', 'benchmark_filename', 'cmm1000_filename');
clear all
addpath(genpath('./src/'),genpath('./params/'))
load('temp_config.mat');

% Run simulation with CMM N_remap = 20
fprintf('Running CMM simulation with N_remap = 20 (%s)...\n', test_case_name);
eval(params_script);
params.N_remap = 20;
params.method = "CMM";
params.Tend = 10;
[params_cmm10, data_cmm10] = Sim(params);

% Extract map stack data for N_remap = 20 case
map_stack_cmm10 = params_cmm10.Map_stack;
Nmaps_cmm10 = params_cmm10.Nmaps;
fprintf('CMM N_remap=20 generated %d maps in the stack\n', Nmaps_cmm10);

% Load benchmark and CMM 1000 data back
load(benchmark_filename, 'params_benchmark', 'data_benchmark');
load(cmm1000_filename, 'params_cmm1000', 'data_cmm1000');

% Find the time step closest to t = 10 for all methods
target_time = 10;
time_benchmark = data_benchmark.time;
time_cmm1000 = data_cmm1000.time;
time_cmm10 = data_cmm10.time;

[~, idx_benchmark] = min(abs(time_benchmark - target_time));
[~, idx_cmm1000] = min(abs(time_cmm1000 - target_time));
[~, idx_cmm10] = min(abs(time_cmm10 - target_time));

actual_time_benchmark = time_benchmark(idx_benchmark);
actual_time_cmm1000 = time_cmm1000(idx_cmm1000);
actual_time_cmm10 = time_cmm10(idx_cmm10);

fprintf('NuFi (benchmark) time: %.4f\n', actual_time_benchmark);
fprintf('CMM N_remap=1000 time: %.4f\n', actual_time_cmm1000);
fprintf('CMM N_remap=20 time:   %.4f\n', actual_time_cmm10);

% Extract distribution functions at time ≈ 10
f_benchmark = data_benchmark.fs(:,:,idx_benchmark);
f_cmm1000 = data_cmm1000.fs(:,:,idx_cmm1000);
f_cmm10 = data_cmm10.fs(:,:,idx_cmm10);

% Calculate error metrics for CMM N_remap=1000 vs NuFi
% Relative errors
error_L2_rel_1000 = norm(f_cmm1000 - f_benchmark, 'fro') / norm(f_benchmark, 'fro');
error_Linf_rel_1000 = max(abs(f_cmm1000(:) - f_benchmark(:))) / max(abs(f_benchmark(:)));
% Absolute errors
error_L2_abs_1000 = norm(f_cmm1000 - f_benchmark, 'fro');
error_Linf_abs_1000 = max(abs(f_cmm1000(:) - f_benchmark(:)));

% Calculate error metrics for CMM N_remap=20 vs NuFi
% Relative errors
error_L2_rel_10 = norm(f_cmm10 - f_benchmark, 'fro') / norm(f_benchmark, 'fro');
error_Linf_rel_10 = max(abs(f_cmm10(:) - f_benchmark(:))) / max(abs(f_benchmark(:)));
% Absolute errors
error_L2_abs_10 = norm(f_cmm10 - f_benchmark, 'fro');
error_Linf_abs_10 = max(abs(f_cmm10(:) - f_benchmark(:)));

% Display results
fprintf('\n=== ERROR ANALYSIS RESULTS (%s) ===\n', test_case_name);
fprintf('Benchmark: NuFi method\n');
fprintf('Error measured at time ≈ %.1f\n\n', target_time);

fprintf('CMM N_remap=1000 vs NuFi:\n');
fprintf('  Absolute L2 error:   %.6e\n', error_L2_abs_1000);
fprintf('  Absolute L∞ error:   %.6e\n', error_Linf_abs_1000);
fprintf('  Relative L2 error:   %.6e\n', error_L2_rel_1000);
fprintf('  Relative L∞ error:   %.6e\n\n', error_Linf_rel_1000);

fprintf('CMM N_remap=20 vs NuFi:\n');
fprintf('  Absolute L2 error:   %.6e\n', error_L2_abs_10);
fprintf('  Absolute L∞ error:   %.6e\n', error_Linf_abs_10);
fprintf('  Relative L2 error:   %.6e\n', error_L2_rel_10);
fprintf('  Relative L∞ error:   %.6e\n', error_Linf_rel_10);
fprintf('================================\n');

% Plot individual maps in the map stack for CMM N_remap = 20
if Nmaps_cmm10 > 0
    fprintf('\nPlotting individual maps from CMM N_remap=20 map stack (%s)...\n', test_case_name);
    
    % Get grid information for plotting
    grid_info = params_cmm10.grids(1); % Assuming we want to plot for first species
    
    % Create a figure for all maps
    figure('Name', sprintf('CMM Map Stack (N_remap=20, %s)', test_case_name), 'Position', [100, 100, 1200, 800]);
    
    % Calculate subplot layout (try to make it roughly square)
    n_cols = ceil(sqrt(Nmaps_cmm10 * 2)); % *2 because we plot X and V components
    n_rows = ceil((Nmaps_cmm10 * 2) / n_cols);
    
    for map_idx = 1:Nmaps_cmm10
        % Extract X and V components for species 1 (index 1)
        Map_X = map_stack_cmm10(:, :, 1, 1, map_idx); % [Nv, Nx, X_component, species_1, map_idx]
        Map_V = map_stack_cmm10(:, :, 2, 1, map_idx); % [Nv, Nx, V_component, species_1, map_idx]
        
        % Plot X component
        subplot(n_rows, n_cols, (map_idx-1)*2 + 1);
        pcolor(grid_info.X, grid_info.V, Map_X);
        shading flat;
        colorbar;
        title(sprintf('Map %d: X-component', map_idx));
        xlabel('x');
        ylabel('v');
        
        % Plot V component
        subplot(n_rows, n_cols, (map_idx-1)*2 + 2);
        pcolor(grid_info.X, grid_info.V, Map_V);
        shading flat;
        colorbar;
        title(sprintf('Map %d: V-component', map_idx));
        xlabel('x');
        ylabel('v');
        
        fprintf('  Plotted Map %d: X range [%.3f, %.3f], V range [%.3f, %.3f]\n', ...
                map_idx, min(Map_X,[],'all'), max(Map_X,[],'all'), ...
                min(Map_V,[],'all'), max(Map_V,[],'all'));
    end
    
    sgtitle(sprintf('CMM Map Stack (N_{remap}=20, %s): %d Maps Generated', test_case_name, Nmaps_cmm10));
    
    % Also create displacement plots (Map - Identity)
    figure('Name', sprintf('CMM Map Displacements (N_remap=20, %s)', test_case_name), 'Position', [150, 150, 1200, 800]);
    
    for map_idx = 1:Nmaps_cmm10
        % Extract X and V components
        Map_X = map_stack_cmm10(:, :, 1, 1, map_idx);
        Map_V = map_stack_cmm10(:, :, 2, 1, map_idx);
        
        % Calculate displacements from identity
        Disp_X = Map_X - grid_info.X;
        Disp_V = Map_V - grid_info.V;
        
        % Plot X displacement
        subplot(n_rows, n_cols, (map_idx-1)*2 + 1);
        pcolor(grid_info.X, grid_info.V, Disp_X);
        shading flat;
        colorbar;
        title(sprintf('Map %d: ΔX = X - X_{id}', map_idx));
        xlabel('x');
        ylabel('v');
        
        % Plot V displacement
        subplot(n_rows, n_cols, (map_idx-1)*2 + 2);
        pcolor(grid_info.X, grid_info.V, Disp_V);
        shading flat;
        colorbar;
        title(sprintf('Map %d: ΔV = V - V_{id}', map_idx));
        xlabel('x');
        ylabel('v');
    end
    
    sgtitle(sprintf('CMM Map Displacements (N_{remap}=20, %s): %d Maps Generated', test_case_name, Nmaps_cmm10));
end

% Plot comparison
figure(1);
subplot(2,3,1);
imagesc(f_benchmark);
title(sprintf('NuFi (Benchmark) - %s', test_case_name));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,2);
imagesc(f_cmm1000);
title(sprintf('CMM N_{remap} = 1000 - %s', test_case_name));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,3);
imagesc(f_cmm10);
title(sprintf('CMM N_{remap} = 20 - %s', test_case_name));
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,4);
imagesc(abs(f_cmm1000 - f_benchmark));
title('Error: CMM 1000 vs NuFi');
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,5);
imagesc(abs(f_cmm10 - f_benchmark));
title('Error: CMM 20 vs NuFi');
xlabel('Velocity index');
ylabel('Position index');
colorbar;

subplot(2,3,6);
imagesc(abs(f_cmm10 - f_cmm1000));
title('Error: CMM 20 vs CMM 1000');
xlabel('Velocity index');
ylabel('Position index');
colorbar;

sgtitle(sprintf('Distribution Function Comparison (%s) at t ≈ %.1f', test_case_name, target_time));

% Save results
results = struct();
results.test_case = test_case_name;
results.benchmark_method = "NuFi";
results.target_time = target_time;
results.actual_time_benchmark = actual_time_benchmark;
results.actual_time_cmm1000 = actual_time_cmm1000;
results.actual_time_cmm10 = actual_time_cmm10;

% Error metrics
results.error_L2_abs_1000 = error_L2_abs_1000;
results.error_Linf_abs_1000 = error_Linf_abs_1000;
results.error_L2_rel_1000 = error_L2_rel_1000;
results.error_Linf_rel_1000 = error_Linf_rel_1000;
results.error_L2_abs_10 = error_L2_abs_10;
results.error_Linf_abs_10 = error_Linf_abs_10;
results.error_L2_rel_10 = error_L2_rel_10;
results.error_Linf_rel_10 = error_Linf_rel_10;

% Distribution functions
results.f_benchmark = f_benchmark;
results.f_cmm1000 = f_cmm1000;
results.f_cmm10 = f_cmm10;

% Map stack data for CMM N_remap = 20
results.map_stack_cmm10 = map_stack_cmm10;
results.Nmaps_cmm10 = Nmaps_cmm10;
results.grid_info_cmm10 = params_cmm10.grids(1); % Store grid info for plotting

results_filename = sprintf('remap_error_results_%s.mat', test_case_name);
save(results_filename, 'results');
fprintf('Results saved to %s\n', results_filename);

% Clean up temporary files (preserve benchmark data)
delete('temp_config.mat');
fprintf('Temporary files cleaned up. Benchmark data preserved in:\n');
fprintf('  NuFi: %s\n', benchmark_filename);
fprintf('  CMM 1000: %s\n', cmm1000_filename);