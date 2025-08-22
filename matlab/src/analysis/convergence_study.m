%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convergence Study Script for CMM Grid Resolution
% Studies how error decreases as grid resolution (Nx, Nv) increases
% Fixed N_remap = 10, Benchmark: NuFi method (most accurate reference)
% Test: CMM with various grid resolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%% Convergence Study Parameters
grid_resolutions = [64, 128, 256, 512, 1024];  % Grid sizes (Nx = Nv)
fixed_N_remap = 10;  % Fixed N_remap value
target_time = 10;    % Time at which to measure error
n_tests = length(grid_resolutions);

% Preallocate results arrays
error_L2_abs = zeros(n_tests, 1);
error_L1_abs = zeros(n_tests, 1);
error_Linf_abs = zeros(n_tests, 1);
error_L2_rel = zeros(n_tests, 1);
error_L1_rel = zeros(n_tests, 1);
error_Linf_rel = zeros(n_tests, 1);
cpu_times = zeros(n_tests, 1);
actual_times = zeros(n_tests, 1);

%% Step 1: Run CMM Benchmark (with finest grid)
fprintf('=== GRID CONVERGENCE STUDY ===\n');
fprintf('Running CMM benchmark simulation with finest grid...\n');
tic;
PARAMS_two_stream;
params.Nx = max(grid_resolutions);  % Use finest grid for benchmark
params.Nv = max(grid_resolutions);
params.N_remap = fixed_N_remap;
params.method = "CMM";
params.Tend = target_time;
[params_benchmark, data_benchmark] = Sim(params);
benchmark_time = toc;

% Extract benchmark data first
time_benchmark = data_benchmark.time;
[~, idx_benchmark] = min(abs(time_benchmark - target_time));
f_benchmark = data_benchmark.fs(:,:,idx_benchmark);
actual_time_benchmark = time_benchmark(idx_benchmark);

% Save benchmark and clear everything
save('benchmark_convergence.mat', 'params_benchmark', 'data_benchmark', 'benchmark_time', ...
     'f_benchmark', 'actual_time_benchmark');
clear all
addpath(genpath('./src/'),genpath('./params/'))

% Restore needed variables
grid_resolutions = [64, 128, 256, 512, 1024];
fixed_N_remap = 10;
target_time = 10;
n_tests = length(grid_resolutions);

% Load benchmark data
load('benchmark_convergence.mat');

fprintf('Benchmark completed in %.2f seconds (grid %dx%d) at time %.4f\n\n', ...
    benchmark_time, max(grid_resolutions), max(grid_resolutions), actual_time_benchmark);

% Initialize results file
results_temp = struct();
results_temp.grid_resolutions = grid_resolutions;
results_temp.fixed_N_remap = fixed_N_remap;
results_temp.target_time = target_time;
results_temp.n_tests = n_tests;
results_temp.error_L2_abs = NaN(n_tests, 1);
results_temp.error_L1_abs = NaN(n_tests, 1);
results_temp.error_Linf_abs = NaN(n_tests, 1);
results_temp.error_L2_rel = NaN(n_tests, 1);
results_temp.error_L1_rel = NaN(n_tests, 1);
results_temp.error_Linf_rel = NaN(n_tests, 1);
results_temp.cpu_times = NaN(n_tests, 1);
save('results_temp.mat', 'results_temp');

%% Step 2: Run CMM with different grid resolutions
fprintf('Running CMM grid convergence tests (N_remap = %d)...\n', fixed_N_remap);
for i = 1:n_tests
    grid_size = grid_resolutions(i);
    fprintf('Test %d/%d: Grid %dx%d... ', i, n_tests, grid_size, grid_size);
    
    % Save loop variables before clearing
    save('loop_vars.mat', 'i', 'grid_size');
    
    % Clear everything and reload parameters
    clear all
    addpath(genpath('./src/'),genpath('./params/'))
    
    % Restore needed variables
    load('results_temp.mat');
    load('benchmark_convergence.mat');
    load('loop_vars.mat');
    grid_resolutions = results_temp.grid_resolutions;
    fixed_N_remap = results_temp.fixed_N_remap;
    target_time = results_temp.target_time;
    n_tests = results_temp.n_tests;
    
    PARAMS_two_stream;
    
    % Set grid and CMM parameters
    params.Nx = grid_size;
    params.Nv = grid_size;
    params.N_remap = fixed_N_remap;
    params.method = "CMM";
    params.Tend = target_time;
    
    % Run simulation with timing
    tic;
    try
        [params_cmm, data_cmm] = Sim(params);
        cpu_time = toc;
        
        % Extract solution at target time
        time_cmm = data_cmm.time;
        [~, idx_cmm] = min(abs(time_cmm - target_time));
        f_cmm = data_cmm.fs(:,:,idx_cmm);
        actual_times(i) = time_cmm(idx_cmm);
        
        % Interpolate CMM solution to benchmark grid for comparison
        if grid_size ~= max(grid_resolutions)
            % Use the actual grids from simulation parameters
            x_cmm = linspace(0, params.Lx, grid_size);
            v_cmm = linspace(-params.Lv/2, params.Lv/2, grid_size);
            x_bench = linspace(0, params.Lx, max(grid_resolutions));  
            v_bench = linspace(-params.Lv/2, params.Lv/2, max(grid_resolutions));
            
            [X_cmm, V_cmm] = meshgrid(x_cmm, v_cmm);
            [X_bench, V_bench] = meshgrid(x_bench, v_bench);
            
            % Interpolate CMM solution to benchmark grid
            f_cmm_interp = interp2(X_cmm, V_cmm, f_cmm, X_bench, V_bench, 'cubic', 0);
            
            % Check for any NaN values from interpolation
            if any(isnan(f_cmm_interp(:)))
                fprintf('Warning: NaN values in interpolation, using nearest neighbor\n');
                f_cmm_interp = interp2(X_cmm, V_cmm, f_cmm, X_bench, V_bench, 'nearest', 0);
            end
        else
            f_cmm_interp = f_cmm;
        end
        
        % Debug: Check array sizes
        fprintf('f_cmm_interp size: %dx%d, f_benchmark size: %dx%d\n', ...
            size(f_cmm_interp,1), size(f_cmm_interp,2), size(f_benchmark,1), size(f_benchmark,2));
        
        % Ensure arrays have same size
        if ~isequal(size(f_cmm_interp), size(f_benchmark))
            error('Size mismatch: f_cmm_interp %s vs f_benchmark %s', ...
                mat2str(size(f_cmm_interp)), mat2str(size(f_benchmark)));
        end
        
        % Calculate errors
        % Absolute errors
        error_L2_abs = norm(f_cmm_interp - f_benchmark, 'fro');
        error_L1_abs = sum(abs(f_cmm_interp(:) - f_benchmark(:)));
        error_Linf_abs = max(abs(f_cmm_interp(:) - f_benchmark(:)));
        
        % Relative errors
        error_L2_rel = error_L2_abs / norm(f_benchmark, 'fro');
        error_L1_rel = error_L1_abs / sum(abs(f_benchmark(:)));
        error_Linf_rel = error_Linf_abs / max(abs(f_benchmark(:)));
        
        fprintf('completed in %.2f s, L2 error = %.2e\n', cpu_time, error_L2_rel);
        
        % Save results to temporary file
        results_temp.error_L2_abs(i) = error_L2_abs;
        results_temp.error_L1_abs(i) = error_L1_abs;
        results_temp.error_Linf_abs(i) = error_Linf_abs;
        results_temp.error_L2_rel(i) = error_L2_rel;
        results_temp.error_L1_rel(i) = error_L1_rel;  
        results_temp.error_Linf_rel(i) = error_Linf_rel;
        results_temp.cpu_times(i) = cpu_time;
        save('results_temp.mat', 'results_temp');
        
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        fprintf('  Error details: %s (line %d in %s)\n', ME.message, ME.stack(1).line, ME.stack(1).name);
        % Results already initialized as NaN in results_temp
    end
end

%% Step 3: Load final results and analyze convergence
load('results_temp.mat');
fprintf('\n=== CONVERGENCE ANALYSIS ===\n');

% Filter out failed cases
error_L2_abs = results_temp.error_L2_abs;
error_L1_abs = results_temp.error_L1_abs;
error_Linf_abs = results_temp.error_Linf_abs;
error_L2_rel = results_temp.error_L2_rel;
error_L1_rel = results_temp.error_L1_rel;
error_Linf_rel = results_temp.error_Linf_rel;
cpu_times = results_temp.cpu_times;
grid_resolutions = results_temp.grid_resolutions;

% Exclude benchmark grid (1024x1024) from convergence analysis since it has zero error
valid_idx = ~isnan(error_L2_rel) & error_L2_rel > 0;
grid_valid = grid_resolutions(valid_idx);
error_L2_valid = error_L2_rel(valid_idx);
error_L1_valid = error_L1_rel(valid_idx);
error_Linf_valid = error_Linf_rel(valid_idx);

if sum(valid_idx) >= 3  % Need at least 3 points for convergence rate
    % Calculate convergence rates using log-log regression
    % Error ∝ (1/h)^p where h = 1/grid_size → log(Error) = -p*log(grid_size) + C
    
    % L2 convergence rate
    p_L2 = polyfit(log(grid_valid), log(error_L2_valid), 1);
    convergence_rate_L2 = -p_L2(1);
    
    % L1 convergence rate  
    p_L1 = polyfit(log(grid_valid), log(error_L1_valid), 1);
    convergence_rate_L1 = -p_L1(1);
    
    % L∞ convergence rate
    p_Linf = polyfit(log(grid_valid), log(error_Linf_valid), 1);
    convergence_rate_Linf = -p_Linf(1);
    
    fprintf('Convergence rates (Error ∝ (grid_size)^(-p)):\n');
    fprintf('  L2 norm:  p = %.3f\n', convergence_rate_L2);
    fprintf('  L1 norm:  p = %.3f\n', convergence_rate_L1);  
    fprintf('  L∞ norm:  p = %.3f\n', convergence_rate_Linf);
    
    % Expected theoretical rates for different interpolation orders:
    fprintf('\nTheoretical convergence rates for spatial discretization:\n');
    fprintf('  2nd order methods: p ≈ 2\n');
    fprintf('  4th order methods: p ≈ 4\n');
    fprintf('  Spectral methods:  p > 6\n');
    
else
    fprintf('Not enough valid data points for convergence analysis\n');
    convergence_rate_L2 = NaN;
    convergence_rate_L1 = NaN;
    convergence_rate_Linf = NaN;
end

%% Step 4: Results Summary
fprintf('\n=== DETAILED RESULTS ===\n');
fprintf('Grid Size\tL2_rel\t\t\tL1_rel\t\t\tLinf_rel\t\tCPU_time\n');
fprintf('---------\t------\t\t\t------\t\t\t--------\t\t--------\n');
for i = 1:n_tests
    if valid_idx(i)
        fprintf('%dx%d\t\t%.3e\t\t%.3e\t\t%.3e\t\t%.2f s\n', ...
            grid_resolutions(i), grid_resolutions(i), error_L2_rel(i), error_L1_rel(i), error_Linf_rel(i), cpu_times(i));
    else
        fprintf('%dx%d\t\tFAILED\t\t\tFAILED\t\t\tFAILED\t\t\tFAILED\n', grid_resolutions(i), grid_resolutions(i));
    end
end

%% Step 5: Plotting
figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Convergence plot
subplot(2,2,1);
if sum(valid_idx) >= 2
    loglog(grid_valid, error_L2_valid, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    loglog(grid_valid, error_L1_valid, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    loglog(grid_valid, error_Linf_valid, 'g^-', 'LineWidth', 2, 'MarkerSize', 8);
    
    % Add theoretical slopes
    if ~isnan(convergence_rate_L2)
        grid_theory = [min(grid_valid), max(grid_valid)];
        error_theory = error_L2_valid(1) * (grid_theory/grid_valid(1)).^(-convergence_rate_L2);
        loglog(grid_theory, error_theory, 'b--', 'LineWidth', 1);
        text(grid_theory(end)/2, error_theory(end)*2, sprintf('slope = %.2f', convergence_rate_L2), 'Color', 'blue');
    end
    
    xlabel('Grid Size (Nx = Nv)');
    ylabel('Relative Error');
    title('CMM Grid Convergence (N_{remap} = 10)');
    legend('L2 norm', 'L1 norm', 'L∞ norm', 'Location', 'northeast');
    grid on;
else
    text(0.5, 0.5, 'Insufficient data for plotting', 'HorizontalAlignment', 'center');
end

% CPU time vs Grid Size
subplot(2,2,2);
if sum(valid_idx) >= 2
    loglog(grid_valid, cpu_times(valid_idx), 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Grid Size (Nx = Nv)');
    ylabel('CPU Time (seconds)');
    title('Computational Cost vs Grid Size');
    grid on;
end

% Error vs CPU time (efficiency plot)
subplot(2,2,3);
if sum(valid_idx) >= 2
    loglog(cpu_times(valid_idx), error_L2_valid, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('CPU Time (seconds)');
    ylabel('L2 Relative Error');
    title('Accuracy vs Computational Cost');
    grid on;
end

% Error distribution
subplot(2,2,4);
if sum(valid_idx) >= 2
    bar(1:sum(valid_idx), [error_L2_valid, error_L1_valid, error_Linf_valid]);
    set(gca, 'XTickLabel', arrayfun(@num2str, grid_valid, 'UniformOutput', false));
    set(gca, 'YScale', 'log');
    xlabel('Grid Size');
    ylabel('Relative Error');
    title('Error Comparison by Norm');
    legend('L2', 'L1', 'L∞', 'Location', 'northeast');
    grid on;
end

sgtitle(sprintf('CMM Grid Convergence Study (N_{remap} = %d) at t = %.1f', fixed_N_remap, target_time));

%% Step 6: Save Results
results = struct();
results.grid_resolutions = grid_resolutions;
results.fixed_N_remap = fixed_N_remap;
results.target_time = target_time;
results.benchmark_time = benchmark_time;
results.benchmark_grid_size = max(grid_resolutions);
results.error_L2_abs = error_L2_abs;
results.error_L1_abs = error_L1_abs;
results.error_Linf_abs = error_Linf_abs;
results.error_L2_rel = error_L2_rel;
results.error_L1_rel = error_L1_rel;
results.error_Linf_rel = error_Linf_rel;
results.cpu_times = cpu_times;
results.actual_times = actual_times;
results.convergence_rate_L2 = convergence_rate_L2;
results.convergence_rate_L1 = convergence_rate_L1;
results.convergence_rate_Linf = convergence_rate_Linf;

save('convergence_study_results.mat', 'results');
fprintf('\nResults saved to convergence_study_results.mat\n');

% Clean up temporary files
delete('benchmark_convergence.mat');
delete('results_temp.mat');
if exist('loop_vars.mat', 'file')
    delete('loop_vars.mat');
end
fprintf('Convergence study completed!\n');