%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison script for NuFi vs CMM methods
% This script runs the same plasma simulation using both methods and
% compares their accuracy and computational performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select test case and configure parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load base parameters (you can change this to test different cases)
PARAMS_two_stream;  % Options: PARAMS_landau_damping, PARAMS_two_stream, PARAMS_ion_acoustic_waves

% Override ending time for comparison
params.Tend = 50;

% Store original parameters for comparison
params_base = params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with NuFi method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation with NuFi method...\n');
params_nufi = params_base;
params_nufi.method = "NuFi";

tic()
[params_nufi_result, data_nufi] = Sim(params_nufi);
time_nufi = toc();
fprintf('NuFi simulation completed in %.2f seconds\n', time_nufi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run simulation with CMM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Running simulation with CMM method...\n');
params_cmm = params_base;
params_cmm.method = "CMM";

tic()
[params_cmm_result, data_cmm] = Sim(params_cmm);
time_cmm = toc();
fprintf('CMM simulation completed in %.2f seconds\n', time_cmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare results and calculate errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n=== PERFORMANCE COMPARISON ===\n');
fprintf('NuFi computation time: %.2f seconds\n', time_nufi);
fprintf('CMM computation time:  %.2f seconds\n', time_cmm);
fprintf('Speedup factor: %.2f\n', time_nufi/time_cmm);

% Extract final distribution functions
fs_nufi_final = data_nufi.fs(:,:,:,end);
fs_cmm_final = data_cmm.fs(:,:,:,end);

% Calculate L2 error in distribution function
error_fs_L2 = sqrt(sum((fs_nufi_final - fs_cmm_final).^2, 'all')) / sqrt(sum(fs_nufi_final.^2, 'all'));

% Calculate maximum absolute error
error_fs_max = max(abs(fs_nufi_final - fs_cmm_final), [], 'all');

% Compare electric field evolution - read from saved CSV files
% Load diagnostic data from CSV files for electric field energy comparison
try
    % Try to read from diagnostic files first
    species_name = params_nufi_result.species_name(1);
    filename_nufi = fullfile(params_nufi_result.data_dir, species_name+".csv");
    filename_cmm = fullfile(params_cmm_result.data_dir, species_name+".csv");
    
    if exist(filename_nufi, 'file') && exist(filename_cmm, 'file')
        table_nufi = readtable(filename_nufi);
        table_cmm = readtable(filename_cmm);
        E_nufi = table_nufi.Epot;
        E_cmm = table_cmm.Epot;
    else
        % Fallback: use electric field data if available
        if isfield(data_nufi, 'Efield') && isfield(data_cmm, 'Efield')
            % Calculate electric field energy from field data
            E_nufi = 0.5 * sum(data_nufi.Efield.^2, 1) * params_base.grids(1).dx;
            E_cmm = 0.5 * sum(data_cmm.Efield.^2, 1) * params_base.grids(1).dx;
        else
            warning('Electric field data not found, skipping E-field comparison');
            E_nufi = [];
            E_cmm = [];
        end
    end
catch
    warning('Could not load electric field data, skipping E-field comparison');
    E_nufi = [];
    E_cmm = [];
end

% Calculate electric field errors if data is available
if ~isempty(E_nufi) && ~isempty(E_cmm)
    % Ensure same length for comparison
    min_len = min(length(E_nufi), length(E_cmm));
    E_nufi = E_nufi(1:min_len);
    E_cmm = E_cmm(1:min_len);
    
    error_E_L2 = sqrt(mean((E_nufi - E_cmm).^2));
    error_E_max = max(abs(E_nufi - E_cmm));
else
    error_E_L2 = NaN;
    error_E_max = NaN;
end

fprintf('\n=== ERROR ANALYSIS ===\n');
fprintf('Distribution function L2 relative error: %.6e\n', error_fs_L2);
fprintf('Distribution function max absolute error: %.6e\n', error_fs_max);
if ~isnan(error_E_L2)
    fprintf('Electric field L2 error: %.6e\n', error_E_L2);
    fprintf('Electric field max error: %.6e\n', error_E_max);
else
    fprintf('Electric field comparison: Data not available\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create comparison plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: Electric field time evolution comparison
figure(1);
clf;
if ~isempty(E_nufi) && ~isempty(E_cmm)
    % Get time data from CSV if available
    try
        species_name = params_nufi_result.species_name(1);
        filename_nufi = fullfile(params_nufi_result.data_dir, species_name+".csv");
        if exist(filename_nufi, 'file')
            table_nufi = readtable(filename_nufi);
            t_save = table_nufi.time;
        else
            t_save = linspace(0, params_base.Tend, length(E_nufi));
        end
    catch
        t_save = linspace(0, params_base.Tend, length(E_nufi));
    end
    
    subplot(2,1,1);
    plot(t_save, E_nufi, 'b-', 'LineWidth', 2, 'DisplayName', 'NuFi');
    hold on;
    plot(t_save, E_cmm, 'r--', 'LineWidth', 2, 'DisplayName', 'CMM');
    xlabel('Time');
    ylabel('Electric Field Energy');
    title('Electric Field Energy Evolution');
    legend('Location', 'best');
    grid on;

    subplot(2,1,2);
    plot(t_save, abs(E_nufi - E_cmm), 'k-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('|E_{NuFi} - E_{CMM}|');
    title('Absolute Error in Electric Field Energy');
    grid on;
    set(gca, 'YScale', 'log');
else
    text(0.5, 0.5, 'Electric field data not available for comparison', ...
         'HorizontalAlignment', 'center', 'FontSize', 14);
    title('Electric Field Comparison - Data Not Available');
end

% Figure 2: Final distribution function comparison
figure(2);
clf;
if params_base.Ns >= 1
    subplot(2,2,1);
    pcolor(params_nufi_result.grids(1).X, params_nufi_result.grids(1).V, fs_nufi_final(:,:,1));
    shading flat;
    colorbar;
    title('NuFi Final Distribution');
    xlabel('Position x');
    ylabel('Velocity v');
    
    subplot(2,2,2);
    pcolor(params_cmm_result.grids(1).X, params_cmm_result.grids(1).V, fs_cmm_final(:,:,1));
    shading flat;
    colorbar;
    title('CMM Final Distribution');
    xlabel('Position x');
    ylabel('Velocity v');
    
    subplot(2,2,3);
    pcolor(params_nufi_result.grids(1).X, params_nufi_result.grids(1).V, ...
           fs_nufi_final(:,:,1) - fs_cmm_final(:,:,1));
    shading flat;
    colorbar;
    title('Difference (NuFi - CMM)');
    xlabel('Position x');
    ylabel('Velocity v');
    
    subplot(2,2,4);
    pcolor(params_nufi_result.grids(1).X, params_nufi_result.grids(1).V, ...
           abs(fs_nufi_final(:,:,1) - fs_cmm_final(:,:,1)));
    shading flat;
    colorbar;
    title('Absolute Difference |NuFi - CMM|');
    xlabel('Position x');
    ylabel('Velocity v');
end



% Figure 3: Performance summary
figure(3);
clf;
methods = {'NuFi', 'CMM'};
times = [time_nufi, time_cmm];
bar(times);
set(gca, 'XTickLabel', methods);
ylabel('Computation Time (seconds)');
title('Performance Comparison');
for i = 1:length(times)
    text(i, times(i) + 0.05*max(times), sprintf('%.2fs', times(i)), ...
         'HorizontalAlignment', 'center');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save comparison results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comparison_results = struct();
comparison_results.time_nufi = time_nufi;
comparison_results.time_cmm = time_cmm;
comparison_results.speedup = time_nufi/time_cmm;
comparison_results.error_fs_L2 = error_fs_L2;
comparison_results.error_fs_max = error_fs_max;
comparison_results.error_E_L2 = error_E_L2;
comparison_results.error_E_max = error_E_max;
comparison_results.params_used = params_base;

% Save to file
save('nufi_cmm_comparison_results.mat', 'comparison_results', 'data_nufi', 'data_cmm', ...
     'params_nufi_result', 'params_cmm_result');

fprintf('\nComparison results saved to: nufi_cmm_comparison_results.mat\n');
fprintf('Plots have been generated in figures 1-4\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n');
fprintf('=========================================\n');
fprintf('           SUMMARY REPORT\n');
fprintf('=========================================\n');
fprintf('Test case: %s\n', params_base.mycase);
fprintf('Grid size: %dx%d (Nx x Nv)\n', params_base.Nx, params_base.Nv);
fprintf('Time steps: %d\n', params_base.Nt_max);
fprintf('Final time: %.2f\n', params_base.Tend);
fprintf('-----------------------------------------\n');
fprintf('NuFi time:     %.2f seconds\n', time_nufi);
fprintf('CMM time:      %.2f seconds\n', time_cmm);
fprintf('Speedup:       %.2fx\n', time_nufi/time_cmm);
fprintf('-----------------------------------------\n');
fprintf('Distribution function L2 error: %.2e\n', error_fs_L2);
if ~isnan(error_E_L2)
    fprintf('Electric field L2 error:       %.2e\n', error_E_L2);
else
    fprintf('Electric field L2 error:       N/A\n');
end
fprintf('=========================================\n'); 