%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplified CMM and NuFi method error comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select test case
PARAMS_two_stream;
%PARAMS_landau_damping;
%PARAMS_ion_acoustic_waves;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run CMM method
fprintf("Running CMM method...\n")
tic()
params.method = "CMM";
[params_cmm, data_cmm] = Sim(params);
t_cmm = toc();
fprintf("CMM method runtime: %2.2f seconds\n", t_cmm)

%% Run NuFi method
fprintf("Running NuFi method...\n")
tic()
params.method = "NuFi";
[params_nufi, data_nufi] = Sim(params);
t_nufi = toc();
fprintf("NuFi method runtime: %2.2f seconds\n", t_nufi)

%% Calculate errors
fprintf("Calculating errors...\n")

% Read diagnostic data
species_name = params_cmm.species_name(1);
cmm_file = fullfile(params_cmm.data_dir, species_name + ".csv");
nufi_file = fullfile(params_nufi.data_dir, species_name + ".csv");

if exist(cmm_file, 'file') && exist(nufi_file, 'file')
    cmm_data = readtable(cmm_file);
    nufi_data = readtable(nufi_file);
    
    % Ensure data length matches
    min_length = min(height(cmm_data), height(nufi_data));
    cmm_data = cmm_data(1:min_length, :);
    nufi_data = nufi_data(1:min_length, :);
    
    % Calculate errors for each diagnostic quantity
    fields = {'Mass', 'Momentum', 'Epot', 'Ekin', 'Etot', 'L2norm'};
    
    error_stats = struct();
    
    for i = 1:length(fields)
        field_name = fields{i};
        
        if ismember(field_name, cmm_data.Properties.VariableNames) && ...
           ismember(field_name, nufi_data.Properties.VariableNames)
            
            cmm_values = cmm_data.(field_name);
            nufi_values = nufi_data.(field_name);
            
            % Absolute error
            absolute_error = abs(cmm_values - nufi_values);
            relative_error = abs(cmm_values - nufi_values) ./ abs(cmm_values);
            
            error_stats.(field_name).absolute_error = absolute_error;
            error_stats.(field_name).relative_error = relative_error;
            error_stats.(field_name).max_absolute = max(absolute_error);
            error_stats.(field_name).mean_absolute = mean(absolute_error);
            error_stats.(field_name).max_relative = max(relative_error);
            error_stats.(field_name).mean_relative = mean(relative_error);
        end
    end
    
    % Plot error analysis
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Total energy error
    subplot(2, 3, 1);
    if isfield(error_stats, 'Etot')
        plot(cmm_data.time, error_stats.Etot.absolute_error, 'b-', 'LineWidth', 2);
        title('Total Energy Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    % Subplot 2: Kinetic energy error
    subplot(2, 3, 2);
    if isfield(error_stats, 'Ekin')
        plot(cmm_data.time, error_stats.Ekin.absolute_error, 'r-', 'LineWidth', 2);
        title('Kinetic Energy Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    % Subplot 3: Potential energy error
    subplot(2, 3, 3);
    if isfield(error_stats, 'Epot')
        plot(cmm_data.time, error_stats.Epot.absolute_error, 'g-', 'LineWidth', 2);
        title('Potential Energy Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    % Subplot 4: Mass error
    subplot(2, 3, 4);
    if isfield(error_stats, 'Mass')
        plot(cmm_data.time, error_stats.Mass.absolute_error, 'm-', 'LineWidth', 2);
        title('Mass Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    % Subplot 5: Momentum error
    subplot(2, 3, 5);
    if isfield(error_stats, 'Momentum')
        plot(cmm_data.time, error_stats.Momentum.absolute_error, 'c-', 'LineWidth', 2);
        title('Momentum Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    % Subplot 6: L2 norm error
    subplot(2, 3, 6);
    if isfield(error_stats, 'L2norm')
        plot(cmm_data.time, error_stats.L2norm.absolute_error, 'k-', 'LineWidth', 2);
        title('L2 Norm Absolute Error');
        xlabel('Time');
        ylabel('Absolute Error');
        grid on;
    end
    
    sgtitle('CMM vs NuFi Error Analysis', 'FontSize', 16);
    
    % Print statistics
    fprintf('\n=== Error Statistics ===\n');
    for i = 1:length(fields)
        field_name = fields{i};
        if isfield(error_stats, field_name)
            fprintf('%s:\n', field_name);
            fprintf('  Max absolute error: %.2e\n', error_stats.(field_name).max_absolute);
            fprintf('  Mean absolute error: %.2e\n', error_stats.(field_name).mean_absolute);
            fprintf('  Max relative error: %.2e\n', error_stats.(field_name).max_relative);
            fprintf('  Mean relative error: %.2e\n', error_stats.(field_name).mean_relative);
            fprintf('\n');
        end
    end
    
    % Save results
    result_dir = sprintf('./data/simple_error_comparison_%s/', params_cmm.mycase);
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    
    % Save error data
    save(fullfile(result_dir, 'error_stats.mat'), 'error_stats', 'params_cmm', 'params_nufi');
    
    % Create error data table
    error_table = table();
    error_table.time = cmm_data.time;
    
    for i = 1:length(fields)
        field_name = fields{i};
        if isfield(error_stats, field_name)
            error_table.(sprintf('%s_absolute_error', field_name)) = error_stats.(field_name).absolute_error;
            error_table.(sprintf('%s_relative_error', field_name)) = error_stats.(field_name).relative_error;
        end
    end
    
    writetable(error_table, fullfile(result_dir, 'error_data.csv'));
    
    fprintf('Error analysis results saved to: %s\n', result_dir);
    
else
    fprintf('Error: Cannot find diagnostic data files\n');
    fprintf('CMM file: %s\n', cmm_file);
    fprintf('NuFi file: %s\n', nufi_file);
end

fprintf("Error analysis completed!\n") 