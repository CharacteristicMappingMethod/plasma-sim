%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conservation Quantities Analysis for Three Methods - Multiple Test Cases
%
% This script analyzes the conservation of mass, momentum, and energy
% for four different numerical methods on selected test cases:
% 1. CMM with N_remap = 20 (Nmap = Nsample)
% 2. CMM_vargrid with N_remap = 20 (Nmap from params file)
% 3. NuFi 
% 4. Predictor-corrector scheme
%
% Available Test Cases:
% - Landau Damping (alpha = 0.05)
% - Two-Stream Instability
%
% Note: The script reads conservation quantities directly from CSV files
% created by the measure.m function during simulation. No .mat files are
% saved or loaded - all data is read directly from CSV files.
%
% Physical Quantities:
%   - Total Mass M(t) = ∬f(x,v,t)dxdv
%   - Total Momentum P(t) = ∬vf(x,v,t)dxdv  
%   - Total Energy E(t) = E_kin(t) + E_pot(t)
%     - E_kin(t) = 1/2 ∬f(x,v,t)|v|²dxdv
%     - E_pot(t) = 1/2 ∫|E(x,t)|²dx
%   - L2 Norm ||f||_2(t) = (∬|f(x,v,t)|²dxdv)^(1/2)
%
% Error Calculations:
%   - Mass Error (Relative): |M(t) - M(0)| / |M(0)|
%   - Momentum Error (Absolute): |P(t) - P(0)|
%   - Energy Error (Relative): |E(t) - E(0)| / |E(0)|
%   - L2 Norm Error (Relative): ||f||_2(t) - ||f||_2(0)| / ||f||_2(0)|
%
% Methods: CMM, CMM_vargrid, NuFi, Predictor-Corrector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS

%% Test Case Selection
fprintf('=== Conservation Quantities Analysis - Three Methods Comparison ===\n');
fprintf('Available test cases:\n');
fprintf('  1. Landau Damping (alpha = 0.05)\n');
fprintf('  2. Two-Stream Instability\n');

% User input for test case selection
while true
    test_case_choice = input('Select test case (1 or 2): ');
    if test_case_choice == 1 || test_case_choice == 2
        break;
    else
        fprintf('Invalid choice. Please enter 1 or 2.\n');
    end
end

% Configure based on test case choice
if test_case_choice == 1
    % Landau Damping configuration
    test_case_name = 'landau_damping';
    test_case_display = 'Landau Damping (α = 0.05)';
    params_script = 'PARAMS_landau_damping';
    target_time = 50.0;
    plot_filename = sprintf('conservation_three_methods_landau_damping_Tend%.0f', target_time);
else
    % Two-Stream Instability configuration
    test_case_name = 'two_stream';
    test_case_display = 'Two-Stream Instability';
    params_script = 'PARAMS_two_stream';
    target_time = 50.0;  % Two-stream typically needs longer time
    plot_filename = sprintf('conservation_three_methods_two_stream_Tend%.0f', target_time);
end

%% Configuration
fprintf('\nSelected test case: %s\n', test_case_display);
fprintf('Methods: CMM, CMM-vargrid, NuFi, Predictor-Corrector\n');
fprintf('Time range: 0 to %.1f\n\n', target_time);

% Analysis parameters
methods = {"CMM", "CMM_vargrid", "NuFi", "predcorr"};
method_names = {"CMM-NuFi", "CMM-vargrid", "NuFi", "Pred.-Corr."};
colors = {'red', 'magenta', 'blue', 'green'};

% Data point spacing for plotting (downsample data points)
delta_idx = 2;

%% Read data directly from CSV files
results = cell(4, 1);

for method_idx = 1:4
    method = methods{method_idx};
    method_name = method_names{method_idx};
    
    fprintf('Reading data for %s (%s)...\n', method_name, test_case_display);
    
    % Construct data directory path based on test case and method
    % Pattern: ./data/{test_case_name}_Tend{Tend}_{method}/
    data_dir = sprintf('./data/%s_Tend%.0f_%s/', test_case_name, target_time, method);
    data_dir_standard = data_dir;  % Keep standard path for new simulations
    
    % Read CSV file
    species_name = "electrons";  % Default species name for single-species cases
    csv_file = fullfile(data_dir, species_name + ".csv");
    
    if exist(csv_file, 'file')
        fprintf('  Found CSV file: %s\n', csv_file);
        
        % Load parameters to compute initial conditions (t=0)
        % We need params to compute initial values but don't run simulation
        eval(params_script);
        params.method = method;
        params.Tend = target_time;
        if strcmp(method, "CMM") || strcmp(method, "CMM_vargrid")
            params.N_remap = 20;
        end
        
        % Initialize simulation parameters (but don't run simulation)
        % We only need grids and initial distribution for t=0 calculation
        [sim_params, ~] = initialize_simulation(params);
        
        % Read conservation quantities from CSV
        method_results = read_conservation_quantities_from_csv(csv_file, sim_params);
        method_results.method = method;
        method_results.method_name = method_name;
        method_results.test_case = test_case_name;
        method_results.test_case_display = test_case_display;
        method_results.csv_file = csv_file;
        
        results{method_idx} = method_results;
        fprintf('  Successfully read %d time points from CSV\n', length(method_results.time_array));
    else
        fprintf('  CSV file not found: %s\n', csv_file);
        fprintf('  Running simulation for %s (%s)...\n', method_name, test_case_display);
        
        % Load appropriate parameters based on test case
        eval(params_script);
        
        % Set method-specific parameters
        params.method = method;
        params.Tend = target_time;
        
        if strcmp(method, "CMM")
            params.N_remap = 20;
            % For CMM method, set Nmap = Nsample
            params.Nmap = params.Nsample;
            fprintf('  Setting Nmap = Nsample for CMM method: Nmap = [%d, %d]\n', params.Nmap(1), params.Nmap(2));
        elseif strcmp(method, "CMM_vargrid")
            params.N_remap = 20;
            % For CMM_vargrid, use Nmap from params file (already set)
            fprintf('  Using Nmap from params file for CMM_vargrid: Nmap = [%d, %d]\n', params.Nmap(1), params.Nmap(2));
        end
        
        % Set data directory to ensure CSV is saved in the expected location
        % Use standard directory path for new simulations (not CMM_vargrid)
        params.data_dir = data_dir_standard;
        
        fprintf('  Running %s simulation (%s, Tend=%.1f)...\n', method_name, test_case_display, params.Tend);
        
        % Run simulation (this will create the CSV file)
        tic_sim = tic();
        [sim_params, sim_data] = Sim(params);
        sim_time = toc(tic_sim);
        
        fprintf('  Simulation completed in %.1f seconds\n', sim_time);
        
        % Update CSV file path to standard directory (where simulation saved it)
        csv_file = fullfile(data_dir_standard, species_name + ".csv");
        
        % Verify CSV file was created
        if exist(csv_file, 'file')
            fprintf('  Reading conservation quantities from newly created CSV file...\n');
            
            % Read conservation quantities from CSV
            method_results = read_conservation_quantities_from_csv(csv_file, sim_params);
            method_results.method = method;
            method_results.method_name = method_name;
            method_results.test_case = test_case_name;
            method_results.test_case_display = test_case_display;
            method_results.csv_file = csv_file;
            
            results{method_idx} = method_results;
            fprintf('  Successfully read %d time points from CSV\n', length(method_results.time_array));
        else
            fprintf('  Error: CSV file was not created after simulation: %s\n', csv_file);
            fprintf('  Skipping %s method\n', method_name);
            results{method_idx} = [];
        end
    end
end

%% Create Publication-Quality Conservation Quantities Plot
fprintf('\nCreating publication-quality conservation quantities plot for %s...\n', test_case_display);

% Publication-quality figure settings
fig = figure('Position', [100, 100, 1000, 1000], 'Units', 'pixels', ...
             'PaperPositionMode', 'auto', 'Color', 'white', 'Visible', 'on');
% Set paper size in inches for publication
set(fig, 'PaperUnits', 'inches', 'PaperSize', [10, 10]);

% Professional color scheme - colorblind-friendly palette
colors = [0.8500, 0.3250, 0.0980;    % Red-orange (CMM)
          0.9290, 0.6940, 0.1250;    % Yellow-orange (CMM_vargrid)
          0.0000, 0.4470, 0.7410;    % Blue (NuFi)
          0.4940, 0.1840, 0.5560];   % Purple (Pred-corr)
          
line_styles = {'-', '--', ':', '-.'};  % Distinct line styles (solid, dashed, dotted, dash-dot)
line_widths = [1.5, 1.5, 1.5, 1.5];
marker_styles = {'none', 'x', '+', 'o'};

% Set publication font properties with larger sizes
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

% Top subplot: Mass Error
subplot(4, 1, 1);
hold on;
h_mass = [];
all_mass_errors = [];
for method_idx = 1:4
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        mass_error = results{method_idx}.relative_mass_error;
        
        % Downsample data points using delta_idx
        idx_plot = 1:delta_idx:length(time_array);
        if idx_plot(end) ~= length(time_array)
            idx_plot = [idx_plot, length(time_array)];  % Always include last point
        end
        time_array_plot = time_array(idx_plot);
        mass_error_plot = mass_error(idx_plot);
        
        % Handle zero/negative values for log scale by setting minimum threshold
        zero_mask = mass_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(mass_error_plot(mass_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            mass_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array_plot, mass_error_plot, ...
                 'Color', colors(method_idx, :), ...
                 'LineStyle', line_styles{method_idx}, ...
                 'LineWidth', line_widths(method_idx), ...
                 'Marker', marker_styles{method_idx}, ...
                 'DisplayName', method_names{method_idx});
        h_mass = [h_mass, h];
        all_mass_errors = [all_mass_errors; mass_error_plot(mass_error_plot > 0)];  % Collect non-zero errors
    end
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 14);
xlim([0, target_time]);
% Adaptive y-limits for mass error
if ~isempty(all_mass_errors)
    min_err = min(all_mass_errors);
    max_err = max(all_mass_errors);
    margin = 2;  % Factor for margin
    ylim([min_err/margin, max_err*margin]);
else
    ylim([1e-16, 1e-2]);  % Fallback if no data
end
ylabel('$|\Delta M|/M_0$', 'Interpreter', 'latex', 'FontSize', 16);
title('(a) Mass', 'FontSize', 16, 'FontWeight', 'normal');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
% No legend for mass subplot - will be placed on energy subplot

% Middle subplot: Momentum Error (Absolute)
subplot(4, 1, 2);
hold on;
h_momentum = [];
all_momentum_errors = [];
for method_idx = 1:4
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        % Handle backward compatibility for field names
        if isfield(results{method_idx}, 'absolute_momentum_error')
            momentum_error = results{method_idx}.absolute_momentum_error;
        else
            % Calculate from raw data if field doesn't exist
            P0 = results{method_idx}.total_momentum(1);
            momentum_error = abs(results{method_idx}.total_momentum - P0);
        end
        
        % Downsample data points using delta_idx
        idx_plot = 1:delta_idx:length(time_array);
        if idx_plot(end) ~= length(time_array)
            idx_plot = [idx_plot, length(time_array)];  % Always include last point
        end
        time_array_plot = time_array(idx_plot);
        momentum_error_plot = momentum_error(idx_plot);
        
        % Handle zero/negative values for log scale
        zero_mask = momentum_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(momentum_error_plot(momentum_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            momentum_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array_plot, momentum_error_plot, ...
                 'Color', colors(method_idx, :), ...
                 'LineStyle', line_styles{method_idx}, ...
                 'LineWidth', line_widths(method_idx), ...
                 'Marker', marker_styles{method_idx}, ...
                 'DisplayName', method_names{method_idx});
        h_momentum = [h_momentum, h];
        all_momentum_errors = [all_momentum_errors; momentum_error_plot(momentum_error_plot > 0)];  % Collect non-zero errors
    end
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 14);
xlim([0, target_time]);
% Adaptive y-limits for momentum error
if ~isempty(all_momentum_errors)
    min_err = min(all_momentum_errors);
    max_err = max(all_momentum_errors);
    margin = 2;  % Factor for margin
    ylim([min_err/margin, max_err*margin]);
else
    ylim([1e-16, 1e-2]);  % Fallback if no data
end
ylabel('$|\Delta P|$', 'Interpreter', 'latex', 'FontSize', 16);
title('(b) Momentum', 'FontSize', 16, 'FontWeight', 'normal');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
% No legend for momentum subplot

% Third subplot: Energy Error
subplot(4, 1, 3);
hold on;
h_energy = [];
all_energy_errors = [];
for method_idx = 1:4
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        energy_error = results{method_idx}.relative_energy_error;
        
        % Downsample data points using delta_idx
        idx_plot = 1:delta_idx:length(time_array);
        if idx_plot(end) ~= length(time_array)
            idx_plot = [idx_plot, length(time_array)];  % Always include last point
        end
        time_array_plot = time_array(idx_plot);
        energy_error_plot = energy_error(idx_plot);
        
        % Handle zero/negative values for log scale
        zero_mask = energy_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(energy_error_plot(energy_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            energy_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array_plot, energy_error_plot, ...
                 'Color', colors(method_idx, :), ...
                 'LineStyle', line_styles{method_idx}, ...
                 'LineWidth', line_widths(method_idx), ...
                 'Marker', marker_styles{method_idx}, ...
                 'DisplayName', method_names{method_idx});
        h_energy = [h_energy, h];
        all_energy_errors = [all_energy_errors; energy_error_plot(energy_error_plot > 0)];  % Collect non-zero errors
    end
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 14);
xlim([0, target_time]);
% Adaptive y-limits for energy error
if ~isempty(all_energy_errors)
    min_err = min(all_energy_errors);
    max_err = max(all_energy_errors);
    margin = 2;  % Factor for margin
    ylim([min_err/margin, max_err*margin]);
else
    ylim([1e-16, 1e-2]);  % Fallback if no data
end
ylabel('$|\Delta E|/E_0$', 'Interpreter', 'latex', 'FontSize', 16);
title('(c) Energy', 'FontSize', 16, 'FontWeight', 'normal');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
% No legend for energy subplot - will be placed on L2norm subplot

% Fourth subplot: L2norm Error
subplot(4, 1, 4);
hold on;
h_l2norm = [];
all_l2norm_errors = [];
for method_idx = 1:4
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        l2norm_error = results{method_idx}.relative_l2norm_error;
        
        % Downsample data points using delta_idx
        idx_plot = 1:delta_idx:length(time_array);
        if idx_plot(end) ~= length(time_array)
            idx_plot = [idx_plot, length(time_array)];  % Always include last point
        end
        time_array_plot = time_array(idx_plot);
        l2norm_error_plot = l2norm_error(idx_plot);
        
        % Handle zero/negative values for log scale
        zero_mask = l2norm_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(l2norm_error_plot(l2norm_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            l2norm_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array_plot, l2norm_error_plot, ...
                 'Color', colors(method_idx, :), ...
                 'LineStyle', line_styles{method_idx}, ...
                 'LineWidth', line_widths(method_idx), ...
                 'Marker', marker_styles{method_idx}, ...
                 'DisplayName', method_names{method_idx});
        h_l2norm = [h_l2norm, h];
        all_l2norm_errors = [all_l2norm_errors; l2norm_error_plot(l2norm_error_plot > 0)];  % Collect non-zero errors
    end
end
hold off;
set(gca, 'YScale', 'log', 'FontSize', 14);
xlim([0, target_time]);
% Adaptive y-limits for L2norm error
if ~isempty(all_l2norm_errors)
    min_err = min(all_l2norm_errors);
    max_err = max(all_l2norm_errors);
    margin = 2;  % Factor for margin
    ylim([min_err/margin, max_err*margin]);
else
    ylim([1e-16, 1e-2]);  % Fallback if no data
end
ylabel('$|\Delta \|f\|_2|/\|f_0\|_2$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Time t', 'FontSize', 16);
title('(d) L2 Norm', 'FontSize', 16, 'FontWeight', 'normal');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
legend(h_l2norm, 'Location', 'southeast', 'FontSize', 14, 'Box', 'off');

% Fine-tune subplot positioning for very tight layout
subplot_positions = [0.13, 0.76, 0.82, 0.18;   % Top subplot (Mass)
                     0.13, 0.55, 0.82, 0.18;   % Second subplot (Momentum)
                     0.13, 0.34, 0.82, 0.18;   % Third subplot (Energy)
                     0.13, 0.07, 0.82, 0.18]; % Bottom subplot (L2norm)

for i = 1:4
    subplot(4, 1, i);
    set(gca, 'Position', subplot_positions(i, :));
end

% Set paper properties for publication
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'inches', ...
    'PaperSize', [pos(3), pos(4)]);

% Display the plot
fprintf('Attempting to display figure...\n');
figure(fig);  % Bring figure to front
set(fig, 'Visible', 'on');  % Ensure figure is visible
drawnow;      % Force immediate display

% Debug information
fprintf('Figure handle: %s\n', class(fig));
fprintf('Figure visibility: %s\n', get(fig, 'Visible'));
fprintf('Figure position: [%.0f, %.0f, %.0f, %.0f]\n', get(fig, 'Position'));
fprintf('Figure should now be displayed on screen.\n');

% Pause briefly to ensure display
pause(0.5);

% Save publication-quality plots using save_fig_tikz
images_dir = './analysis/images';

% Create analysis/images directory if it doesn't exist
if ~exist(images_dir, 'dir')
    mkdir(images_dir);
    fprintf('Created %s directory\n', images_dir);
end

% Save using save_fig_tikz (generates PNG and TEX files)
plot_filename_pub = fullfile(images_dir, [plot_filename '_publication']);
if exist('save_fig_tikz', 'file')
    save_fig_tikz(plot_filename_pub, gcf);
    fprintf('Publication-quality plots saved using save_fig_tikz:\n');
    fprintf('  TEX (TikZ): %s.tex\n', plot_filename_pub);
    fprintf('  PNG (raster, 600 DPI): %s.png\n', plot_filename_pub);
else
    fprintf('Warning: save_fig_tikz function not found. Using standard MATLAB print.\n');
    % Fallback to standard saving
    print(gcf, [plot_filename_pub '.png'], '-dpng', '-r600');
    print(gcf, [plot_filename_pub '.eps'], '-depsc2', '-r600');
    fprintf('  PNG (raster): %s.png\n', plot_filename_pub);
    fprintf('  EPS (vector): %s.eps\n', plot_filename_pub);
end

%% Publication Summary
fprintf('\n=== Publication-Ready Analysis Summary ===\n');
fprintf('Test case: %s\n', test_case_display);
fprintf('Time range: 0 to %.1f plasma periods\n', target_time);
fprintf('Methods compared: %d\n', length(methods));
fprintf('Data point spacing (delta_idx): %d\n', delta_idx);
if exist('save_fig_tikz', 'file')
    fprintf('Publication formats: TEX (TikZ), PNG (600 DPI)\n');
else
    fprintf('Publication formats: PNG (600 DPI), EPS\n');
end

for method_idx = 1:4
    if ~isempty(results{method_idx})
        method_results = results{method_idx};
        fprintf('\n%s:\n', method_names{method_idx});
        fprintf('  Final time: %.2f\n', method_results.time_array(end));
        fprintf('  Time steps: %d\n', length(method_results.time_array));
        fprintf('  Final Mass Error (rel): %.2e\n', method_results.relative_mass_error(end));
        
        % Handle backward compatibility for momentum error field
        if isfield(method_results, 'absolute_momentum_error')
            final_momentum_error = method_results.absolute_momentum_error(end);
        else
            P0 = method_results.total_momentum(1);
            final_momentum_error = abs(method_results.total_momentum(end) - P0);
        end
        fprintf('  Final Momentum Error (abs): %.2e\n', final_momentum_error);
        
        fprintf('  Final Energy Error (rel): %.2e\n', method_results.relative_energy_error(end));
        
        % Handle backward compatibility for L2norm error field
        if isfield(method_results, 'relative_l2norm_error')
            final_l2norm_error = method_results.relative_l2norm_error(end);
        else
            L2norm0 = method_results.l2norm(1);
            final_l2norm_error = abs(method_results.l2norm(end) - L2norm0) / abs(L2norm0);
        end
        fprintf('  Final L2norm Error (rel): %.2e\n', final_l2norm_error);
        if isfield(method_results, 'csv_file')
            fprintf('  Data source: %s\n', method_results.csv_file);
        end
    else
        fprintf('\n%s: No data available\n', method_names{method_idx});
    end
end

fprintf('\nPublication-ready analysis complete!\n');
fprintf('All plots saved in multiple formats suitable for journal submission.\n');

%% Helper function to read conservation quantities from CSV file
function method_results = read_conservation_quantities_from_csv(csv_file, sim_params)
    fprintf('    Reading conservation quantities from CSV file...\n');
    fprintf('    CSV file: %s\n', csv_file);
    
    % Read CSV file
    data_table = readtable(csv_file);
    
    % Extract quantities from CSV (these start from t = dt, not t = 0)
    time_array_csv = data_table.time;
    total_mass_csv = data_table.Mass;
    total_momentum_csv = data_table.Momentum;
    kinetic_energy_csv = data_table.Ekin;
    potential_energy_csv = data_table.Epot;
    total_energy_csv = data_table.Etot;
    l2norm_csv = data_table.L2norm;
    
    % Calculate initial conditions (t=0) from initial distribution
    % This is missing from the CSV because measure() is called after time stepping
    fprintf('    Computing initial conditions (t=0)...\n');
    
    % Get initial distribution function and grid
    grid = sim_params.grids(1);
    fini = sim_params.fini{1};
    f0 = fini(grid.X, grid.V);
    
    % Calculate initial quantities
    M0 = sum(f0, "all") * grid.dx * grid.dv;
    P0 = sum(f0 .* grid.V, "all") * grid.dx * grid.dv;
    Ekin0 = 0.5 * sum(f0 .* (grid.V.^2), "all") * grid.dx * grid.dv;
    L2norm0 = sum((abs(f0).^2) .* grid.dv(:), "all") * grid.dx;
    
    % Initial electric field (from Poisson equation with initial distribution)
    % vPoisson expects fs as 3D array, so add species dimension
    f0_3d = zeros(size(f0, 1), size(f0, 2), 1);
    f0_3d(:, :, 1) = f0;
    Efield0 = vPoisson(f0_3d, sim_params.grids, sim_params.charge);
    Epot0 = 0.5 * sum(Efield0.^2) * grid.dx;
    
    % Use simulation's first mass value as reference for mass-conserving methods
    % This fixes the horizontal line issue for NuFi and predictor-corrector
    M0_sim = total_mass_csv(1);  % First simulation mass value
    P0_sim = total_momentum_csv(1);  % First simulation momentum value
    E0_sim = total_energy_csv(1);  % First simulation energy value
    L2norm0_sim = l2norm_csv(1);  % First simulation L2norm value
    
    % Combine initial conditions with CSV data using simulation reference values
    time_array = [0; time_array_csv];
    total_mass = [M0_sim; total_mass_csv];  % Use simulation M0
    total_momentum = [P0_sim; total_momentum_csv];  % Use simulation P0  
    total_energy = [E0_sim; total_energy_csv];  % Use simulation E0
    l2norm = [L2norm0_sim; l2norm_csv];  % Use simulation L2norm0
    
    % Also store kinetic and potential energy (use computed t=0 values)
    kinetic_energy = [Ekin0; kinetic_energy_csv];
    potential_energy = [Epot0; potential_energy_csv];
    
    fprintf('    Successfully read %d time points from CSV + 1 initial condition\n', length(time_array_csv));
    fprintf('    Total time points: %d (including t=0)\n', length(time_array));
    fprintf('    Using simulation reference values: M0=%.6e, P0=%.6e, E0=%.6e\n', M0_sim, P0_sim, E0_sim);
    
    % Calculate initial values for error computation
    M0 = total_mass(1);
    P0 = total_momentum(1);
    E0 = total_energy(1);
    L2norm0 = l2norm(1);
    
    % Calculate relative errors (except momentum which is absolute)
    relative_mass_error = abs(total_mass - M0) / abs(M0);
    absolute_momentum_error = abs(total_momentum - P0);  % Absolute error for momentum
    relative_energy_error = abs(total_energy - E0) / abs(E0);
    relative_l2norm_error = abs(l2norm - L2norm0) / abs(L2norm0);
    
    % Store results
    method_results = struct();
    method_results.time_array = time_array;
    method_results.total_mass = total_mass;
    method_results.total_momentum = total_momentum;
    method_results.kinetic_energy = kinetic_energy;
    method_results.potential_energy = potential_energy;
    method_results.total_energy = total_energy;
    method_results.l2norm = l2norm;
    method_results.relative_mass_error = relative_mass_error;
    method_results.absolute_momentum_error = absolute_momentum_error;
    method_results.relative_energy_error = relative_energy_error;
    method_results.relative_l2norm_error = relative_l2norm_error;
    method_results.initial_values = [M0, P0, E0, L2norm0];
    
    fprintf('    Initial values: M0=%.6e, P0=%.6e, E0=%.6e, L2norm0=%.6e\n', M0, P0, E0, L2norm0);
    fprintf('    Final errors: Mass=%.2e (rel), Momentum=%.2e (abs), Energy=%.2e (rel), L2norm=%.2e (rel)\n', ...
        relative_mass_error(end), absolute_momentum_error(end), relative_energy_error(end), relative_l2norm_error(end));
end
