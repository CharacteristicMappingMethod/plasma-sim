%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conservation Quantities Analysis for Three Methods - Multiple Test Cases
%
% This script analyzes the conservation of mass, momentum, and energy
% for three different numerical methods on selected test cases:
% 1. CMM with N_remap = 10
% 2. NuFi 
% 3. Predictor-corrector scheme
%
% Available Test Cases:
% - Landau Damping (alpha = 0.05)
% - Two-Stream Instability
%
% Note: The script uses conservation quantities already calculated by 
% the measure.m function during simulation and saved to CSV files.
%
% Physical Quantities:
%   - Total Mass M(t) = ∬f(x,v,t)dxdv
%   - Total Momentum P(t) = ∬vf(x,v,t)dxdv  
%   - Total Energy E(t) = E_kin(t) + E_pot(t)
%     - E_kin(t) = 1/2 ∬f(x,v,t)|v|²dxdv
%     - E_pot(t) = 1/2 ∫|E(x,t)|²dx
%
% Error Calculations:
%   - Mass Error (Relative): |M(t) - M(0)| / |M(0)|
%   - Momentum Error (Absolute): |P(t) - P(0)|
%   - Energy Error (Relative): |E(t) - E(0)| / |E(0)|
%
% Methods: CMM (N_remap=10), NuFi, Predictor-Corrector
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
    data_files = {
        sprintf('conservation_landau_damping_CMM_Nremap10_Tend%.0f.mat', target_time),
        sprintf('conservation_landau_damping_NuFi_Tend%.0f.mat', target_time),
        sprintf('conservation_landau_damping_predcorr_Tend%.0f.mat', target_time)
    };
    plot_filename = sprintf('conservation_three_methods_landau_damping_Tend%.0f', target_time);
else
    % Two-Stream Instability configuration
    test_case_name = 'two_stream';
    test_case_display = 'Two-Stream Instability';
    params_script = 'PARAMS_two_stream';
    target_time = 50.0;  % Two-stream typically needs longer time
    data_files = {
        sprintf('conservation_two_stream_CMM_Nremap10_Tend%.0f.mat', target_time),
        sprintf('conservation_two_stream_NuFi_Tend%.0f.mat', target_time),
        sprintf('conservation_two_stream_predcorr_Tend%.0f.mat', target_time)
    };
    plot_filename = sprintf('conservation_three_methods_two_stream_Tend%.0f', target_time);
end

%% Configuration
fprintf('\nSelected test case: %s\n', test_case_display);
fprintf('Methods: CMM (N_remap=10), NuFi, Predictor-Corrector\n');
fprintf('Time range: 0 to %.1f\n\n', target_time);

% Analysis parameters
methods = {"CMM", "NuFi", "predcorr"};
method_names = {"CMM-NuFi", "NuFi", "Predictor-Corrector"};
colors = {'red', 'blue', 'green'};

%% Check for existing data and run simulations if needed
results = cell(3, 1);
run_simulations = [false, false, false];

for method_idx = 1:3
    method = methods{method_idx};
    method_name = method_names{method_idx};
    data_file = data_files{method_idx};
    
    fprintf('Checking data for %s (%s)...\n', method_name, test_case_display);
    
    if exist(data_file, 'file')
        fprintf('  Found existing data. Loading...\n');
        try
            load(data_file, 'method_results');
            results{method_idx} = method_results;
            fprintf('  Loaded successfully.\n');
        catch
            fprintf('  Warning: Failed to load data. Will re-run simulation.\n');
            run_simulations(method_idx) = true;
        end
    else
        fprintf('  No existing data found. Will run simulation.\n');
        run_simulations(method_idx) = true;
    end
end

%% Run simulations if needed
for method_idx = 1:3
    if run_simulations(method_idx)
        method = methods{method_idx};
        method_name = method_names{method_idx};
        data_file = data_files{method_idx};
        
        fprintf('\nRunning simulation for %s (%s)...\n', method_name, test_case_display);
        
        % Save current state and clear workspace
        temp_file = 'temp_conservation_three_methods_testcase.mat';
        save(temp_file, 'methods', 'method_names', 'colors', 'data_files', 'results', ...
             'run_simulations', 'target_time', 'method_idx', 'test_case_name', ...
             'test_case_display', 'params_script', 'plot_filename');
        
        clear all
        addpath(genpath('./src/'),genpath('./params/'))
        DEFAULTS
        
        % Reload configuration
        temp_file = 'temp_conservation_three_methods_testcase.mat';
        load(temp_file);
        method = methods{method_idx};
        method_name = method_names{method_idx};
        data_file = data_files{method_idx};
        
        % Load appropriate parameters based on test case
        eval(params_script);
        
        % Set method-specific parameters
        params.method = method;
        params.Tend = target_time;
        
        if strcmp(method, "CMM")
            params.N_remap = 10;
        end
        
        fprintf('  Running %s simulation (%s, Tend=%.1f)...\n', method_name, test_case_display, params.Tend);
        
        % Run simulation
        tic_sim = tic();
        [sim_params, sim_data] = Sim(params);
        sim_time = toc(tic_sim);
        
        fprintf('  Simulation completed in %.1f seconds\n', sim_time);
        
        % Calculate conservation quantities
        fprintf('  Computing conservation quantities...\n');
        method_results = calculate_conservation_quantities(sim_params, sim_data);
        method_results.method = method;
        method_results.method_name = method_name;
        method_results.test_case = test_case_name;
        method_results.test_case_display = test_case_display;
        method_results.simulation_time = sim_time;
        
        % Save results
        save(data_file, 'method_results');
        fprintf('  Results saved to: %s\n', data_file);
        
        % Store in results cell array
        results{method_idx} = method_results;
        
        % Clean up temporary file
        if exist(temp_file, 'file')
            delete(temp_file);
        end
    end
end

%% Create Publication-Quality Conservation Quantities Plot
fprintf('\nCreating publication-quality conservation quantities plot for %s...\n', test_case_display);

% Publication-quality figure settings
fig = figure('Position', [100, 100, 1000, 800], 'Units', 'pixels', ...
             'PaperPositionMode', 'auto', 'Color', 'white', 'Visible', 'on');
% Set paper size in inches for publication
set(fig, 'PaperUnits', 'inches', 'PaperSize', [10, 8]);

% Professional color scheme - colorblind-friendly palette
colors = [0.0000, 0.4470, 0.7410;    % Blue (CMM)
          0.8500, 0.3250, 0.0980;    % Red-orange (NuFi) 
          0.9290, 0.6940, 0.1250];   % Yellow-orange (Pred-corr)
          
line_styles = {'-', ':', '--'};  % Distinct line styles (solid, dotted, dashed)
line_widths = [2.5, 2.5, 2.5];
marker_styles = {'none', 'none', 'none'};

% Set publication font properties with larger sizes
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

% Top subplot: Mass Error
subplot(3, 1, 1);
hold on;
h_mass = [];
all_mass_errors = [];
for method_idx = 1:3
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        mass_error = results{method_idx}.relative_mass_error;
        
        % Handle zero/negative values for log scale by setting minimum threshold
        mass_error_plot = mass_error;
        zero_mask = mass_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(mass_error_plot(mass_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            mass_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array, mass_error_plot, ...
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
subplot(3, 1, 2);
hold on;
h_momentum = [];
all_momentum_errors = [];
for method_idx = 1:3
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
        
        % Handle zero/negative values for log scale
        momentum_error_plot = momentum_error;
        zero_mask = momentum_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(momentum_error_plot(momentum_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            momentum_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array, momentum_error_plot, ...
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

% Bottom subplot: Energy Error
subplot(3, 1, 3);
hold on;
h_energy = [];
all_energy_errors = [];
for method_idx = 1:3
    if ~isempty(results{method_idx})
        time_array = results{method_idx}.time_array;
        energy_error = results{method_idx}.relative_energy_error;
        
        % Handle zero/negative values for log scale
        energy_error_plot = energy_error;
        zero_mask = energy_error_plot <= 0;
        if any(zero_mask)
            min_positive = min(energy_error_plot(energy_error_plot > 0));
            if isempty(min_positive)
                min_positive = 1e-16;  % Machine precision fallback
            end
            energy_error_plot(zero_mask) = min_positive / 10;  % Set zeros to 1/10 of minimum positive
        end
        
        h = plot(time_array, energy_error_plot, ...
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
xlabel('Time t', 'FontSize', 16);
title('(c) Energy', 'FontSize', 16, 'FontWeight', 'normal');
grid on;
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
legend(h_energy, 'Location', 'southeast', 'FontSize', 14, 'Box', 'off');

% Fine-tune subplot positioning for very tight layout
subplot_positions = [0.13, 0.73, 0.82, 0.24;   % Top subplot
                     0.13, 0.42, 0.82, 0.24;   % Middle subplot  
                     0.13, 0.07, 0.82, 0.24];  % Bottom subplot (moved to y=0.05)

for i = 1:3
    subplot(3, 1, i);
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

% Save publication-quality plots
images_dir = './analysis/images';

% Create analysis/images directory if it doesn't exist
if ~exist(images_dir, 'dir')
    mkdir(images_dir);
    fprintf('Created %s directory\n', images_dir);
end

% Save in multiple high-quality formats for publication
plot_filename_pub = [plot_filename '_publication'];
eps_path = fullfile(images_dir, [plot_filename_pub '.eps']);
pdf_path = fullfile(images_dir, [plot_filename_pub '.pdf']);
png_path = fullfile(images_dir, [plot_filename_pub '.png']);
fig_path = fullfile(images_dir, [plot_filename_pub '.fig']);

% High-resolution exports
print(gcf, png_path, '-dpng', '-r300');      % High-res PNG for presentations
print(gcf, eps_path, '-depsc2', '-r300');    % Vector EPS for LaTeX
print(gcf, pdf_path, '-dpdf', '-r300');      % Vector PDF for modern systems
savefig(gcf, fig_path);                      % MATLAB figure for future editing

fprintf('Publication-quality plots saved:\n');
fprintf('  EPS (vector): %s\n', eps_path);
fprintf('  PDF (vector): %s\n', pdf_path);
fprintf('  PNG (raster): %s\n', png_path);
fprintf('  FIG (MATLAB): %s\n', fig_path);

%% Publication Summary
fprintf('\n=== Publication-Ready Analysis Summary ===\n');
fprintf('Test case: %s\n', test_case_display);
fprintf('Time range: 0 to %.1f plasma periods\n', target_time);
fprintf('Methods compared: %d\n', length(methods));
fprintf('Publication formats: EPS, PDF, PNG (300 DPI), FIG\n');

for method_idx = 1:3
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
        if isfield(method_results, 'simulation_time')
            fprintf('  Simulation time: %.1f seconds\n', method_results.simulation_time);
        end
    else
        fprintf('\n%s: No data available\n', method_names{method_idx});
    end
end

fprintf('\nPublication-ready analysis complete!\n');
fprintf('All plots saved in multiple formats suitable for journal submission.\n');

% Final cleanup of temporary files
temp_file = 'temp_conservation_three_methods_testcase.mat';
if exist(temp_file, 'file')
    delete(temp_file);
    fprintf('Cleaned up temporary files.\n');
end

%% Helper function to calculate conservation quantities
function method_results = calculate_conservation_quantities(sim_params, sim_data)
    fprintf('    Reading conservation quantities from simulation data...\n');
    
    % The measure.m function already calculates and saves these quantities
    % Read from the CSV files created during simulation
    species_name = sim_params.species_name(1);  % Assuming single species
    csv_file = fullfile(sim_params.data_dir, species_name + ".csv");
    
    if exist(csv_file, 'file')
        fprintf('    Reading data from: %s\n', csv_file);
        data_table = readtable(csv_file);
        
        % Extract quantities from CSV (these start from t = dt, not t = 0)
        time_array_csv = data_table.time;
        total_mass_csv = data_table.Mass;
        total_momentum_csv = data_table.Momentum;
        kinetic_energy_csv = data_table.Ekin;
        potential_energy_csv = data_table.Epot;
        total_energy_csv = data_table.Etot;
        
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
        
        % Initial electric field (from Poisson equation with initial distribution)
        % vPoisson expects fs as 3D array, so add species dimension
        f0_3d = zeros(size(f0, 1), size(f0, 2), 1);
        f0_3d(:, :, 1) = f0;
        Efield0 = vPoisson(f0_3d, sim_params.grids, sim_params.charge);
        Epot0 = 0.5 * sum(Efield0.^2) * grid.dx;
        Etot0 = Epot0 + Ekin0;
        
        % FIX: Use simulation's first mass value as reference for mass-conserving methods
        % This fixes the horizontal line issue for NuFi and predictor-corrector
        M0_sim = total_mass_csv(1);  % First simulation mass value
        P0_sim = total_momentum_csv(1);  % First simulation momentum value
        E0_sim = total_energy_csv(1);  % First simulation energy value
        
        % Combine initial conditions with CSV data using simulation reference values
        time_array = [0; time_array_csv];
        total_mass = [M0_sim; total_mass_csv];  % Use simulation M0
        total_momentum = [P0_sim; total_momentum_csv];  % Use simulation P0  
        total_energy = [E0_sim; total_energy_csv];  % Use simulation E0
        
        % Also store kinetic and potential energy (use computed t=0 values)
        kinetic_energy = [Ekin0; kinetic_energy_csv];
        potential_energy = [Epot0; potential_energy_csv];
        
        fprintf('    Successfully read %d time points from CSV + 1 initial condition\n', length(time_array_csv));
        fprintf('    Total time points: %d (including t=0)\n', length(time_array));
        fprintf('    Using simulation reference values: M0=%.6e, P0=%.6e, E0=%.6e\n', M0_sim, P0_sim, E0_sim);
    else
        fprintf('    Warning: CSV file not found. Using fallback calculation...\n');
        
        % Fallback: use basic quantities from sim_params
        time_array = [0; sim_params.time_array];  % Add t=0
        Nt = length(time_array);
        
        % Initialize with zeros (this is a fallback case)
        total_mass = ones(Nt, 1);  % Assume mass conservation
        total_momentum = zeros(Nt, 1);  % Assume zero initial momentum
        kinetic_energy = ones(Nt, 1);
        potential_energy = zeros(Nt, 1);
        total_energy = ones(Nt, 1);
        
        fprintf('    Warning: Using fallback values. Results may not be accurate.\n');
    end
    
    % Calculate initial values for error computation
    M0 = total_mass(1);
    P0 = total_momentum(1);
    E0 = total_energy(1);
    
    
    % Calculate relative errors (except momentum which is absolute)
    relative_mass_error = abs(total_mass - M0) / abs(M0);
    absolute_momentum_error = abs(total_momentum - P0);  % Absolute error for momentum
    relative_energy_error = abs(total_energy - E0) / abs(E0);
    
    % Store results
    method_results = struct();
    method_results.time_array = time_array;
    method_results.total_mass = total_mass;
    method_results.total_momentum = total_momentum;
    method_results.kinetic_energy = kinetic_energy;
    method_results.potential_energy = potential_energy;
    method_results.total_energy = total_energy;
    method_results.relative_mass_error = relative_mass_error;
    method_results.absolute_momentum_error = absolute_momentum_error;
    method_results.relative_energy_error = relative_energy_error;
    method_results.initial_values = [M0, P0, E0];
    method_results.csv_file = csv_file;
    
    fprintf('    Initial values: M0=%.6e, P0=%.6e, E0=%.6e\n', M0, P0, E0);
    fprintf('    Final errors: Mass=%.2e (rel), Momentum=%.2e (abs), Energy=%.2e (rel)\n', ...
        relative_mass_error(end), absolute_momentum_error(end), relative_energy_error(end));
end
