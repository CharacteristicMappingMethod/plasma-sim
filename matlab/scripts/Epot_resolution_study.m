% Run or analyze Landau damping simulations with different space
% resolutions Nx=Nv = 32, 128, 256 for both CMM and NuFi methods If CSV
% files exist, reads data from them; otherwise runs simulations Plot
% potential energy over time in log scale
clear all; clc; close all;

% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS

case_name = "two_stream";
%case_name = "landau_damping";

if case_name == "landau_damping"
    fprintf('Selected: Two Stream Instability\n');
    % Load Landau damp         ing parameters
    PARAMS_landau_damping;
    % Define resolutions and methods to test
    resolutions = [32, 128, 256];
    resolution_map = [32];
    methods = ["NuFi", "CMM", "predcorr", "CMM_vargrid"];
        methods = ["NuFi", "CMM_vargrid"];
    Tend = 100;
elseif case_name == "two_stream"
    fprintf('Selected: Two Stream Instability\n');
    PARAMS_two_stream;
    resolutions = [ 128*2];
    resolution_map = [64];
    methods = ["NuFi", "CMM", "predcorr", "predcorr_3rd","CMM_vargrid"];
    Tend = 100;
else
    error("wrong case_name")
end


use_existing = 1;

% Storage for results
results = struct();
results.resolutions = resolutions;
results.methods = methods;
results.Tend = Tend;
results.epot_data = cell(length(resolutions), length(methods));
results.time_data = cell(length(resolutions), length(methods));

fprintf('Starting study...\n');
fprintf('Resolutions: %s\n', mat2str(resolutions));
fprintf('Methods: %s\n', strjoin(methods, ', '));
fprintf('End time: %d\n', Tend);
%% Run simulations or read existing data for each resolution and method combination
for i_res = 1:length(resolutions)
    for i_method = 1:length(methods)
        Nx = resolutions(i_res);
        Nv = resolutions(i_res);
        method = methods(i_method);

        fprintf('\n=== Processing: Nx=Nv=%d, Method=%s ===\n', Nx, method);

        % Create data directory path for this simulation
        data_dir = sprintf("../data/%s_resolution_study_N%d_%s", case_name,Nx, method);
        if method == "predcorr" 

        end
        csv_file = fullfile(data_dir, "electrons.csv");

        % Check if CSV file already exists
        if exist(csv_file, 'file') && use_existing
            fprintf('Loading existing potential energy data from: %s\n', csv_file);
            data_table = readtable(csv_file);
            time = data_table.time;
            Epot = data_table.Epot;

            fprintf('Successfully loaded data: %d time points\n', length(time));
        else
            fprintf('CSV file not found. Running simulation...\n');

            % Set up parameters for this simulation
            params_current = params;
            params_current.Nx = Nx;
            params_current.Nv = Nv;
            params_current.method = method;
            params_current.Tend = Tend;
            params_current.data_dir = data_dir;

            if method == "CMM_vargrid"
                params_current.Nsample = [Nx, Nx];
                params_current.Nmap = resolution_map;
            end
            % Create data directory if it doesn't exist
            if ~exist(data_dir, 'dir')
                mkdir(data_dir);
            end

            % Run simulation
            tic;
            [params_sim, data_sim] = Sim(params_current);
            sim_time = toc;

            fprintf('Simulation completed in %.2f seconds\n', sim_time);

            % If Epot is not directly available, try to read from CSV that
            % should have been created
            if exist(csv_file, 'file')
                fprintf('Reading potential energy from newly created CSV file\n');
                data_table = readtable(csv_file);
                time = data_table.time;
                Epot = data_table.Epot;
            else
                fprintf('Warning: Could not find potential energy data after simulation\n');
                time = [];
                Epot = [];
            end

        end

        % Store results
        results.epot_data{i_res, i_method} = Epot;
        results.time_data{i_res, i_method} = time;

        fprintf('Completed: Nx=Nv=%d, Method=%s\n', Nx, method);
    end
end

%% Create log plot of potential energy over time
fprintf('\nCreating potential energy plots...\n');

% Create figure for potential energy
fig1 = figure('Position', [100, 100, 1200, 800]);

% Define colors and line styles
colors = {'b', 'r', 'g', 'k', 'y'};
line_styles = {'-', 'none', ':','--','none'};
markers = {'none', 's', 'o','none','x'};
% Darker colors for markers (second and third methods)
marker_colors = {[0, 0, 0.6], [0.6, 0, 0], [0, 0.6, 0],  [0.6, 0.6, 0.6]}; % Darker blue, red, green

% Theoretical damping rate
gamma_E = 0.153359*2;

% Plot all methods and resolutions in the same plot
hold on;

for i_method = 1:length(methods)
    method = methods(i_method);

    for i_res = 1:length(resolutions)
        Nx = resolutions(i_res);
        Epot = results.epot_data{i_res, i_method};
        time = results.time_data{i_res, i_method};

        if ~isempty(Epot) && ~isempty(time)
            % Plot in log scale with different line styles for methods
            if method == "CMM"  % Second method (CMM) - use darker marker colors
                semilogy(time, abs(Epot), 'Color', colors{i_method}, ...
                    'LineStyle', line_styles{i_method}, 'LineWidth',1,...
                    'Marker', markers{i_method}, 'MarkerSize',1,...
                    'MarkerFaceColor', marker_colors{i_method}, 'MarkerEdgeColor', marker_colors{i_method}, ...
                    'DisplayName', sprintf('$\\Nsample=%d$ (CMM-NuFi)', Nx));
            elseif method=="predcorr"  % Third method (predcorr) - use dotted lines with markers
                semilogy(time-params.dt/2, abs(Epot), 'Color', colors{i_method}, ...
                    'LineStyle', line_styles{i_method}, 'LineWidth',1,...
                    'Marker', markers{i_method}, 'MarkerSize',1,...
                    'MarkerFaceColor', marker_colors{i_method}, 'MarkerEdgeColor', marker_colors{i_method}, ...
                    'DisplayName', sprintf('$\\Nsample=%d$ (%s)', Nx, method));
            elseif method=="NuFi"  % First method (NuFi) - use regular colors
                semilogy(time, abs(Epot), 'Color', colors{i_method}, ...
                    'LineStyle', line_styles{i_method},'LineWidth',1, ...
                    'Marker', markers{i_method}, ...
                    'DisplayName', sprintf('$\Nsample=%d$ (%s)', Nx, method));
            elseif method=="CMM_vargrid"  % First method (NuFi) - use regular colors
                semilogy(time(1:4:end), abs(Epot(1:4:end)), 'Color', colors{i_method}, ...
                    'LineStyle', line_styles{i_method},'LineWidth',1, ...
                    'Marker', markers{i_method}, ...
                    'DisplayName', sprintf('$\\Nsample=%d$ $\\Nmapgrid=%d$ (CMM-NuFi)', Nx, resolution_map));
            end
        else
            fprintf('Warning: No data available for Nx=%d, Method=%s\n', Nx, method);
        end
    end
end

% Add theoretical decay line only for landau_damping case
if case_name == "landau_damping"
    time_theory = linspace(0, Tend, 1000);
    initial_epot = 0.000298599;
    epot_theory = initial_epot * exp(-gamma_E * (time_theory-2.7));
    semilogy(time_theory, epot_theory, 'k--', ...
        'DisplayName', sprintf('Theory: $\\exp(-\\gamma_E t)$, $\\gamma_E = %.6f$', gamma_E));
end

xlabel('Time $t$', 'Interpreter', 'latex');
ylabel('$E_\text{pot}(t)$', 'Interpreter', 'latex');
legend('Location', 'southeastoutside');
grid on;
xlim([0, Tend]);
set(gca, 'YScale', 'log');

% Save the potential energy plot
fig_name = "../images/"+case_name+"_resolution_study_Tend"+num2str(Tend);
print(fig_name + ".png", '-dpng', '-r600');

% Also save as TikZ if the function is available
if exist('save_fig_tikz', 'file')
    save_fig_tikz(fig_name);
end

fprintf('Potential energy plot saved as: %s\n', fig_name);

