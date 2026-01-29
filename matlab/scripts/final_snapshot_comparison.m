% Final Snapshot Comparison Script
% Compares final distribution functions of different numerical methods
% for plasma physics simulations
%
% This script loads simulation data and creates a comparison plot
% showing the final distribution functions for NuFi, CMM, and PredCorr methods
% If data is not found, it will run the simulations automatically.
%
% Usage: Run this script to generate comparison plots for any test case
clear all; clc; close all;

% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS

PARAMS_two_stream;
params_base = params;
% Configuration
case_name = "two_stream";
Nsample = [1024];
Nmap = [64];

methods = ["NuFi", "CMM", "predcorr"];
Tend = 100;

fprintf('=== Final Snapshot Comparison ===\n');
fprintf('Test case: %s\n', case_name);
fprintf('Resolutions: %s\n', mat2str(Nsample));
fprintf('Methods: %s\n', strjoin(methods, ', '));
fprintf('End time: %d\n', Tend);


%% Load final snapshots for each method and resolution
fprintf('\nLoading final snapshots...\n');

% Storage for results
final_snapshots = cell(length(Nsample), length(methods));
grid_info = cell(length(Nsample), length(methods));

for i_res = 1:length(Nsample)
    for i_method = 1:length(methods)
        Nx = Nsample(i_res);
        method = methods(i_method);
        
        fprintf('Loading: Nx=%d, Method=%s\n', Nx, method);
        
        % Create data directory path for this simulation
        %data_dir = sprintf("../data/%s_zoom_study_N%d_N%d_%s", case_name,  Nsample(1), Nmap(1), method);
        data_dir = sprintf("../data/%s_Tend%d_%s", case_name, Tend, method);
        config_file = fullfile(data_dir, "config_data.mat");
        
        if exist(config_file, 'file')
            load(config_file, 'data', 'params');
            
            % Get the final snapshot (last time sample)
            final_snapshots{i_res, i_method} = data.fs(:,:,end,1); % Last time sample, first species
            grid_info{i_res, i_method} = params.grids(1);
            fprintf('  Successfully loaded data\n');
        else
            fprintf('  Config file not found - running simulation...\n');
            
            % Run simulation for this method and resolution
            run_simulation_for_comparison(case_name, Nsample, Nmap, method, Tend, params_base);
            
            % Try to load the data again
            if exist(config_file, 'file')
                load(config_file, 'data', 'params');
                final_snapshots{i_res, i_method} = data.fs(:,:,end,1);
                grid_info{i_res, i_method} = params.grids(1);
                fprintf('  Simulation completed and data loaded\n');
            else
                fprintf('  Error: Simulation failed to create data file\n');
                final_snapshots{i_res, i_method} = [];
                grid_info{i_res, i_method} = [];
            end
        end
    end
end

%% Create comparison plots for each resolution
for i_res = 1:length(Nsample)
    Nx = Nsample(i_res);
    fprintf('\nCreating comparison plot for Nx=%d...\n', Nx);
    
    % Create figure for this resolution
    fig = figure('Position', [200, 200, 1500, 500]);
    
    % Create subplots for each method
    for i_method = 1:length(methods)
        subplot(1, 3, i_method);
        
        if ~isempty(final_snapshots{i_res, i_method})
            % Get grid coordinates from the loaded parameters
            current_grid = grid_info{i_res, i_method};
            x_coords = current_grid.x;  % x coordinates
            v_coords = current_grid.v;  % v coordinates
            
            % Plot the distribution function with proper coordinates
            imagesc(x_coords, v_coords, final_snapshots{i_res, i_method});
            colorbar;
            %colormap('hot');
            
            % Set axis labels and title
            if i_method == 1
                ylabel('$v$', 'Interpreter', 'latex');
            else
                set(gca, 'YTick', []);
            end
            xlabel('$x$', 'Interpreter', 'latex');
            title(sprintf('\\textbf{%s}', methods{i_method}), 'FontSize', 20, 'Interpreter', 'latex');
            
            axis tight;
        else
            % If no data available, show a message
            text(0.5, 0.5, 'Data not available', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', 14);
            title(sprintf('\\textbf{%s}', methods{i_method}), 'FontSize', 20, 'Interpreter', 'latex');
        end
    end
    
    
    % Save the comparison plot
    fig_name = "../images/"+case_name+"_final_snapshots_N"+num2str(Nx)+"_Tend"+num2str(Tend);
    print(fig_name + ".png", '-dpng', '-r600');

    % Also save as TikZ if the function is available
    if exist('save_fig_tikz', 'file')
        save_fig_tikz(fig_name);
    end

    fprintf('Comparison plot saved as: %s\n', fig_name);
end

fprintf('\n=== Final Snapshot Comparison Complete ===\n');
fprintf('All comparison plots saved in ../images/\n');

%% Create zoom video focusing on [xf, vf] 
fprintf('\n=== Creating Zoom Video ===\n');

% Focus point
xf = 6;
vf = 2;

% Get the first resolution and check if we have NuFi and CMM data
i_res = 1;
Nx = Nsample(i_res);

% Check if we have data for the required methods
nufi_idx = find(methods == "NuFi");
cmm_idx = find(methods == "CMM");
predcorr_idx = find(methods == "predcorr");

% We need to reload the full parameter structures for zoom function
% Load NuFi parameters
data_dir_nufi = sprintf("../data/%s_Tend%d_%s", case_name, Tend, "NuFi");
config_file_nufi = fullfile(data_dir_nufi, "config_data.mat");
if exist(config_file_nufi, 'file')
    load(config_file_nufi, 'params');
    params_nufi = params;
    params_nufi.Efield_list = params.Efield_list;
    dom_nufi = params_nufi.grids(1).dom;
else
    fprintf('Error: NuFi data not found for zoom\n');
    return;
end

% Load CMM parameters
data_dir_cmm = sprintf("../data/%s_Tend%d_%s", case_name, Tend, "CMM");
config_file_cmm = fullfile(data_dir_cmm, "config_data.mat");
if exist(config_file_cmm, 'file')
    load(config_file_cmm, 'params');
    params_cmm = params;
    dom_cmm = params_cmm.grids(1).dom;
else
    fprintf('Error: CMM data not found for zoom\n');
    return;
end
%
% Define zoom sequence parameters
num_frames = 150;  % Number of frames in the zoom sequence
zoom_start = 1.0;  % Starting zoom factor (1.0 = original domain)
zoom_end = 0.005;   % Ending zoom factor (smaller = more zoomed in)

% Create zoom factors (exponential progression for smooth zoom)
zoom_factors = logspace(log10(zoom_start), log10(zoom_end), num_frames);
t
% Create output directory for zoom images
zoom_output_dir = "../images/zoom_" + case_name;
if ~exist(zoom_output_dir, 'dir')
    mkdir(zoom_output_dir);
    fprintf('Created zoom output directory: %s\n', zoom_output_dir);
end

fprintf('Creating %d zoom sequence images...\n', num_frames);

%% Static 4x4 zoom comparison figure
fprintf('\n=== Creating 4x4 zoom comparison figure ===\n');

num_zoom_levels = 4;
zoom_scales = [1.0, 0.25, 0.0625, 0.01];
zoom_titles = ["Full domain", "Zoom 1", "Zoom 2", "Zoom 3"];
N_zoom_plot = 512;

domain_x = [dom_nufi(1), dom_nufi(3)];
domain_v = [dom_nufi(2), dom_nufi(4)];

zoom_windows = cell(1, num_zoom_levels);
for idx = 1:num_zoom_levels
    if idx == 1
        zoom_windows{idx} = struct( ...
            'x_min', domain_x(1), ...
            'x_max', domain_x(2), ...
            'v_min', domain_v(1), ...
            'v_max', domain_v(2));
        continue;
    end
    width_x = diff(domain_x) * zoom_scales(idx);
    width_v = diff(domain_v) * zoom_scales(idx);

    x_min = max(domain_x(1), xf - width_x / 2);
    x_max = min(domain_x(2), xf + width_x / 2);
    v_min = max(domain_v(1), vf - width_v / 2);
    v_max = min(domain_v(2), vf + width_v / 2);

    zoom_windows{idx} = struct( ...
        'x_min', x_min, ...
        'x_max', x_max, ...
        'v_min', v_min, ...
        'v_max', v_max);
end

current_grid = grid_info{i_res, predcorr_idx};
x_coords = current_grid.x;
v_coords = current_grid.v;
f_predcorr = final_snapshots{i_res, predcorr_idx};

nufi_zoom_data = cell(1, num_zoom_levels);
cmm_zoom_data = cell(1, num_zoom_levels);
predcorr_zoom_data = cell(1, num_zoom_levels);

for idx = 1:num_zoom_levels
    x_lin = linspace(zoom_windows{idx}.x_min, zoom_windows{idx}.x_max, N_zoom_plot);
    v_lin = linspace(zoom_windows{idx}.v_min, zoom_windows{idx}.v_max, N_zoom_plot);
    [Xz, Vz] = meshgrid(x_lin, v_lin);

    try
        nufi_zoom_data{idx} = zoom(params_nufi, Xz, Vz);
    catch
        nufi_zoom_data{idx} = nan(size(Xz));
    end

    try
        cmm_zoom_data{idx} = zoom(params_cmm, Xz, Vz);
    catch
        cmm_zoom_data{idx} = nan(size(Xz));
    end

    predcorr_zoom_data{idx} = interp2(x_coords, v_coords, f_predcorr, Xz, Vz, 'linear');
end

fig_zoom_grid = figure('Position', [50, 50, 1900, 1200]);
tiled = tiledlayout(fig_zoom_grid, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
%colormap(tiled, turbo);

method_names = ["PredCorr", "NuFi", "CMM"];
method_data = {predcorr_zoom_data, nufi_zoom_data, cmm_zoom_data};

% Calculate global min/max for each method across all zoom levels
clim_per_method = cell(1, numel(method_names));
for method_idx = 1:numel(method_names)
    all_values = [];
    for col_idx = 1:num_zoom_levels
        data = method_data{method_idx}{col_idx};
        if ~isempty(data) && ~all(isnan(data(:)))
            all_values = [all_values; data(:)];
        end
    end
    if ~isempty(all_values)
        clim_per_method{method_idx} = [min(all_values), max(all_values)];
    else
        clim_per_method{method_idx} = [0, 1];
    end
end

for method_idx = 1:numel(method_names)
    for col_idx = 1:num_zoom_levels
        tile_position = num_zoom_levels * (method_idx - 1) + col_idx;
        ax = nexttile(tiled, tile_position);

        x_vals = linspace(zoom_windows{col_idx}.x_min, zoom_windows{col_idx}.x_max, N_zoom_plot);
        v_vals = linspace(zoom_windows{col_idx}.v_min, zoom_windows{col_idx}.v_max, N_zoom_plot);
        imagesc(x_vals, v_vals, method_data{method_idx}{col_idx});
        set(ax, 'YDir', 'normal');
        clim(ax, clim_per_method{method_idx});

        if method_idx == 1
            title(ax, sprintf('\\textbf{%s}', zoom_titles(col_idx)), 'Interpreter', 'latex', 'FontSize', 16);
        end

        if col_idx == 1
            ylabel(ax, sprintf('\\textbf{%s}', method_names(method_idx)), 'Interpreter', 'latex');
        end

        if method_idx < numel(method_names)
            set(ax, 'XTick', []);
        else
            xlabel(ax, '$x$', 'Interpreter', 'latex');
        end

        % Draw zoom indicator rectangles for all methods
        if col_idx < num_zoom_levels
            hold(ax, 'on');
            next_window = zoom_windows{col_idx + 1};
            rectangle(ax, 'Position', [next_window.x_min, next_window.v_min, ...
                next_window.x_max - next_window.x_min, next_window.v_max - next_window.v_min], ...
                'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
            hold(ax, 'off');
        end

        if col_idx == num_zoom_levels
            colorbar(ax);
        end
    end
end

%sgtitle(tiled, sprintf('\\textbf{%s final snapshot zooms}', case_name), 'Interpreter', 'latex', 'FontSize', 20);

zoom_grid_name = "../images/" + case_name + "_zoom_grid";
exportgraphics(fig_zoom_grid, zoom_grid_name + ".png", 'Resolution', 500);

%% ZOOM SEQUENCE
% Create figure once and reuse it
fig = figure('Position', [100, 100, 1600, 400]);

% Initialize subplots with initial data
% NuFi subplot
subplot(1, 3, 1);
% Get initial zoom data (full domain)
x_initial = linspace(dom_nufi(1), dom_nufi(3), 256);
v_initial = linspace(dom_nufi(2), dom_nufi(4), 256);
[X_initial, V_initial] = meshgrid(x_initial, v_initial);
f_nufi_initial = zoom(params_nufi, X_initial, V_initial);
imagesc(x_initial, v_initial, f_nufi_initial);
set(gca, 'YDir', 'normal');
colorbar;
title('\textbf{NuFi}', 'FontSize', 20,'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16);
ylabel('$v$', 'FontSize', 16);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% CMM subplot
subplot(1, 3, 2);
% Get initial zoom data (full domain)
x_initial = linspace(dom_cmm(1), dom_cmm(3), 256);
v_initial = linspace(dom_cmm(2), dom_cmm(4), 256);
[X_initial, V_initial] = meshgrid(x_initial, v_initial);
f_cmm_initial = zoom(params_cmm, X_initial, V_initial);
imagesc(x_initial, v_initial, f_cmm_initial);
set(gca, 'YDir', 'normal');
colorbar;
title('\textbf{CMM-NuFi}', 'FontSize', 20,'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16);
ylabel('$v$', 'FontSize', 16);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% PredCorr subplot
subplot(1, 3, 3);
% Get the full domain for predcorr
current_grid = grid_info{i_res, predcorr_idx};
x_coords = current_grid.x;
v_coords = current_grid.v;
f_predcorr = final_snapshots{i_res, predcorr_idx};

imagesc(x_coords, v_coords, f_predcorr);
set(gca, 'YDir', 'normal');
colorbar;
title('\textbf{PredCorr}', 'FontSize', 20,'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16);
ylabel('$v$', 'FontSize', 16);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% Get overall title handle
h_title = sgtitle('\textbf{Initial View}', 'FontSize', 20, 'Interpreter', 'latex');
%
for i = [1,1:num_frames]    
    % Calculate current zoom level
    current_zoom = zoom_factors(i);
    
    % Calculate zoomed domain size
    zoom_width_x = (dom_nufi(3) - dom_nufi(1)) * current_zoom;
    zoom_width_v = (dom_nufi(4) - dom_nufi(2)) * current_zoom;
    
    % Calculate center point for this frame (faster recentering to focus point)
    progress = (i - 1) / (num_frames - 1);  % 0 to 1
    
    % Use a power function to make recentering faster (reaches focus point earlier)
    % This makes the window move to focus point in first 1/3 of frames, then stays there
    recenter_progress = min(1, progress * 3);  % Reaches 1 at frame num_frames/3
    
    center_x = (1 - recenter_progress) * (dom_nufi(1) + dom_nufi(3)) / 2 + recenter_progress * xf;
    center_v = (1 - recenter_progress) * (dom_nufi(2) + dom_nufi(4)) / 2 + recenter_progress * vf;
    
    % Calculate zoomed domain boundaries
    x_min = center_x - zoom_width_x / 2;
    x_max = center_x + zoom_width_x / 2;
    v_min = center_v - zoom_width_v / 2;
    v_max = center_v + zoom_width_v / 2;
    
    % Create zoomed grid
    N_zoom = 1024;  % Resolution for zoomed view
    x_zoom = linspace(x_min, x_max, N_zoom);
    v_zoom = linspace(v_min, v_max, N_zoom);
    [X_zoom, V_zoom] = meshgrid(x_zoom, v_zoom);
    
    % Update NuFi plot
    subplot(1, 3, 1);
    %cla; % Clear current axes
    try
        f_nufi_zoom = zoom(params_nufi, X_zoom, V_zoom);
        imagesc(x_zoom, v_zoom, f_nufi_zoom);
        set(gca, 'YDir', 'normal');
        colorbar; % Recreate colorbar
        title('\textbf{NuFi}', 'FontSize', 20,'Interpreter', 'latex');
        xlabel('$x$', 'FontSize', 16);
        ylabel('$v$', 'FontSize', 16);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
    catch
        text(0.5, 0.5, 'NuFi zoom failed', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 16);
    end
    
    % Update CMM plot
    subplot(1, 3, 2);
    %cla; % Clear current axes
    try
        f_cmm_zoom = zoom(params_cmm, X_zoom, V_zoom);
        imagesc(x_zoom, v_zoom, f_cmm_zoom);
        set(gca, 'YDir', 'normal');
        colorbar; % Recreate colorbar
        title('\textbf{CMM-NuFi}', 'FontSize', 20,'Interpreter', 'latex');
        xlabel('$x$', 'FontSize', 16);
        ylabel('$v$', 'FontSize', 16);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
    catch
        text(0.5, 0.5, 'CMM zoom failed', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 16);
    end
    
    % Update PredCorr reference plot
    subplot(1, 3, 3);
    %cla; % Clear current axes
    try
        % Get the full domain for predcorr
        current_grid = grid_info{i_res, predcorr_idx};
        x_coords = current_grid.x;
        v_coords = current_grid.v;
        f_predcorr = final_snapshots{i_res, predcorr_idx};
        
        imagesc(x_coords, v_coords, f_predcorr);
        set(gca, 'YDir', 'normal');
        colorbar; % Recreate colorbar
        title('\textbf{PredCorr}', 'FontSize', 20,'Interpreter', 'latex');
        xlabel('$x$', 'FontSize', 16);
        ylabel('$v$', 'FontSize', 16);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        
        % Add a rectangle to show the zoom region
        hold on;
        rectangle('Position', [x_min, v_min, x_max-x_min, v_max-v_min], ...
                 'EdgeColor', 'white', 'LineWidth', 2, 'LineStyle', '--');
        hold off;
    catch
        text(0.5, 0.5, 'PredCorr plot failed', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 16);
    end
    
    % Update overall title
    set(h_title, 'String', sprintf('\\textbf{Zoom factor: %.1f}', 1/current_zoom));
    h_title.Interpreter="latex";
    % Save combined image
    filename = sprintf('%s/two_stream_fokus_%03d.png', zoom_output_dir, i);
    print(fig, filename, '-dpng', '-r300');
    
    
    % Progress indicator
    if mod(i, 5) == 0 || i == num_frames
        fprintf('Progress: %d/%d frames completed\n', i, num_frames);
    end
end

fprintf('Zoom sequence completed! Images saved to: %s\n', zoom_output_dir);
fprintf('Focus point: (%.1f, %.1f)\n', xf, vf);
fprintf('Zoom range: %.3f to %.3f\n', zoom_start, zoom_end);

%% Function to run simulation for a specific method and resolution
function run_simulation_for_comparison(case_name, Nsample, Nmap, method, Tend, base_params)
    % Create a copy of the current parameters
    params_copy = base_params;
    
    % Set method-specific parameters
    params_copy.method = method;
    params_copy.Nsample = Nsample;
    params_copy.Nmap = Nmap;
    params_copy.Tend = Tend;
    
    % Create data directory
    data_dir = sprintf("../data/%s_zoom_study_N%d_N%d_%s", case_name, Nsample(1), Nmap(1), method);
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end
    
    % Set data directory in parameters
    params_copy.data_dir = data_dir;
    
    fprintf('Running simulation: Nx=%d, Method=%s, Tend=%d\n', Nsample(1), method, Tend);
    
    % Run the simulation
    tic;
    [params, data] = Sim(params_copy);
    sim_time = toc;
    
    fprintf('Simulation completed in %.2f seconds\n', sim_time);
    
    % Save the simulation data
    config_file = fullfile(data_dir, "config_data.mat");
    save(config_file, 'params', 'data');
    fprintf('Data saved to %s\n', config_file);
end
