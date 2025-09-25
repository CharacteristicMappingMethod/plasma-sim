% Run canonical drive Keen waves simulation and analyze results
clear all; clc; close all;


% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS
% Load canonical drive parameters
drive = "canonical";
if drive == "canonical"
    PARAMS_keen_waves_canonical;
else
    PARAMS_keen_waves_weak;
end

start_from_existing = 1;
% Check if data already exists
params.data_dir = "../data/keen_waves_"+drive+"_method_" +params.method + "_dt0.25/";
data_file = params.data_dir+"/data_Tend"+num2str(params.Tend)+".mat";

if exist(data_file, 'file') && start_from_existing
    fprintf('Loading existing simulation data...\n');
    load(data_file, 'params', 'data');
    fs_final = data.fs(:,:,end);
    time_final = data.time(end);
    fprintf('Data loaded successfully!\n');
    fprintf("data: "+ data_file)
else
    fprintf("Starting Keen waves "+ drive+ " drive simulation...\n");
    
    % Run simulation
    tic;
    [params, data] = Sim(params);
    sim_time = toc;
    
    fprintf('Simulation completed in %.2f seconds\n', sim_time);
    
    % Save data
    save(data_file, 'params', 'data');
    fprintf('Data saved to %s\n', data_file);
    
    % Load the final data
    fs_final = data.fs(:,:,end);
    time_final = data.time(end);
end

% Get grid information
dom = params.grids(1).dom;
x = params.grids(1).x;
v = params.grids(1).v;
X = params.grids(1).X;
V = params.grids(1).V;


% Compute initial condition
f0 = params.f0(X, V);

% Compute f - f_0
f_deviation = fs_final - f0;

%% Create figure
fig_name = sprintf("../images/%s_Tend%d_dt%.1f_Nremap%d", params.mycase, params.Tend, params.dt, params.N_remap);
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Initial distribution function
%subplot(2, 1, 1);
imagesc(x, v, f_deviation);
set(gca, 'YDir','normal');
shading flat;
colorbar;
colormap(jet)
title('$\delta f = f - f_0$');
xlabel('$x$'); ylabel('$v$');
save_fig_tikz(fig_name+"_full_domain")
% Compute charge density
dx = x(2) - x(1);
rho = sum(fs_final, 1) * dx;
rho0 = sum(f0, 1) * dx;

%% Create zoomed view around v in [1.2, 1.6]

if drive == "weak"
v_zoom_range = [1.2,1.6];
myclim = [-0.024,0.024];
else
v_zoom_range = [0.375, 2.25];
myclim = [-0.215,0.20];
end
v_zoom = linspace(v_zoom_range(1), v_zoom_range(2), 1048);
x_zoom = linspace(dom(1), dom(3), 1048);
[X_zoom, V_zoom] = meshgrid(x_zoom, v_zoom);

% Use zoom function to get detailed view
f_deviation_zoom = zoom(params, X_zoom, V_zoom);

% Plot 2: Zoomed view of f_deviation
fig = figure(5)
fig.Position = [100, 100, 1200, 800];
imagesc(x_zoom, v_zoom, f_deviation_zoom-params.f0(X_zoom,V_zoom));
set(gca, 'YDir','normal');
colormap("jet")

colorbar;
clim(myclim)
title('Zoom');
xlabel('$x$'); ylabel('$v$');


save_fig_tikz(fig_name+"_zoom")


%% compute density

