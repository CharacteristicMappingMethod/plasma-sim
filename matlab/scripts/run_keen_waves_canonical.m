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
params.method = "CMM";
%%

start_from_existing = 1;
% Check if data already exists
params.data_dir = "../data/keen_waves_"+drive+"_method_" +params.method+"_dt0.05";
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

%% Get grid information
dom = params.grids(1).dom;
x = params.grids(1).X;
v = params.grids(1).v;
if params.method=="CMM_vargrid"
X = params.grids(1).Xsample_grid;
V = params.grids(1).Vsample_grid;
else
    X = params.grids(1).X;
    V = params.grids(1).V;
end


% Compute initial condition
f0 = params.f0(X, V);

% Compute f - f_0
f_deviation = fs_final- f0;

%% Create figure
fig_name = sprintf("../images/%s_Tend%d_dt%.3f_Nremap%d_%s", params.mycase, params.Tend, params.dt, params.N_remap,params.method);
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Initial distribution function
%subplot(2, 1, 1);
x = X(1,:);
v = V(:,1);

uimagesc(x, v, f_deviation);
set(gca, 'YDir','normal');
hold on
dN =2;
dNx = 2;
%mesh(X(1:dN:end,1:dNx:end), V(1:dN:end,1:dNx:end), 0*f_deviation(1:dN:end,1:dNx:end), 'EdgeColor', 'k', 'FaceAlpha', 0);
Xequi= params.grids(1).X;
Vequi = params.grids(1).V;

mesh(Xequi(1:dN:end,1:dNx:end), Vequi(1:dN:end,1:dNx:end), 0*Xequi(1:dN:end,1:dNx:end), 'EdgeColor', 'k', 'FaceAlpha', 0);
%shading flat;
colorbar;

colormap(jet)
title('$\delta f = f - f_0$', 'FontSize', 26);
xlabel('$x$', 'FontSize', 24); 
ylabel('$v$', 'FontSize', 24);
ylim([-5,5])

% Set axis font size
set(gca, 'FontSize', 32);

% Save as PNG with 600 DPI
print(fig_name+"_full_domain.png", '-dpng', '-r600');

% Save as TikZ (your original)
%save_fig_tikz(fig_name+"_full_domain")

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


figure(6)
species_name = params.species_name(1);  % Assuming single species
csv_file = fullfile(params.data_dir, species_name + ".csv");

fprintf('Loading conservation data from: %s\n', csv_file);
data_table = readtable(csv_file);

time = data_table.time;
rho_modes(:,1) = data_table.rho_1;
rho_modes(:,2) = data_table.rho_2;
rho_modes(:,3) = data_table.rho_3;
rho_modes(:,4) = data_table.rho_4;
rho_modes(:,5) = data_table.rho_5;
plot(time(1:1:end), rho_modes(1:1:end,:),"LineWidth",1);
title('Density Modes - Drive Frequency Component');
xlabel('Time $t$', 'Interpreter', 'latex');
ylabel('$|\hat{\rho}(k, t)|$', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
xlim([0,time(end)])
legend("$k=1$","$k=2$","$k=3$","$k=4$","$k=5$")
% Save the plot
save_fig_tikz(fig_name + "_density_modes");


%% Make a video of zooming into the fine structures with focus point (xstar, vstar)

xstar = params.Lx/2;
vstar = sum(abs(v_zoom_range))/2;

% Create output directory
output_dir = sprintf("../images/%s_zoom_sequence", params.mycase);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Define zoom sequence parameters
num_frames = 50;  % Number of frames in the sequence
zoom_start = 0.7;  % Starting zoom factor (smaller = more zoomed out)
zoom_end = 1.0;    % Ending zoom factor (1.0 = original zoom level)

% Create zoom factors (exponential progression for smooth zoom)
zoom_factors = linspace(zoom_start, zoom_end, num_frames);

% Define the focus point for zooming
xstar = params.Lx/2;
vstar = sum(abs(v_zoom_range))/2;

fprintf('Creating %d zoom sequence images...\n', num_frames);

for i = 1:num_frames
    % Calculate current zoom level
    current_zoom = zoom_factors(i);
    
    % Calculate zoomed domain
    zoom_width = (v_zoom_range(2) - v_zoom_range(1)) * current_zoom;
    zoom_height = (dom(3) - dom(1)) * current_zoom;
    
    % Center the zoom around the focus point
    v_zoom_min = vstar - zoom_width/2;
    v_zoom_max = vstar + zoom_width/2;
    x_zoom_min = xstar - zoom_height/2;
    x_zoom_max = xstar + zoom_height/2;
    
    % Create zoomed grid
    v_zoom_current = linspace(v_zoom_min, v_zoom_max, 256);
    x_zoom_current = linspace(x_zoom_min, x_zoom_max, 256);
    [X_zoom_current, V_zoom_current] = meshgrid(x_zoom_current, v_zoom_current);
    
    % Get zoomed data
    f_deviation_zoom_current = zoom(params, X_zoom_current, V_zoom_current) - params.f0(X_zoom_current, V_zoom_current);
    
    % Create figure
    fig = figure('Position', [100, 100, 1200, 800]);
    imagesc(x_zoom_current, v_zoom_current, f_deviation_zoom_current);
    set(gca, 'YDir', 'normal');
    colormap('jet');
    colorbar;
    clim(myclim);
    
    % Add title with zoom level
    title(sprintf('Zoom Level: %.2f', current_zoom), 'FontSize', 20);
    xlabel('$x$', 'FontSize', 18);
    ylabel('$v$', 'FontSize', 18);
    set(gca, 'FontSize', 16);
    
    % Save image with zero-padded numbering
    filename = sprintf('%s/zoom_%03d.png', output_dir, i);
    print(fig, filename, '-dpng', '-r300');
    
    % Close figure to save memory
    close(fig);
    
    % Progress indicator
    if mod(i, 10) == 0
        fprintf('Progress: %d/%d frames completed\n', i, num_frames);
    end
end

fprintf('Zoom sequence completed! Images saved to: %s\n', output_dir);