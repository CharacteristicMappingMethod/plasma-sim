% Run canonical drive Keen waves simulation and analyze results
clear all; clc; close all;


% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS
% Load canonical drive parameters
drive = "weak";
if drive == "canonical"
    PARAMS_keen_waves_canonical;
else
    PARAMS_keen_waves_weak;
end
params.method = "CMM_vargrid";
%%

start_from_existing = 1;
% Check if data already exists
params.data_dir = "../data/keen_waves_"+drive+"_method_" +params.method+"3";
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
[f_deviation_shift] = shift_field_periodic(f_deviation, 280);
uimagesc(x, v, f_deviation_shift);
set(gca, 'YDir','normal');
hold on
dN =2;
dNx = 2;
%mesh(X(1:dN:end,1:dNx:end), V(1:dN:end,1:dNx:end), 0*f_deviation(1:dN:end,1:dNx:end), 'EdgeColor', 'k', 'FaceAlpha', 0);
%Xequi= params.grids(1).X;
%Vequi = params.grids(1).V;

%mesh(Xequi(1:dN:end,1:dNx:end), Vequi(1:dN:end,1:dNx:end), 0*Xequi(1:dN:end,1:dNx:end), 'EdgeColor', 'k', 'FaceAlpha', 0);
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
save_fig_tikz(fig_name+"_full_domain")

%% Create zoomed view around v in [1.2, 1.6]

if drive == "weak"
    v_zoom_range = [1.2,1.6];
    myclim = [-0.024,0.024];
else
    v_zoom_range = [0.375, 2.25];
    myclim = [-0.215,0.20];
end

if  ~ isfield(params.grids(1), 'map') && params.method=="CMM_vargrid"
    params.grids(1).map = params.grids(1);
end

dom = params.grids(1).dom;
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
Nskip = 10;

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

% Load reference data
reference_dir = fullfile('KEEN_reference_data/',drive);
reference_modes = cell(5, 1);
reference_time = cell(5, 1);

for k = 1:5
    ref_file = fullfile(reference_dir, sprintf('mode%d.csv', k));
    if exist(ref_file, 'file')
        ref_data = readtable(ref_file);
        reference_time{k} = ref_data.x;
        reference_modes{k} = ref_data.y;
        fprintf('Loaded reference data for mode %d\n', k);
    else
        fprintf('Warning: Reference file %s not found\n', ref_file);
        reference_time{k} = [];
        reference_modes{k} = [];
    end
end

% Calculate rescaling factor based on mode 1 maximum values
if ~isempty(reference_modes{1}) && ~isempty(rho_modes(:,1))
    sim_max_mode1 = max(abs(rho_modes(:,1)));
    ref_max_mode1 = max(abs(reference_modes{1}));
    rescaling_factor = sim_max_mode1 / ref_max_mode1;
    fprintf('Rescaling factor (sim_max/ref_max): %.6f\n', rescaling_factor);

    % Apply rescaling to all reference modes
    rho_modes = rho_modes/rescaling_factor;

else
    fprintf('Warning: Cannot calculate rescaling factor - missing data\n');
    rescaling_factor = 1;
    end

    % Plot simulation and reference data in the same loop
    markers = {'o', 's', '^', 'd', 'v'};
    colors = {'red', 'green', 'blue', 'magenta', 'black'};
    hold on;

    for k = 1:5
        % Plot simulation data
        plot(time(1:Nskip:end), rho_modes(1:Nskip:end,k), "LineWidth", 1, ...
            'DisplayName', sprintf('$k=%d$ (CMM-NuFi)', k), 'Color', colors{k});

        % Plot reference data with same color
        if ~isempty(reference_modes{k})
            plot(reference_time{k}, reference_modes{k}, ':', 'Marker', markers{k}, ...
                'MarkerSize', 1, 'LineWidth', 1, 'Color', colors{k}, ...
                'DisplayName', sprintf('$k=%d$ (Afeyan)', k));
        end
    end

    %title('Density Modes - Drive Frequency Component (Simulation vs Reference)');
    xlabel('Time $t$', 'Interpreter', 'latex');
    ylabel('Fourier mode $\hat{\rho}(k, t)$', 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', 14);
    xlim([0,time(end)])
    legend('Location', 'best');
    hold off;

    % Save the plot
    save_fig_tikz(fig_name + "_density_modes_comparison");


    %% Make a video of zooming into the fine structures with focus point (xstar, vstar)

    xstar = params.Lx/2-0.55;
    vstar = sum(abs(v_zoom_range))/2+0.35;
    Nsample = 1024;

    % Create output directory
    output_dir = sprintf("../images/%s_zoom_sequence", params.mycase);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        fprintf('Created output directory: %s\n', output_dir);
    end

    % Define zoom sequence parameters
    num_frames = 4;  % Number of frames in the sequence

    % Create zoom factors separately for x and v directions
    zoom_factors_x = [1, 0.3, 0.1, 0.01, 0.005];
    zoom_factors_v = [1, 0.2, 0.03, 0.006, 0.0006];

    fprintf('Creating %d zoom sequence images...\n', num_frames);

    for i = 1:num_frames
        % Calculate current zoom levels for x and v separately
        current_zoom_x = zoom_factors_x(i);
        current_zoom_v = zoom_factors_v(i);

        % Calculate zoomed domain using separate zoom factors
        zoom_width = (dom(4) - dom(2)) * current_zoom_v;
        zoom_height = (dom(3) - dom(1)) * current_zoom_x;

        % Center the zoom around the focus point
        v_zoom_min = vstar - zoom_width/2;
        v_zoom_max = vstar + zoom_width/2;
        x_zoom_min = xstar - zoom_height/2-params.Lx/2;
        x_zoom_max = xstar + zoom_height/2-params.Lx/2;

        % Create zoomed grid
        v_zoom_current = linspace(v_zoom_min, v_zoom_max, Nsample);
        x_zoom_current = linspace(x_zoom_min, x_zoom_max, Nsample);
        [X_zoom_current, V_zoom_current] = meshgrid(x_zoom_current, v_zoom_current);

        % Get zoomed data
        f_deviation_zoom_current = zoom(params, X_zoom_current, V_zoom_current) - params.f0(X_zoom_current, V_zoom_current);

        % Create figure
        fig = figure('Position', [100, 100, 1200, 800]);
        imagesc(x_zoom_current+params.Lx/2, v_zoom_current, f_deviation_zoom_current);
        set(gca, 'YDir', 'normal');
        colormap('jet');
        colorbar;
        clim(myclim);
        hold on;

        % Draw dashed frame showing next zoom window (if not the last frame)
        if i < num_frames
            % Calculate next zoom window using separate zoom factors
            next_zoom_x = zoom_factors_x(i+1);
            next_zoom_v = zoom_factors_v(i+1);
            next_zoom_width = (dom(4) - dom(2)) * next_zoom_v;
            next_zoom_height = (dom(3) - dom(1)) * next_zoom_x;
            next_v_zoom_min = vstar - next_zoom_width/2;
            next_v_zoom_max = vstar + next_zoom_width/2;
            next_x_zoom_min = xstar - next_zoom_height/2;% - params.Lx/2;
            next_x_zoom_max = xstar + next_zoom_height/2;% - params.Lx/2;

            % Draw dashed rectangle for next zoom window
            % Only draw if the next window is within the current view
            if next_x_zoom_min-params.Lx/2 >= x_zoom_min && next_x_zoom_max-params.Lx/2 <= x_zoom_max && ...
                    next_v_zoom_min >= v_zoom_min && next_v_zoom_max <= v_zoom_max
                rectangle('Position', [next_x_zoom_min, next_v_zoom_min, ...
                    next_x_zoom_max - next_x_zoom_min, ...
                    next_v_zoom_max - next_v_zoom_min], ...
                    'EdgeColor', 'black', 'LineWidth', 4, 'LineStyle', '--');
            end
        end
        hold off;

        % Add title with zoom levels for x and v
        %title(sprintf('Zoom Level: x=%.4f, v=%.4f', current_zoom_x, current_zoom_v), 'FontSize', 20);
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