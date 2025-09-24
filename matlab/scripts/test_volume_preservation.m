close all;
clear all;
clc

% TEST_VOLUME_PRESERVATION Simple volume preservation test
%
% Tests volume preservation for different resolutions using a single interpolation method.
% Maps phi(x,v) = (x + sin(pi*v), v) on domain (x,v) ∈ [0,2π] × [-6,6].

fprintf('=== Simple Volume Preservation Test ===\n\n');
% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS

%% Parameters
Lx = 2*pi;  % Spatial domain [0, 2π]
Lv = 12;    % Velocity domain [-6, 6] -> [0, 12] after shift
v_shift = 6; % Shift to make velocity domain [0, 12]

% Test parameters
resolutions = [1024, 32, 64, 128, 256, 512];
max_compositions = 20;
amplitude = 0.1; % Small amplitude for the map

% Colors for plotting
colors = {'b', 'r', 'g', 'm', 'c', 'k'};
markers = {'o', 's', '^', 'd', 'v', 'p'};

%% Test Map Definition
fprintf('Test map: phi(x,v) = (x + sin(pi*v), v)\n');
fprintf('Domain: (x,v) ∈ [0,2π] × [-6,6] (shifted to [0,2π] × [0,12])\n\n');

%% Test Volume Preservation
fprintf('=== Volume Preservation Test ===\n');

% Storage for results
results = struct();
results.resolutions = resolutions;
results.max_errors = zeros(length(resolutions), max_compositions);
results.mean_errors = zeros(length(resolutions), max_compositions);

% Test each resolution
for res_idx = 1:length(resolutions)
    resolution = resolutions(res_idx);
    Nx = resolution;
    Nv = resolution;
    
    fprintf('Testing resolution %dx%d...\n', Nx, Nv);
    
    % Create grid
    [grid] = make_periodic_grid(Lx,Lv,Nx,Nv);
    
    % Create params structure
    params = struct();
    params.Lx = Lx;
    params.Lv = Lv;
    V = grid.V;
    X = grid.X;
    % Create maps for composition
    V_original = V - v_shift;
    Maps = zeros(Nv, Nx, 2, max_compositions);
    for i = 1:max_compositions
        Maps(:,:,1,i) = X + amplitude * sin(pi * V_original);
        Maps(:,:,2,i) = V;
    end
    
        % Test composition with increasing number of maps
        for n_maps = 1:max_compositions
            % Compose n_maps maps
            Maps_subset = Maps(:,:,:,1:n_maps);
        
            Map_composed = compose_maps_numerical(Maps_subset, grid, params);
             
            if res_idx ==1 
                Map_composed_ref(:,:,:,n_maps) = Map_composed;
                X_ref = X;
                V_ref = V;
                grid_ref = grid;
            end
           [X_star, V_star] = evaluate_map(Map_composed, grid, params, X_ref, V_ref);
          
             results.l2_errors(res_idx, n_maps) = norm(X_star - Map_composed_ref(:,:,1, n_maps), 'fro') + norm(V_star - Map_composed_ref(:,:,2,n_maps), 'fro');
            
                
            % Compute volume preservation error
            [detJ, ~, ~, ~, ~] = jacobian_determinant(X_star, V_star, grid_ref, 'fourier');
            
            incomp_error = abs(log(detJ));
            max_error = max(incomp_error(:));
            mean_error = mean(incomp_error(:));
            
            results.max_errors(res_idx, n_maps) = max_error;
            results.mean_errors(res_idx, n_maps) = mean_error;
            
            if mod(n_maps, 5) == 0 || n_maps == 1 || n_maps == max_compositions
                fprintf('  %d maps: max error = %.2e, mean error = %.2e\n', ...
                    n_maps, max_error, mean_error);
            end
        end
        
  
    fprintf('\n');
end

%% Visualization
fprintf('=== Results Visualization ===\n');

% Plot: Error accumulation with map composition for different resolutions
fig3 = figure(56);
for res_idx = 1:length(resolutions)
    resolution = resolutions(res_idx);
   
    semilogy(1:max_compositions, results.max_errors(res_idx, :), ...
        [colors{res_idx} markers{res_idx} '-'], 'LineWidth', 2, ...
        'DisplayName', sprintf('%dx%d', resolution, resolution));
    hold on;
end
xlabel('Number of Composed Maps');
ylabel('Max Incompressibility Error');
title('Error Accumulation with Map Composition');
legend('Location', 'northwest');
grid on;
hold off
save_fig_tikz("../images/volumepreserving_upsample", fig3)


%% Plot L2 error

fig2 = figure(69);
nl=0;
for n_maps = [1, max_compositions]
    nl = nl+1;
    loglog(resolutions(2:end), results.l2_errors(2:end, n_maps),markers{nl})
    hold on;
    
    legend_entries{nl} = sprintf("%d maps", n_maps);
end

loglog(sort(resolutions),1e8*(1./sort(resolutions)).^(4),'k--')
legend_entries{end+1}="$\mathcal{O}(N^-4)$";
xlabel('1/N');
ylabel('l2 error');
legend(legend_entries);
xlim([resolutions(2)-2, resolutions(end)+30])
grid on;


save_fig_tikz("../images/volumepreserving_l2error", fig2)