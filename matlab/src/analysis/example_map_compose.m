%%% MATLAB Script to Define and Compose 2D Maps - Periodic Displacement Fields

clear; clc; close all;

%% 1. Configuration and Grid Initialization
% ========== USER PARAMETERS ==========
N_MAPS = 4;  % Change this to any number of maps you want
VISUALIZE_INDIVIDUAL_MAPS = true;  % Set to false to skip individual map plots

% Set up grid using CMM conventions: x in [0, Lx), v in [-Lv, Lv)
Nx = 2^7;
Nv = 2^7;
Lx = 2*pi;  % Spatial domain length (typical for plasma simulations)
Lv = 6;     % Velocity domain half-width

% Create grid using CMM convention
x = (0:Nx-1) * Lx / Nx;  % x in [0, Lx)
v = (0:Nv-1) * 2*Lv / Nv - Lv;  % v in [-Lv, Lv)
[X, V] = meshgrid(x, v);

% Grid spacing
dx = x(2) - x(1);
dv = v(2) - v(1);

%% 2. Create N Maps Dynamically with PERIODIC DISPLACEMENT FIELDS
fprintf('Creating %d maps with periodic displacement fields...\n', N_MAPS);

% Initialize map storage
Maps_cell = cell(N_MAPS, 1);
Maps = zeros(Nv, Nx, 2, N_MAPS);

% Define map parameters ensuring PERIODIC DISPLACEMENTS
% For X domain [0, 2π]: frequencies must be 2π-periodic -> freq_x = n (integer)
% For V domain [-Lv, Lv]: frequencies must be 2Lv-periodic -> freq_v = n*π/Lv
map_params = [
    % [amp1, freq1_x, freq1_v, amp2, freq2_x, func1, func2]
    % All frequencies chosen to ensure periodic displacement fields
    
    % === HIGH AMPLITUDE MAPS (periodic) ===
    0.8, 1, pi/Lv, 1.5, 1, 1, 2;          % Map 1: Large amp, fundamental modes
    0.6, 2, 2*pi/Lv, 0.4, 3, 2, 1;        % Map 2: Higher harmonics
    
    % === ASYMMETRIC MAPS (periodic) ===
    0.1, 1, pi/Lv, 2.0, 1, 1, 1;          % Map 3: FIXED - now uses π/Lv (fundamental V mode)
    1.2, 1, pi/Lv, 0.05, 2, 2, 2;         % Map 4: Huge X, tiny V displacement
    
    % === MULTI-HARMONIC MAPS (periodic) ===
    0.3, 4, 3*pi/Lv, 0.3, 5, 1, 2;        % Map 5: High harmonics in both directions
    0.4, 1, pi/Lv, 0.6, 1, 2, 1;          % Map 6: FIXED - now uses π/Lv instead of π/(4*Lv)
    
    % === MIXED HARMONIC PATTERNS (periodic) ===
    0.5, 3, 2*pi/Lv, 0.5, 2, 1, 1;        % Map 7: 3rd and 2nd harmonics
    0.5, 2, 3*pi/Lv, 0.5, 4, 2, 2;        % Map 8: Even harmonics
    
    % === VELOCITY-DEPENDENT X MAPS (periodic) ===
    0.7, 1, 4*pi/Lv, 0.2, 2, 1, 2;        % Map 9: Strong V-dependence in X
    0.2, 3, 2*pi/Lv, 0.8, 1, 2, 1;        % Map 10: 3rd harmonic X, fundamental V
    
    % === FINE STRUCTURE MAPS (periodic) ===
    0.3, 6, 5*pi/Lv, 0.3, 7, 1, 1;        % Map 11: High mode numbers
    0.4, 8, pi/Lv, 0.2, 6, 2, 2;          % Map 12: Very high X harmonic
    
    % === FUNDAMENTAL + HARMONIC COMBINATIONS (periodic) ===
    0.25, 1, pi/Lv, 0.35, 1, 1, 2;        % Map 13: Pure fundamentals
    0.3, 2, 2*pi/Lv, 0.4, 3, 2, 1;        % Map 14: 2nd and 3rd harmonics
    
    % === MODERATE HARMONICS (periodic) ===
    0.2, 3, 3*pi/Lv, 0.3, 4, 1, 1;        % Map 15: 3rd and 4th harmonics
    0.35, 5, pi/Lv, 0.25, 2, 2, 2;        % Map 16: 5th and 2nd harmonics
    
    % === LOW AMPLITUDE, HIGH HARMONICS (periodic) ===
    0.02, 10, 6*pi/Lv, 0.03, 12, 1, 2;    % Map 17: Very high harmonics, tiny amplitude
    0.01, 15, 4*pi/Lv, 0.015, 8, 2, 1;    % Map 18: Ultra high harmonics
    
    % === MIXED SCALE COMBINATIONS (periodic) ===
    0.05, 12, pi/Lv, 0.8, 1, 1, 1;        % Map 19: High X harmonic, fundamental V
    0.9, 1, 5*pi/Lv, 0.04, 10, 2, 2;      % Map 20: Fundamental X, high V harmonic
];

% Create each map with GUARANTEED PERIODIC DISPLACEMENTS
for i = 1:N_MAPS
    % Get parameters for this map
    amp1 = map_params(i,1);
    freq1_x = map_params(i,2);         % Integer for 2π-periodicity in X
    freq1_v = map_params(i,3);         % Multiple of π/Lv for 2Lv-periodicity in V
    amp2 = map_params(i,4);
    freq2_x = map_params(i,5);         % Integer for 2π-periodicity in X
    func1 = map_params(i,6);
    func2 = map_params(i,7);
    
    % Apply function choices with PERIODIC DISPLACEMENT FIELDS
    if func1 == 1  % sin-sin
        % ΔX = amp1 * sin(freq1_x * X) * sin(freq1_v * V) is periodic in both X and V
        Map_X = X + amp1 * sin(freq1_x * X) .* sin(freq1_v * V);
    else  % cos-cos
        % ΔX = amp1 * cos(freq1_x * X) * cos(freq1_v * V) is periodic in both X and V
        Map_X = X + amp1 * cos(freq1_x * X) .* cos(freq1_v * V);
    end
    
    if func2 == 1  % sin
        % ΔV = amp2 * sin(freq2_x * X) is periodic in X
        Map_V = V + amp2 * sin(freq2_x * X);
    else  % cos
        % ΔV = amp2 * cos(freq2_x * X) is periodic in X
        Map_V = V + amp2 * cos(freq2_x * X);
    end
    
    % Store the map
    Maps_cell{i} = zeros(Nv, Nx, 2);
    Maps_cell{i}(:,:,1) = Map_X;
    Maps_cell{i}(:,:,2) = Map_V;
    Maps(:,:,:,i) = Maps_cell{i};
    
    % Visualize DISPLACEMENT FIELDS if requested
    if VISUALIZE_INDIVIDUAL_MAPS
        % Calculate displacement fields (guaranteed to be periodic)
        Disp_X = Map_X - X;  % Periodic X displacement
        Disp_V = Map_V - V;  % Periodic V displacement
        
        % Verify periodicity (should be very close to zero)
        period_error_x = max(abs(Disp_X(:,end) - Disp_X(:,1)));  % X periodicity
        period_error_v_left = max(abs(Disp_V(1,:) - Disp_V(end,:)));  % V periodicity (left-right)
        
        figure(i);
        sgtitle(sprintf('Map %d: PERIODIC Displacement Fields (amp1=%.2f, nx=%d, amp2=%.2f, nx=%d)', ...
                       i, amp1, freq1_x, amp2, freq2_x));
        
        subplot(1, 2, 1);
        pcolor(X, V, Disp_X); shading flat; colorbar;
        title(sprintf('ΔX (periodic, max=%.3f, period_err=%.2e)', max(abs(Disp_X),[],'all'), period_error_x));
        xlabel('x'); ylabel('v');
        
        subplot(1, 2, 2);
        pcolor(X, V, Disp_V); shading flat; colorbar;
        title(sprintf('ΔV (periodic, max=%.3f, period_err=%.2e)', max(abs(Disp_V),[],'all'), period_error_v_left));
        xlabel('x'); ylabel('v');
        
        % Print displacement and periodicity statistics
        fprintf('Map %d: Max |ΔX| = %.3e, Max |ΔV| = %.3e, Period errors: X=%.2e, V=%.2e\n', ...
                i, max(abs(Disp_X),[],'all'), max(abs(Disp_V),[],'all'), period_error_x, period_error_v_left);
    end
end


%% 3. Compose N Maps: Map_final = Map1 o Map2 o ... o Map_N
fprintf('Calculating the exact analytical composition of %d maps...\n', N_MAPS);

% Start with identity (initial grid positions)
X_current = X;
V_current = V;

% Apply maps in reverse order: Map_N first, then Map_{N-1}, ..., Map_1 last
for i = N_MAPS:-1:1
    fprintf('  Applying Map %d analytically...\n', i);
    
    % Get parameters for this map
    amp1 = map_params(i,1);
    freq1_x = map_params(i,2);
    freq1_v = map_params(i,3);
    amp2 = map_params(i,4);
    freq2_x = map_params(i,5);
    func1 = map_params(i,6);
    func2 = map_params(i,7);
    
    % Apply the transformation to current coordinates
    if func1 == 1  % sin-sin
        X_new = X_current + amp1 * sin(freq1_x * X_current) .* sin(freq1_v * V_current);
    else  % cos-cos
        X_new = X_current + amp1 * cos(freq1_x * X_current) .* cos(freq1_v * V_current);
    end
    
    if func2 == 1  % sin
        V_new = V_current + amp2 * sin(freq2_x * X_current);
    else  % cos
        V_new = V_current + amp2 * cos(freq2_x * X_current);
    end
    
    % Update current coordinates
    X_current = X_new;
    V_current = V_new;
end

% Store the final result
Map_final_exact = zeros(Nv, Nx, 2);
Map_final_exact(:,:,1) = X_current;
Map_final_exact(:,:,2) = V_current;

%% 4. Compose the Maps Numerically for Comparison
fprintf('Calculating the numerical composition of %d maps...\n', N_MAPS);

% Create grid structure (matching CMM conventions)
grid.X = X;
grid.V = V;
grid.x = x;
grid.v = v;
grid.Lx = Lx;
grid.Lv = Lv;
grid.dx = dx;
grid.dv = dv;

% Create params structure
params.Lx = Lx;
params.Lv = Lv;

% Call numerical composition function
addpath('../interpolation');  % Add path to interpolation function
addpath('../interpolation/lagrange');  % Add path to lagrange interpolation
Map_final_numerical = compose_maps_numerical(Maps, grid, params);

%% 5. Compare Analytical vs Numerical Results - SHOWING DISPLACEMENT FIELDS
fprintf('Comparing analytical and numerical results for %d-map composition...\n', N_MAPS);

% Calculate displacement fields for final composed maps
Final_Disp_X_exact = Map_final_exact(:,:,1) - X;
Final_Disp_V_exact = Map_final_exact(:,:,2) - V;
Final_Disp_X_numerical = Map_final_numerical(:,:,1) - X;
Final_Disp_V_numerical = Map_final_numerical(:,:,2) - V;

% Check final periodicity
final_period_error_x_exact = max(abs(Final_Disp_X_exact(:,end) - Final_Disp_X_exact(:,1)));
final_period_error_x_numerical = max(abs(Final_Disp_X_numerical(:,end) - Final_Disp_X_numerical(:,1)));

% Calculate absolute differences
diff_x = abs(Final_Disp_X_exact - Final_Disp_X_numerical);
diff_v = abs(Final_Disp_V_exact - Final_Disp_V_numerical);

% Visualize the comparison - DISPLACEMENT FIELDS
figure(N_MAPS + 1);
sgtitle(sprintf('%d-Map Composition: PERIODIC Displacement Fields (Analytical vs Numerical)', N_MAPS));

subplot(2, 3, 1);
pcolor(X, V, Final_Disp_X_exact); shading flat; colorbar;
title(sprintf('Analytical ΔX (period_err=%.2e)', final_period_error_x_exact));
xlabel('x'); ylabel('v');

subplot(2, 3, 2);
pcolor(X, V, Final_Disp_X_numerical); shading flat; colorbar;
title(sprintf('Numerical ΔX (period_err=%.2e)', final_period_error_x_numerical));
xlabel('x'); ylabel('v');

subplot(2, 3, 3);
pcolor(X, V, diff_x); shading flat; colorbar;
title('|Difference| in ΔX');
xlabel('x'); ylabel('v');

subplot(2, 3, 4);
pcolor(X, V, Final_Disp_V_exact); shading flat; colorbar;
title('Analytical ΔV');
xlabel('x'); ylabel('v');

subplot(2, 3, 5);
pcolor(X, V, Final_Disp_V_numerical); shading flat; colorbar;
title('Numerical ΔV');
xlabel('x'); ylabel('v');

subplot(2, 3, 6);
pcolor(X, V, diff_v); shading flat; colorbar;
title('|Difference| in ΔV');
xlabel('x'); ylabel('v');

%% 6. Calculate Comprehensive Error Norms (L1, L2, L∞) and Relative Errors
fprintf('\n=== COMPREHENSIVE ERROR ANALYSIS ===\n');

% Calculate error norms for displacement differences
error_X = Final_Disp_X_numerical - Final_Disp_X_exact;
error_V = Final_Disp_V_numerical - Final_Disp_V_exact;

% L1 (Manhattan) norm: sum of absolute values
L1_error_X = sum(abs(error_X(:))) * dx * dv;  % Scaled by grid spacing for proper integration
L1_error_V = sum(abs(error_V(:))) * dx * dv;

% L2 (Euclidean) norm: RMS scaled by domain area
L2_error_X = sqrt(sum(error_X(:).^2) * dx * dv);
L2_error_V = sqrt(sum(error_V(:).^2) * dx * dv);

% L∞ (Maximum) norm: maximum absolute value
Linf_error_X = max(abs(error_X(:)));
Linf_error_V = max(abs(error_V(:)));

% Calculate norms of analytical solution for relative errors
L1_analytical_X = sum(abs(Final_Disp_X_exact(:))) * dx * dv;
L1_analytical_V = sum(abs(Final_Disp_V_exact(:))) * dx * dv;

L2_analytical_X = sqrt(sum(Final_Disp_X_exact(:).^2) * dx * dv);
L2_analytical_V = sqrt(sum(Final_Disp_V_exact(:).^2) * dx * dv);

Linf_analytical_X = max(abs(Final_Disp_X_exact(:)));
Linf_analytical_V = max(abs(Final_Disp_V_exact(:)));

% Calculate relative errors (avoid division by zero)
rel_L1_error_X = L1_error_X / max(L1_analytical_X, eps);
rel_L1_error_V = L1_error_V / max(L1_analytical_V, eps);

rel_L2_error_X = L2_error_X / max(L2_analytical_X, eps);
rel_L2_error_V = L2_error_V / max(L2_analytical_V, eps);

rel_Linf_error_X = Linf_error_X / max(Linf_analytical_X, eps);
rel_Linf_error_V = Linf_error_V / max(Linf_analytical_V, eps);



% Print comprehensive error statistics
fprintf('\n--- ABSOLUTE ERRORS ---\n');
fprintf('L1 Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', L1_error_X, L1_error_V);
fprintf('L2 Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', L2_error_X, L2_error_V);
fprintf('L∞ Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', Linf_error_X, Linf_error_V);

fprintf('\n--- RELATIVE ERRORS ---\n');
fprintf('Relative L1 Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', rel_L1_error_X, rel_L1_error_V);
fprintf('Relative L2 Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', rel_L2_error_X, rel_L2_error_V);
fprintf('Relative L∞ Errors:\n');
fprintf('  ΔX: %.6e,  ΔV: %.6e\n', rel_Linf_error_X, rel_Linf_error_V);

fprintf('\n=== ERROR ANALYSIS COMPLETE ===\n');

%% 7. CONVERGENCE ANALYSIS - Error vs Grid Resolution
fprintf('\n=== STARTING CONVERGENCE ANALYSIS ===\n');

% Define grid resolutions to test (powers of 2)
N_resolutions = [2^5, 2^6, 2^7, 2^8];  % 32, 64, 128, 256
n_res = length(N_resolutions);

% Store L∞ errors for each resolution
Linf_errors_X = zeros(n_res, 1);
Linf_errors_V = zeros(n_res, 1);
grid_spacings = zeros(n_res, 1);

% Use same map parameters but test on different grids
fprintf('Testing convergence on %d different grid resolutions...\n', n_res);

for res_idx = 1:n_res
    Nx_test = N_resolutions(res_idx);
    Nv_test = N_resolutions(res_idx);
    
    fprintf('  Resolution %d/%d: %dx%d grid...\n', res_idx, n_res, Nx_test, Nv_test);
    
    % Create test grid with same domain but different resolution
    x_test = (0:Nx_test-1) * Lx / Nx_test;
    v_test = (0:Nv_test-1) * 2*Lv / Nv_test - Lv;
    [X_test, V_test] = meshgrid(x_test, v_test);
    
    dx_test = x_test(2) - x_test(1);
    dv_test = v_test(2) - v_test(1);
    grid_spacings(res_idx) = sqrt(dx_test^2 + dv_test^2);  % Combined grid spacing
    
    % Create test grid structure
    grid_test.X = X_test;
    grid_test.V = V_test;
    grid_test.x = x_test;
    grid_test.v = v_test;
    grid_test.Lx = Lx;
    grid_test.Lv = Lv;
    grid_test.dx = dx_test;
    grid_test.dv = dv_test;
    
    % Generate maps for this resolution
    Maps_test = zeros(Nv_test, Nx_test, 2, N_MAPS);
    
    for i = 1:N_MAPS
        % Get parameters for this map (same as before)
        amp1 = map_params(i,1);
        freq1_x = map_params(i,2);
        freq1_v = map_params(i,3);
        amp2 = map_params(i,4);
        freq2_x = map_params(i,5);
        func1 = map_params(i,6);
        func2 = map_params(i,7);
        
        % Apply function choices
        if func1 == 1  % sin-sin
            Map_X_test = X_test + amp1 * sin(freq1_x * X_test) .* sin(freq1_v * V_test);
        else  % cos-cos
            Map_X_test = X_test + amp1 * cos(freq1_x * X_test) .* cos(freq1_v * V_test);
        end
        
        if func2 == 1  % sin
            Map_V_test = V_test + amp2 * sin(freq2_x * X_test);
        else  % cos
            Map_V_test = V_test + amp2 * cos(freq2_x * X_test);
        end
        
        Maps_test(:,:,1,i) = Map_X_test;
        Maps_test(:,:,2,i) = Map_V_test;
    end
    
    % Analytical composition for this resolution
    X_current_test = X_test;
    V_current_test = V_test;
    
    for i = N_MAPS:-1:1
        amp1 = map_params(i,1);
        freq1_x = map_params(i,2);
        freq1_v = map_params(i,3);
        amp2 = map_params(i,4);
        freq2_x = map_params(i,5);
        func1 = map_params(i,6);
        func2 = map_params(i,7);
        
        if func1 == 1  % sin-sin
            X_new_test = X_current_test + amp1 * sin(freq1_x * X_current_test) .* sin(freq1_v * V_current_test);
        else  % cos-cos
            X_new_test = X_current_test + amp1 * cos(freq1_x * X_current_test) .* cos(freq1_v * V_current_test);
        end
        
        if func2 == 1  % sin
            V_new_test = V_current_test + amp2 * sin(freq2_x * X_current_test);
        else  % cos
            V_new_test = V_current_test + amp2 * cos(freq2_x * X_current_test);
        end
        
        X_current_test = X_new_test;
        V_current_test = V_new_test;
    end
    
    Map_final_exact_test = zeros(Nv_test, Nx_test, 2);
    Map_final_exact_test(:,:,1) = X_current_test;
    Map_final_exact_test(:,:,2) = V_current_test;
    
    % Numerical composition for this resolution
    Map_final_numerical_test = compose_maps_numerical(Maps_test, grid_test, params);
    
    % Calculate errors
    Final_Disp_X_exact_test = Map_final_exact_test(:,:,1) - X_test;
    Final_Disp_V_exact_test = Map_final_exact_test(:,:,2) - V_test;
    Final_Disp_X_numerical_test = Map_final_numerical_test(:,:,1) - X_test;
    Final_Disp_V_numerical_test = Map_final_numerical_test(:,:,2) - V_test;
    
    error_X_test = Final_Disp_X_numerical_test - Final_Disp_X_exact_test;
    error_V_test = Final_Disp_V_numerical_test - Final_Disp_V_exact_test;
    
    % Store L∞ error norms
    Linf_errors_X(res_idx) = max(abs(error_X_test(:)));
    Linf_errors_V(res_idx) = max(abs(error_V_test(:)));
    
    fprintf('    Grid: %dx%d, L∞ errors: ΔX=%.3e, ΔV=%.3e\n', ...
            Nx_test, Nv_test, Linf_errors_X(res_idx), Linf_errors_V(res_idx));
end

% Calculate convergence rates
fprintf('\nCalculating convergence rates...\n');

% Fit power law: error = C * h^p, where h is grid spacing, p is convergence rate
% log(error) = log(C) + p * log(h)
log_h = log(grid_spacings);

% Calculate L∞ convergence rates
conv_rates = struct();

% L∞ convergence rates
p_Linf_X = polyfit(log_h, log(Linf_errors_X), 1);
p_Linf_V = polyfit(log_h, log(Linf_errors_V), 1);
conv_rates.Linf_X = p_Linf_X(1);
conv_rates.Linf_V = p_Linf_V(1);

% Create convergence plot
figure(N_MAPS + 2);
set(gcf, 'Position', [100, 100, 600, 500]);
sgtitle(sprintf('%d-Map Composition: L∞ Convergence Analysis', N_MAPS), 'FontSize', 14);

% L∞ norm convergence vs number of grid points
loglog(N_resolutions, Linf_errors_X, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(N_resolutions, Linf_errors_V, 'bs-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Grid Points'); ylabel('L∞ Error');
title('L∞ Norm Convergence');
legend('ΔX', 'ΔV', 'Location', 'best');
grid on;

% Print convergence results
fprintf('\n=== CONVERGENCE ANALYSIS RESULTS ===\n');
fprintf('Grid resolutions tested: ');
for i = 1:n_res
    fprintf('%dx%d ', N_resolutions(i), N_resolutions(i));
end
fprintf('\n\n');

fprintf('Measured L∞ Convergence Rates:\n');
fprintf('  ΔX: %6.2f\n', conv_rates.Linf_X);
fprintf('  ΔV: %6.2f\n', conv_rates.Linf_V);
fprintf('\n');

% Interpretation
fprintf('Interpretation:\n');
avg_rate = mean([conv_rates.Linf_X, conv_rates.Linf_V]);
if avg_rate > 3.5
    fprintf('✓ Excellent convergence observed (close to O(h⁴) theoretical)\n');
elseif avg_rate > 2.5
    fprintf('✓ Good convergence observed\n');
elseif avg_rate > 1.5
    fprintf('○ Moderate convergence observed\n');
else
    fprintf('⚠ Slower than expected convergence - may indicate:\n');
    fprintf('  - Interpolation limitations\n');
    fprintf('  - Map complexity effects\n');
    fprintf('  - Accumulated composition errors\n');
end

fprintf('\n=== CONVERGENCE ANALYSIS COMPLETE ===\n');