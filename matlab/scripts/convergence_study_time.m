%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Step Convergence Study Script
% Measures how L∞ error decreases as time step decreases
% Time steps: 2^-7, 2^-6, 2^-5, 2^-4 (reference dt = 2^-8)
% Grid fixed at 1024x1024, measured at time t = 10
% CMM method with N_remap = 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('../src/'),genpath('../params/'))
DEFAULTS

PARAMS_two_stream;

% Configuration
dt_values = [2^-9,2^-8,2^-7, 2^-6];  % Time steps to test
dt_ref = 2^-10;                    % Reference time step
grid_size = 2^7;                       % Fixed grid size
target_time = 0.5;
N_remap_value = 5;
params.Tend         = target_time;                 % Simulation end time
params.dt_save      = target_time;  
params.Nsample = [grid_size, grid_size];                % number of grid points in the sample grid
params.Nmap = [grid_size, grid_size];  
params.method = "NuFi";

%% compute reference solution

params.dt = dt_ref;
[params_ref, data_ref] = Sim(params);


%% compute comparisons
i = 1;
for dt = dt_values
    params.dt = dt;
    [params_comp{i}, data_comp{i}] = Sim(params);
    i = i +1; 
end

%% Extract final time distribution functions
% Reference solution (finest time step)
f_ref = data_ref.fs(:,:,end,1);  % Final time, first species

% Comparison solutions
f_comp = cell(length(dt_values), 1);
for i = 1:length(dt_values)
    f_comp{i} = data_comp{i}.fs(:,:,end,1);  % Final time, first species
end

%% Calculate L_infinity norm errors
L_inf_errors = zeros(length(dt_values), 1);

fprintf('\n=== Error Analysis ===\n');
fprintf('%-15s %-20s\n', 'dt', 'L∞ Error');
fprintf('%-15s %-20s\n', '---', '--------');

for i = 1:length(dt_values)
    % Calculate absolute error
    error_abs = abs(f_comp{i} - f_ref);
    
    % L_infinity norm: maximum absolute error
    L_inf_errors(i) = max(error_abs(:));
    
    fprintf('%-15.6e %-20.6e\n', dt_values(i), L_inf_errors(i));
end

%% Plot error vs dt
figure('Name', 'Time Step Convergence Study', 'Position', [100, 100, 800, 600]);

% Log-log plot for convergence study
loglog(dt_values, L_inf_errors, 'o-', 'LineWidth', 2, 'MarkerSize', 10, ...
       'MarkerFaceColor', 'auto', 'MarkerEdgeColor', 'auto');
hold on;

% Optional: Add reference line showing expected convergence rate
% For first-order methods, error ~ dt, so we can fit a line
if all(L_inf_errors > 0)
    % Fit log(error) = log(C) + p*log(dt)
    log_dt = log(dt_values);
    log_error = log(L_inf_errors);
    coeffs = polyfit(log_dt, log_error, 1);
    convergence_rate = coeffs(1);
    log_C = coeffs(2);
    C_constant = exp(log_C);
    
    % Plot fitted line
    dt_fit = linspace(min(dt_values), max(dt_values), 100);
    error_fit = C_constant * (dt_fit).^convergence_rate;
    loglog(dt_fit, error_fit, '--r', 'LineWidth', 1.5);
    
    % Add legend with convergence rate
    legend(sprintf('Computed errors (rate = %.2f)', convergence_rate), ...
           'Location', 'best', 'FontSize', 12);
    
    fprintf('\nConvergence rate: %.2f\n', convergence_rate);
    fprintf('Fitted model: Error = %.3e * dt^%.2f\n', C_constant, convergence_rate);
else
    legend('Computed errors', 'Location', 'best', 'FontSize', 12);
end

hold off;

% Formatting
xlabel('Time Step (dt)', 'FontSize', 14);
ylabel('L∞ Error', 'FontSize', 14);
title('Time Step Convergence Study: L∞ Error vs dt', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

% Set axis limits for better visualization
xlim([min(dt_values)*0.8, max(dt_values)*1.2]);
if all(L_inf_errors > 0)
    ylim([min(L_inf_errors)*0.5, max(L_inf_errors)*2]);
end

fprintf('\n=== Convergence Study Complete ===\n');