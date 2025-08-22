%% Convergence Study for 2D Periodic Interpolation Methods
% This script copies the test from test_interp2d_periodic.m and calculates
% convergence rates for different interpolation schemes
clc
clear all
close all

% Add paths
addpath(genpath('./src/'))

%% Parameters copied from test_interp2d_periodic.m
Ns = 2.^(4:8);
orders = [3, 10]; % Lagrange orders
schemes = {'bspline', 'lagrange-bary', 'lagrange-normal'};
colors = {'-o', '-^', '-s'};

errors = zeros(length(Ns), length(schemes), length(orders));
times = zeros(length(Ns), length(schemes), length(orders));
convergence_rates = zeros(length(schemes), length(orders));

% Test function and evaluation grid
f = @(x,y) sin(x) .* cos(y) + 0.5 * cos(2*x - y);
[x_eval, y_eval] = meshgrid(linspace(0,2*pi,200), linspace(0,2*pi,200));
x_eval = x_eval(:); y_eval = y_eval(:);
f_exact = f(x_eval, y_eval);

fprintf('Starting convergence study...\n');

%% Main convergence study loop (copied from test_interp2d_periodic.m)
for k = 1:length(Ns)
    N = Ns(k);
    dx = 2*pi / N;
    dy = 2*pi / N;
    xgrid = linspace(0, 2*pi - dx, N);
    ygrid = linspace(0, 2*pi - dy, N);
    [X, Y] = meshgrid(xgrid, ygrid);
    F = f(X, Y);

    fprintf('N = %d: ', N);

    % --- B-spline ---
    opts = struct('scheme', 'bspline', 'degree', 3, 'use_mex', true);
    tic;
    F_interp_bspl = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
    times(k,1,1) = toc;
    errors(k,1,1) = max(abs(F_interp_bspl(:) - f_exact(:)));
    fprintf('B-spline done, ');

    % --- Lagrange barycentric ---
    for o = 1:length(orders)
        opts = struct('scheme', 'lagrange-bary', 'order', orders(o), 'use_mex', true);
        tic;
        F_interp_lagr_bary = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
        times(k,2,o) = toc;
        errors(k,2,o) = max(abs(F_interp_lagr_bary(:) - f_exact(:)));
    end
    fprintf('Lagrange-bary done, ');

    % --- Lagrange normal ---
    for o = 1:length(orders)
        opts = struct('scheme', 'lagrange-normal', 'order', orders(o), 'use_mex', true);
        tic;
        F_interp_lagr_norm = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
        times(k,3,o) = toc;
        errors(k,3,o) = max(abs(F_interp_lagr_norm(:) - f_exact(:)));
    end
    fprintf('Lagrange-normal done\n');
end

%% Calculate convergence rates
fprintf('\n=== Convergence Rate Calculation ===\n');

% Use least squares fit to log(error) vs log(h) to find slope
h_values = 2*pi ./ Ns;  % Grid spacing

for s = 1:length(schemes)
    if s == 1
        % B-spline only has one configuration
        valid_errors = errors(:,s,1);
        valid_errors = valid_errors(valid_errors > 0);  % Remove zeros
        if length(valid_errors) >= 2
            log_h = log(h_values(1:length(valid_errors)));
            log_err = log(valid_errors);
            p = polyfit(log_h, log_err, 1);
            convergence_rates(s,1) = -p(1);  % Negative because error decreases with h
            fprintf('B-spline (degree 3): convergence rate = %.2f\n', convergence_rates(s,1));
        end
    else
        % Lagrange methods have multiple orders
        for o = 1:length(orders)
            valid_errors = errors(:,s,o);
            valid_errors = valid_errors(valid_errors > 0);  % Remove zeros
            if length(valid_errors) >= 2
                log_h = log(h_values(1:length(valid_errors)));
                log_err = log(valid_errors);
                p = polyfit(log_h, log_err, 1);
                convergence_rates(s,o) = -p(1);  % Negative because error decreases with h
                fprintf('%s (order %d): convergence rate = %.2f\n', schemes{s}, orders(o), convergence_rates(s,o));
            end
        end
    end
end

%% Plot accuracy (convergence) - copied and enhanced
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1);
for s = 1:length(schemes)
    if s == 1
        loglog(Ns, errors(:,s,1), colors{s}, 'LineWidth', 2, 'MarkerSize', 8, ...
               'DisplayName', sprintf('B-spline (deg 3), rate=%.2f', convergence_rates(s,1))); hold on;
    else
        for o = 1:length(orders)
            loglog(Ns, errors(:,s,o), colors{s}, 'LineWidth', 2, 'MarkerSize', 8, ...
                   'DisplayName', sprintf('%s (ord %d), rate=%.2f', schemes{s}, orders(o), convergence_rates(s,o)));
        end
    end
end

% Add theoretical convergence lines
h_ref = 2*pi ./ Ns;
for rate = [2, 4, 6, 8, 10]
    ref_line = 1e-2 * (h_ref/h_ref(1)).^rate;
    loglog(Ns, ref_line, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(Ns(end), ref_line(end), sprintf('O(h^{%d})', rate), 'FontSize', 10);
end

xlabel('Grid Points N'); 
ylabel('Max Error');
title('Convergence of 2D Periodic Interpolation Methods');
legend('Location', 'southwest', 'FontSize', 10);
grid on;

%% Plot performance
subplot(2,1,2);
for s = 1:length(schemes)
    if s == 1
        loglog(Ns, times(:,s,1), colors{s}, 'LineWidth', 2, 'MarkerSize', 8, ...
               'DisplayName', 'B-spline (deg 3)'); hold on;
    else
        for o = 1:length(orders)
            loglog(Ns, times(:,s,o), colors{s}, 'LineWidth', 2, 'MarkerSize', 8, ...
                   'DisplayName', sprintf('%s (order %d)', schemes{s}, orders(o)));
        end
    end
end
xlabel('Grid Points N'); 
ylabel('Evaluation Time (s)');
title('Performance of 2D Periodic Interpolation Methods');
legend('Location', 'northwest', 'FontSize', 10);
grid on;

%% Save results
results = struct();
results.Ns = Ns;
results.errors = errors;
results.times = times;
results.convergence_rates = convergence_rates;
results.schemes = schemes;
results.orders = orders;

save('interpolation_convergence_results.mat', 'results');

%% Summary table
fprintf('\n=== Summary Table ===\n');
fprintf('Method\t\t\tOrder\tConv. Rate\tMin Error\tMax Time (s)\n');
fprintf('--------------------------------------------------------\n');

% B-spline
fprintf('B-spline\t\t3\t%.2f\t\t%.2e\t%.4f\n', ...
        convergence_rates(1,1), min(errors(:,1,1)), max(times(:,1,1)));

% Lagrange methods
for s = 2:length(schemes)
    for o = 1:length(orders)
        fprintf('%-15s\t%d\t%.2f\t\t%.2e\t%.4f\n', ...
                schemes{s}, orders(o), convergence_rates(s,o), ...
                min(errors(:,s,o)), max(times(:,s,o)));
    end
end

fprintf('\nResults saved to: interpolation_convergence_results.mat\n');
fprintf('Convergence study completed!\n');