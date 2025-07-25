%% Test and compare interp2d_periodic for B-spline, Lagrange barycentric, and Lagrange normal
clc
clear all
close all

Ns = 2.^(4:8);
orders = [2, 4]; % Lagrange orders
schemes = {'bspline', 'lagrange-bary', 'lagrange-normal'};
colors = {'-o', '-^', '-s'};

errors = zeros(length(Ns), length(schemes), length(orders));
times = zeros(length(Ns), length(schemes), length(orders));

% Test function and evaluation grid
f = @(x,y) sin(x) .* cos(y) + 0.5 * cos(2*x - y);
[x_eval, y_eval] = meshgrid(linspace(0,2*pi,200), linspace(0,2*pi,200));
x_eval = x_eval(:); y_eval = y_eval(:);
f_exact = f(x_eval, y_eval);

for k = 1:length(Ns)
    N = Ns(k);
    dx = 2*pi / N;
    dy = 2*pi / N;
    xgrid = linspace(0, 2*pi - dx, N);
    ygrid = linspace(0, 2*pi - dy, N);
    [X, Y] = meshgrid(xgrid, ygrid);
    F = f(X, Y);

    % --- B-spline ---
    opts = struct('scheme', 'bspline', 'degree', 3, 'use_mex', true);
    tic;
    F_interp_bspl = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
    times(k,1,1) = toc;
    errors(k,1,1) = max(abs(F_interp_bspl(:) - f_exact(:)));

    % --- Lagrange barycentric ---
    for o = 1:length(orders)
        opts = struct('scheme', 'lagrange-bary', 'order', orders(o), 'use_mex', true);
        tic;
        F_interp_lagr_bary = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
        times(k,2,o) = toc;
        errors(k,2,o) = max(abs(F_interp_lagr_bary(:) - f_exact(:)));
    end

    % --- Lagrange normal ---
    for o = 1:length(orders)
        opts = struct('scheme', 'lagrange-normal', 'order', orders(o), 'use_mex', true);
        tic;
        F_interp_lagr_norm = interp2d_periodic(x_eval, y_eval, xgrid, ygrid, F, opts);
        times(k,3,o) = toc;
        errors(k,3,o) = max(abs(F_interp_lagr_norm(:) - f_exact(:)));
    end
end

%% Plot accuracy (convergence)
figure;
for s = 1:length(schemes)
    if s == 1
        loglog(Ns, errors(:,s,1), colors{s}, 'LineWidth', 2, 'DisplayName', 'B-spline (deg 3)'); hold on;
    else
        for o = 1:length(orders)
            loglog(Ns, errors(:,s,o), colors{s}, 'LineWidth', 2, 'DisplayName', sprintf('%s (order %d)', schemes{s}, orders(o)));
        end
    end
end
xlabel('N'); ylabel('Max Error');
title('Convergence of 2D Periodic Interpolation Methods');
legend('Location', 'southwest');
grid on;

%% Plot performance
figure;
for s = 1:length(schemes)
    if s == 1
        loglog(Ns, times(:,s,1), colors{s}, 'LineWidth', 2, 'DisplayName', 'B-spline (deg 3)'); hold on;
    else
        for o = 1:length(orders)
            loglog(Ns, times(:,s,o), colors{s}, 'LineWidth', 2, 'DisplayName', sprintf('%s (order %d)', schemes{s}, orders(o)));
        end
    end
end
xlabel('N'); ylabel('Evaluation Time (s)');
title('Performance of 2D Periodic Interpolation Methods');
legend('Location', 'northwest');
grid on;