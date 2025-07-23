
clear; close all; clc;

% Function to interpolate
f = @(x) sin(2*pi*x).*cos(6*pi*x)+1;

% Interpolation order (number of stencil points)
order = 4;

% Grid resolutions to test
Ns = 2.^(4:10);  % e.g., [16, 32, 64, 128, 256, 512]
errors = zeros(size(Ns));

% Reference solution: very fine grid
N_ref = 2^14;
x_ref = ((0:N_ref-1) + 0.5) / N_ref;  % offset points to avoid grid nodes
f_ref = f(x_ref);

for k = 1:length(Ns)
    N = Ns(k);
    xgrid = (0:N-1) / N;           % periodic coarse grid
    ygrid = f(xgrid);              % function values at coarse grid
    f_interp = lagrange_local_interp_periodic(x_ref, xgrid, ygrid, order);
    
    % Error compared to exact solution
    errors(k) = max(abs(f_interp - f_ref));
end

% Estimate convergence rate with linear fit
logN = log10(Ns);
logE = log10(errors);
p = polyfit(logN, logE, 1);
logE_fit = polyval(p, logN);

% Plot convergence
figure;
loglog(Ns, errors, '-o', 'LineWidth', 2); hold on;
loglog(Ns, 10.^logE_fit, '--r', 'LineWidth', 1.5);
xlabel('Number of grid points (N)');
ylabel('Max interpolation error');
title(['Convergence of Local Lagrange Interpolation (order = ', num2str(order), ')']);
legend('Error', ['Fit slope = ', num2str(p(1), '%.2f')], 'Location', 'southwest');
grid on;


