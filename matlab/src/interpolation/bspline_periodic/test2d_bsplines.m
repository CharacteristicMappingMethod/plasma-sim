%% Convergence test for 2D periodic B-spline interpolation (cubic)
% Tests both MATLAB and MEX implementations

Ns = 2.^(4:8);
errors_matlab = zeros(length(Ns), 1);
errors_mex = zeros(length(Ns), 1);
degree_x = 3; degree_y = 3;

% Exact function and evaluation grid
f = @(x,y) sin(x) .* cos(y) + 0.5 * cos(2*x - y);
[x_eval, y_eval] = meshgrid(linspace(0,2*pi,200), linspace(0,2*pi,200));
x_eval = x_eval(:);
y_eval = y_eval(:);
f_exact = f(y_eval, x_eval);

% Check if MEX function is available
mex_available = exist('bspline_periodic_eval_2d_mex_cpp', 'file') == 3;
if mex_available
    fprintf('MEX implementation available - testing both MATLAB and MEX\n');
else
    fprintf('MEX implementation not available - testing MATLAB only\n');
    fprintf('Run compile_mex to build the MEX function\n');
end

for k = 1:length(Ns)
    N = Ns(k);
    dx = 2*pi / N;
    dy = 2*pi / N;
    xgrid = linspace(0, 2*pi - dx, N);
    ygrid = linspace(0, 2*pi - dy, N);
    [X, Y] = meshgrid(xgrid, ygrid);
    F = f(X, Y);

    % Create interpolant using MATLAB function
    interp_matlab = periodic_bspline_interpolant_2d(F, degree_x, degree_y, 'matlab');
    F_interp_matlab = interp_matlab.evaluate(x_eval, y_eval);
    errors_matlab(k) = max(abs(F_interp_matlab - f_exact));

    % Create interpolant using MEX function if available
    if mex_available
        interp_mex = periodic_bspline_interpolant_2d(F, degree_x, degree_y, 'mex');
        F_interp_mex = interp_mex.evaluate(x_eval, y_eval);
        errors_mex(k) = max(abs(F_interp_mex - f_exact));

        % Check agreement between MATLAB and MEX
        max_diff = max(abs(F_interp_matlab - F_interp_mex));
        if max_diff > 1e-12
            fprintf('Warning: MATLAB and MEX results differ by %e for N=%d\n', max_diff, N);
        end
    end
end

figure(99);
subplot(1,2,1);
pcolor(reshape(F_interp_mex,200,200));
title('MEX');
shading flat;
colorbar;
subplot(1,2,2);
pcolor(reshape(F_interp_matlab,200,200));
title('MATLAB');
shading flat;
colorbar;


% Plot convergence
figure;
loglog(Ns, errors_matlab, '-o', 'LineWidth', 2, 'DisplayName', 'MATLAB');
hold on;
p_matlab = polyfit(log(Ns), log(errors_matlab)', 1);
loglog(Ns, exp(polyval(p_matlab,log(Ns))), '--', 'LineWidth', 1.5, 'DisplayName', sprintf('MATLAB fit, slope %.2f', p_matlab(1)));

if mex_available
    loglog(Ns, errors_mex, '-s', 'LineWidth', 2, 'DisplayName', 'MEX');
    p_mex = polyfit(log(Ns), log(errors_mex)', 1);
    loglog(Ns, exp(polyval(p_mex,log(Ns))), ':', 'LineWidth', 1.5, 'DisplayName', sprintf('MEX fit, slope %.2f', p_mex(1)));
end

xlabel('N'); ylabel('Max Error');
title('Convergence of 2D Periodic Cubic B-spline Interpolation');
legend('Location', 'southwest');
grid on;

% Performance comparison (if MEX is available)
if mex_available
    fprintf('\nPerformance comparison:\n');
    fprintf('Testing with N=64, 1000 evaluation points...\n');

    N_test = 64;
    dx_test = 2*pi / N_test;
    dy_test = 2*pi / N_test;
    xgrid_test = linspace(0, 2*pi - dx_test, N_test);
    ygrid_test = linspace(0, 2*pi - dy_test, N_test);
    [X_test, Y_test] = meshgrid(xgrid_test, ygrid_test);
    F_test = f(X_test, Y_test);
    interp_matlab_test = periodic_bspline_interpolant_2d(F_test, degree_x, degree_y, 'matlab');
    interp_mex_test = periodic_bspline_interpolant_2d(F_test, degree_x, degree_y, 'mex');

    % Create test evaluation points
    n_eval = 1000;
    x_eval_test = 2*pi * rand(n_eval, 1);
    y_eval_test = 2*pi * rand(n_eval, 1);

    % Time MATLAB implementation
    tic;
    for i = 1:10
        F_interp_matlab_test = interp_matlab_test.evaluate(x_eval_test, y_eval_test);
    end
    time_matlab = toc / 10;

    % Time MEX implementation
    tic;
    for i = 1:10
        F_interp_mex_test = interp_mex_test.evaluate(x_eval_test, y_eval_test);
    end
    time_mex = toc / 10;

    fprintf('MATLAB time: %.4f seconds\n', time_matlab);
    fprintf('MEX time: %.4f seconds\n', time_mex);
    fprintf('Speedup: %.2fx\n', time_matlab / time_mex);

    % Verify results match
    max_diff_test = max(abs(F_interp_matlab_test - F_interp_mex_test));
    if max_diff_test < 1e-12
        fprintf('✓ MATLAB and MEX results match (max diff: %e)\n', max_diff_test);
    else
        fprintf('✗ MATLAB and MEX results differ (max diff: %e)\n', max_diff_test);
    end
end