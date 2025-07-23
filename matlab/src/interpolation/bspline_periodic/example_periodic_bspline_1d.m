% Example script to test 1D periodic B-spline interpolation

n = 32;
x = linspace(0, 2*pi, n+1); x(end) = [];
fun = @(x) sin(x) + 0.5*cos(2*x); % Example periodic data
u = fun(x);
degree = 3; % Cubic spline

% Construct the periodic B-spline interpolant
interp = periodic_bspline_interpolant_1d(u, degree);

% Query points for evaluation
xq = linspace(0, 2*pi, 200);
vq = interp.evaluate(xq);

% Plot the results
figure;
plot(x, u, 'o', 'DisplayName', 'Data nodes'); hold on;
plot(xq, vq, '-', 'LineWidth', 1.5, 'DisplayName', 'Periodic B-spline Interpolant');
plot(xq, fun(xq), '--', 'LineWidth', 1.5, 'DisplayName', 'True function');

legend;
xlabel('x'); ylabel('u(x)');
title('1D Periodic B-spline Interpolation Example');
grid on;