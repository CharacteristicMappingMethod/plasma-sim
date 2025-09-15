function test_interp1d_periodic()
% Comprehensive test for 1D periodic interpolation
% Tests both B-spline and Lagrange methods

clear; close all; clc;

fprintf('Testing 1D Periodic Interpolation\n');
fprintf('================================\n\n');

% Test functions
test_functions = {
    % @(x) sin(2*pi*x) + 0.5*cos(4*pi*x), 'sin(2πx) + 0.5cos(4πx)';
     @(x) exp(sin(2*pi*x)), 'exp(sin(2πx))';
     @(x) sin(2*pi*x).*cos(6*pi*x) + 1, 'sin(2πx)cos(6πx) + 1';
    %@(x) sin(8*pi*x) + cos(12*pi*x), 'sin(8πx) + cos(12πx)';
};

% Interpolation schemes to test
schemes = {'bspline', 'lagrange-bary', 'lagrange-normal'};

% Grid resolutions to test
Ns = 2.^(4:10);  % [16, 32, 64, 128, 256]

% Reference solution: very fine grid
N_ref = 2^12;
x_ref = ((0:N_ref-1) + 0.5) / N_ref;  % offset points to avoid grid nodes

% Initialize storage for all results
all_errors = cell(size(test_functions, 1), length(schemes));
all_times = cell(size(test_functions, 1), length(schemes));
all_convergence_rates = zeros(size(test_functions, 1), length(schemes));

% Test each function
for func_idx = 1:size(test_functions, 1)
    func = test_functions{func_idx, 1};
    func_name = test_functions{func_idx, 2};
    
    fprintf('Testing function: %s\n', func_name);
    fprintf('----------------------------------------\n');
    
    % Evaluate reference solution
    f_ref = func(x_ref);
    
    % Test each scheme
    for scheme_idx = 1:length(schemes)
        scheme = schemes{scheme_idx};
        fprintf('\nScheme: %s\n', scheme);
        
        errors = zeros(size(Ns));
        times = zeros(size(Ns));
        
        for k = 1:length(Ns)
            N = Ns(k);
            xgrid = (0:N-1) / N;           % periodic coarse grid
            fgrid = func(xgrid);            % function values at coarse grid
            
            % Set up interpolation options
            opts = struct();
            opts.scheme = scheme;
            opts.order = 3;  % for Lagrange

            opts.degree = 3; % for B-spline
            opts.use_mex = true;
            opts.cache = true;
            
            % Time the interpolation
            tic;
            f_interp = interp1d_periodic(x_ref, xgrid, fgrid, opts);
            times(k) = toc;
            
            % Error compared to exact solution
            errors(k) = max(abs(f_interp - f_ref));
        end
        
        % Store results
        all_errors{func_idx, scheme_idx} = errors;
        all_times{func_idx, scheme_idx} = times;
        
        % Estimate convergence rate
        logN = log10(Ns);
        logE = log10(errors);
        p = polyfit(logN, logE, 1);
        all_convergence_rates(func_idx, scheme_idx) = p(1);
        
        fprintf('  Convergence rate: %.2f\n', p(1));
        fprintf('  Final error (N=%d): %.2e\n', Ns(end), errors(end));
        fprintf('  Average time (N=%d): %.4f ms\n', Ns(end), times(end)*1000);
    end
    
    fprintf('\n');
end

% Create comparison plots - only error convergence and timing
figure('Position', [100, 100, 1200, 500]);

% Colors and markers for different schemes
colors = {'r', 'g', 'b'};
markers = {'o', 's', '^'};
linestyles = {'-', '--', ':', '-.'};

% Plot 1: Error convergence for all functions and schemes
subplot(1, 2, 1);
for func_idx = 1:size(test_functions, 1)
    for scheme_idx = 1:length(schemes)
        errors = all_errors{func_idx, scheme_idx};
        loglog(Ns, errors, 'Color', colors{scheme_idx}, 'Marker', markers{scheme_idx}, ...
               'LineStyle', linestyles{func_idx}, 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
    end
end
xlabel('Number of grid points (N)');
ylabel('Max interpolation error');
title('Error Convergence - All Methods');
legend_entries = {};
for func_idx = 1:size(test_functions, 1)
    for scheme_idx = 1:length(schemes)
        legend_entries{end+1} = sprintf('%s - %s', test_functions{func_idx, 2}, schemes{scheme_idx});
    end
end
legend(legend_entries, 'Location', 'southwest', 'FontSize', 10);
grid on;

% Plot 2: Timing comparison for all functions and schemes
subplot(1, 2, 2);
for func_idx = 1:size(test_functions, 1)
    for scheme_idx = 1:length(schemes)
        times = all_times{func_idx, scheme_idx};
        semilogx(Ns, times*1000, 'Color', colors{scheme_idx}, 'Marker', markers{scheme_idx}, ...
                 'LineStyle', linestyles{func_idx}, 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
    end
end
xlabel('Number of grid points (N)');
ylabel('Interpolation time (ms)');
title('Timing Comparison - All Methods');
legend(legend_entries, 'Location', 'northwest', 'FontSize', 10);
grid on;

% Test periodic boundary conditions
fprintf('Testing periodic boundary conditions...\n');
test_periodic_boundary();

% Test interpolation accuracy at specific points
fprintf('Testing interpolation accuracy...\n');
test_interpolation_accuracy();

fprintf('All tests completed!\n');

end

function test_periodic_boundary()
% Test that interpolation respects periodic boundary conditions

N = 64;
xgrid = (0:N-1) / N;
f = @(x) sin(2*pi*x) + cos(4*pi*x);
fgrid = f(xgrid);

% Test points at boundaries
x_test = [0, 1, 0.5, 0.25, 0.75];
opts = struct('scheme', 'bspline', 'degree', 3);

f_interp = interp1d_periodic(x_test, xgrid, fgrid, opts);
f_exact = f(x_test);

fprintf('  Boundary test errors: %.2e\n', max(abs(f_interp - f_exact)));

% Test that f(0) = f(1) for periodic functions
assert(abs(f_interp(1) - f_interp(2)) < 1e-10, 'Periodic boundary condition failed!');
end

function test_interpolation_accuracy()
% Test interpolation accuracy for known functions

N = 32;
xgrid = (0:N-1) / N;

% Test polynomial that should be exact for cubic B-splines
f = @(x) x.^3 - 2*x.^2 + x + 1;
fgrid = f(xgrid);

x_test = linspace(0, 1, 100);
opts = struct('scheme', 'bspline', 'degree', 3);

f_interp = interp1d_periodic(x_test, xgrid, fgrid, opts);
f_exact = f(x_test);

fprintf('  Cubic polynomial test error: %.2e\n', max(abs(f_interp - f_exact)));
end
