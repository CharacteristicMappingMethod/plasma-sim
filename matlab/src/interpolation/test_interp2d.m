%% Test script for 2D periodic Lagrange interpolation and convergence evaluation

orders = [2, 4];
Ns = 2.^(4:8);
errors = zeros(length(Ns), length(orders));

% Exact function and evaluation grid
f = @(x,y) sin(2*pi*x) .* cos(2*pi*y);
[x_eval, y_eval] = meshgrid(linspace(0,1,200), linspace(0,1,200));
x_eval = x_eval(:);
y_eval = y_eval(:);
f_exact = f(x_eval, y_eval);

for o = 1:length(orders)
    order = orders(o);
    for k = 1:length(Ns)
        N = Ns(k);
        xgrid = (0:N-1)/N;
        ygrid = (0:N-1)/N;
        [X, Y] = meshgrid(xgrid, ygrid);
        F = f(X, Y);
        F_interp = lagrange2d_local_interp_periodic(x_eval, y_eval, xgrid, ygrid, F, order);
        errors(k, o) = max(abs(F_interp - f_exact));
    end
end

% Plot convergence
figure;
loglog(Ns, errors, '-o');
hold on;
legend_entries = cell(1, length(orders));
for o = 1:length(orders)
    p = polyfit(log(Ns), log(errors(:,o))', 1);
    loglog(Ns, exp(polyval(p,log(Ns))), '--');
    legend_entries{o} = sprintf('Order %d, slope %.2f', orders(o), p(1));
end
xlabel('N'); ylabel('Max Error');
title('Convergence of 2D Periodic Lagrange Interpolation');
legend(legend_entries); grid on;
