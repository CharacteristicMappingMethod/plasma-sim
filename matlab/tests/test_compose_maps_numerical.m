%% Test script for compose_maps_numerical
clear; close all; clc;

%% Add paths
addpath('../src');
addpath('../');

%% Setup grid and parameters
Nx = 64;
Nv = 64;
Lx = 2*pi;
Lv = 2*pi;

grid.x = linspace(0, Lx, Nx+1); grid.x = grid.x(1:end-1);
grid.v = linspace(0, Lv, Nv+1); grid.v = grid.v(1:end-1);
grid.dx = grid.x(2) - grid.x(1);
grid.dv = grid.v(2) - grid.v(1);
[grid.X, grid.V] = meshgrid(grid.x, grid.v);
grid.Lx = Lx;
grid.Lv = Lv;

params.Nmaps = 2;

%% Define two trigonometric maps
% Amplitudes for the displacement
a = 0.1; 
b = 0.1;
c = 0.1;
d = 0.1;
k = 1; % Wavenumber

% Preallocate Maps array
Maps = zeros(Nv, Nx, 2, params.Nmaps);

% Map 1: (x,v) -> (x + a*sin(k*v), v + b*sin(k*x))
Maps(:,:,1,1) = grid.X + a * sin(k * grid.V);
Maps(:,:,2,1) = grid.V + b * sin(k * grid.X);

% Map 2: (x,v) -> (x + c*sin(k*v), v + d*sin(k*x))
Maps(:,:,1,2) = grid.X + c * sin(k * grid.V);
Maps(:,:,2,2) = grid.V + d * sin(k * grid.X);


%% Compute analytical composition (Map1 o Map2)
% The function applies maps in reverse order, so we compute Map1(Map2(x,v))
% Map2(x,v) = (x', v')
x_prime = grid.X + c * sin(k * grid.V);
v_prime = grid.V + d * sin(k * grid.X);

% Map1(x',v') = (x'', v'')
x_double_prime = x_prime + a * sin(k * v_prime);
v_double_prime = v_prime + b * sin(k * x_prime);

Map_composed_analytical = zeros(Nv, Nx, 2);
Map_composed_analytical(:,:,1) = x_double_prime;
Map_composed_analytical(:,:,2) = v_double_prime;


%% Compute numerical composition
[Map_composed_numeric] = compose_maps_numerical(Maps, grid, params);

%% Compare analytical and numerical results
error_x = max(abs(Map_composed_analytical(:,:,1) - Map_composed_numeric(:,:,1)), [], 'all');
error_v = max(abs(Map_composed_analytical(:,:,2) - Map_composed_numeric(:,:,2)), [], 'all');

fprintf('Max error in X component: %e\n', error_x);
fprintf('Max error in V component: %e\n', error_v);

% Define a tolerance
tolerance = 1e-4; 

% Assert that the error is within tolerance
assert(error_x < tolerance, 'Numerical composition for X-component is not accurate enough.');
assert(error_v < tolerance, 'Numerical composition for V-component is not accurate enough.');

fprintf('Test passed! The numerical composition is close to the analytical one.\n');

%% Plotting results for visual comparison
figure;
subplot(2, 2, 1);
imagesc(grid.x, grid.v, Map_composed_analytical(:,:,1));
title('Analytical X component');
xlabel('x'); ylabel('v');
colorbar;
axis xy;

subplot(2, 2, 2);
imagesc(grid.x, grid.v, Map_composed_numeric(:,:,1));
title('Numerical X component');
xlabel('x'); ylabel('v');
colorbar;
axis xy;

subplot(2, 2, 3);
imagesc(grid.x, grid.v, Map_composed_analytical(:,:,2));
title('Analytical V component');
xlabel('x'); ylabel('v');
colorbar;
axis xy;

subplot(2, 2, 4);
imagesc(grid.x, grid.v, Map_composed_numeric(:,:,2));
title('Numerical V component');
xlabel('x'); ylabel('v');
colorbar;
axis xy;

figure;
subplot(1, 2, 1);
imagesc(grid.x, grid.v, abs(Map_composed_analytical(:,:,1) - Map_composed_numeric(:,:,1)));
title('Error in X component');
xlabel('x'); ylabel('v');
colorbar;
axis xy;

subplot(1, 2, 2);
imagesc(grid.x, grid.v, abs(Map_composed_analytical(:,:,2) - Map_composed_numeric(:,:,2)));
title('Error in V component');
xlabel('x'); ylabel('v');
colorbar;
axis xy; 