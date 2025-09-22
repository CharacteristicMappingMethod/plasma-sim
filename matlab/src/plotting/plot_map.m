function plot_map(X, V, grid_struct)
% PLOT_MAP Plot coordinate lines of a map transformation
%
% PLOT_MAP(X, V, grid_struct) plots the coordinate lines showing how the original
% grid is transformed by the map (X, V).
%
% Inputs:
%   X, V    - Mapped coordinate arrays [Nv, Nx]
%   grid_struct - Grid structure containing original coordinates
%             Required fields: x, v, X, V, Lx, Lv
%
% Example usage with CMM data:
%   X = params.coordinate_maps.X(:,:,1);
%   V = params.coordinate_maps.V(:,:,1);
%   plot_map(X, V, params.grids(1));

% Validate input dimensions
if ~isequal(size(X), size(V))
    error('X and V must have the same dimensions');
end

if ~isequal(size(X), size(grid_struct.X))
    error('X and V must have the same dimensions as grid_struct.X and grid_struct.V');
end

% Get original grid coordinates
x_orig = grid_struct.x;
v_orig = grid_struct.v;

% Create figure
figure;
hold on;

% Plot original grid lines (dashed)
nlines = 10;
x_indices = round(linspace(1, length(x_orig), nlines));
v_indices = round(linspace(1, length(v_orig), nlines));

% Plot mapped coordinate lines (solid)
for i = 1:length(x_indices)
    idx = x_indices(i);
    plot(X(idx, :), V(idx, :), '-', 'Color', 'b', 'LineWidth', 1.5);
    plot(X(idx, :)-grid_struct.Lx, V(idx, :), '-', 'Color', 'b', 'LineWidth', 1.5);
    plot(X(idx, :)+grid_struct.Lx, V(idx, :), '-', 'Color', 'b', 'LineWidth', 1.5);
end

for i = 1:length(v_indices)
    idx = v_indices(i);
    plot(X(:, idx), V(:, idx), '-', 'Color', 'b', 'LineWidth', 1.5);
    plot(X(:, idx)-grid_struct.Lx, V(:, idx), '-', 'Color', 'b', 'LineWidth', 1.5);
    plot(X(:, idx)+grid_struct.Lx, V(:, idx), '-', 'Color', 'b', 'LineWidth', 1.5);
end

% Set plot properties
xlabel('$x$');
ylabel('$v$');
title('Coordinate Map');
box on;


% Set axis limits
xlim([0, grid_struct.Lx]);
ylim([-grid_struct.Lv, grid_struct.Lv]);


hold off;

end
