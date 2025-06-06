function [Map_composed_numeric] = compose_maps_numerical(Maps, grid, Params)
% COMPOSE_MAPS_NUMERICAL Numerically compose a sequence of maps using displacement interpolation
%
% SYNTAX:
%   Map_composed = compose_maps_numerical(Maps, grid, Params)
%
% INPUTS:
%   Maps            - 4D array [Nv, Nx, 2, N_maps] containing the maps
%   grid            - Structure containing grid information:
%                     grid.X, grid.V - 2D meshgrid arrays for coordinates
%                     grid.Lx, grid.Lv - Domain sizes
%                     grid.dx, grid.dv - Grid spacings
%   Params          - Structure containing parameters

% OUTPUT:
%   Map_composed_numeric - 3D array [Nv, Nx, 2] containing the composed map


order = 4;

% Get number of maps
N_maps = size(Maps, 4);

x = grid.x;
v = grid.v; 
% Initialize with identity
current_numeric_x = grid.X;
current_numeric_v = grid.V;

% Apply each map in reverse order using interpolation
for i = N_maps:-1:1
    fprintf('Applying Map%d numerically...\n', i);
    
    % Get current map's displacement field
    Delta_X = Maps(:,:,1,i) - grid.X;
    Delta_V = Maps(:,:,2,i) - grid.V;
    
    % Current query points  
    Query_X = current_numeric_x;
    Query_V = current_numeric_v;
    

    
    % Interpolate displacement at query points
    Interpolated_Delta_X = lagrange2d_local_interp_periodic(Query_X, Query_V, x, v, Delta_X, order);
    Interpolated_Delta_V = lagrange2d_local_interp_periodic(Query_X, Query_V, x, v, Delta_V, order);
    Interpolated_Delta_X = reshape(Interpolated_Delta_X, size(Query_X));
    Interpolated_Delta_V = reshape(Interpolated_Delta_V, size(Query_V));
    
    % Apply transformation: new_point = old_point + displacement
    current_numeric_x = Query_X + Interpolated_Delta_X;
    current_numeric_v = Query_V + Interpolated_Delta_V;
end

Map_composed_numeric(:,:,1) = current_numeric_x;
Map_composed_numeric(:,:,2) = current_numeric_v;

end 