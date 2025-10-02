function [Map_composed_numeric] = compose_maps_numerical(Maps, grid, params)
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

% Get number of maps
N_maps = params.Nmaps;

x = grid.x;
v = grid.v; 

% Initialize with identity
current_numeric_x = grid.X;
current_numeric_v = grid.V;

% Apply each map in reverse order using interpolation
for i = N_maps:-1:1
    
    % Get current map's displacement field
    Delta_X = Maps(:,:,1,i) - grid.X;
    Delta_V = Maps(:,:,2,i) - grid.V;
    
    % Current query points  
    Query_X = current_numeric_x;
    Query_V = current_numeric_v;
 
    
    %  Coordinate transformation to ensure interpolation function receives positive values

    Query_V_shifted = Query_V + grid.Lv;
    v_shifted = v + grid.Lv;
    
    %  Use correct coordinates to call interpolation function
    twopi = 2*pi;
    Interpolated_Delta_X = interp2d_periodic(Query_X/params.Lx*twopi, Query_V_shifted/params.Lv*pi, x/params.Lx*twopi, v_shifted/params.Lv*pi, Delta_X); % all coordinates must be rescaled to [0,2pi]
    Interpolated_Delta_V = interp2d_periodic(Query_X/params.Lx*twopi, Query_V_shifted/params.Lv*pi, x/params.Lx*twopi, v_shifted/params.Lv*pi, Delta_V);
    
    % Reshape interpolated results
    Interpolated_Delta_X = reshape(Interpolated_Delta_X, size(Query_X));
    Interpolated_Delta_V = reshape(Interpolated_Delta_V, size(Query_V));
    
    % Apply transformation: new_point = old_point + displacement
    current_numeric_x = Query_X + Interpolated_Delta_X;
    current_numeric_v = Query_V + Interpolated_Delta_V;
end

Map_composed_numeric(:,:,1) = current_numeric_x;
Map_composed_numeric(:,:,2) = current_numeric_v;

end 