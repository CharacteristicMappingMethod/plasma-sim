function [X_eval, V_eval] = evaluate_map(Maps, grid, params, X_query, V_query)
% EVALUATE_MAP Evaluate a sequence of maps at given coordinate points
%
% SYNTAX:
%   [X_eval, V_eval] = evaluate_map(Maps, grid, params, X_query, V_query)
%
% INPUTS:
%   Maps            - 4D array [Nv, Nx, 2, N_maps] containing the maps
%   grid            - Structure containing grid information:
%                     grid.X, grid.V - 2D meshgrid arrays for coordinates
%                     grid.Lx, grid.Lv - Domain sizes
%   params          - Structure containing parameters
%   X_query, V_query - Query points where to evaluate the map (can be arrays)
%
% OUTPUTS:
%   X_eval, V_eval  - Evaluated coordinates after applying all maps

% Get number of maps
N_maps = params.Nmaps;

% Initialize with query points
current_x = X_query;
current_v = V_query;

% Apply each map in reverse order using interpolation
for i = N_maps:-1:1
    
    % Get current map's displacement field
    Delta_X = Maps(:,:,1,i) - grid.X;
    Delta_V = Maps(:,:,2,i) - grid.V;
    
    % Current query points  
    Query_X = current_x;
    Query_V = current_v;
    
    % Coordinate transformation to ensure interpolation function receives positive values
    Query_V_shifted = Query_V + grid.Lv;
    v_shifted = grid.v + grid.Lv;
    
    % Use correct coordinates to call interpolation function
    twopi = 2*pi;
    Interpolated_Delta_X = interp2d_periodic(Query_X/params.Lx*twopi, Query_V_shifted/params.Lv*pi, ...
                                            grid.x/params.Lx*twopi, v_shifted/params.Lv*pi, Delta_X, params.opt_interp);
    Interpolated_Delta_V = interp2d_periodic(Query_X/params.Lx*twopi, Query_V_shifted/params.Lv*pi, ...
                                            grid.x/params.Lx*twopi, v_shifted/params.Lv*pi, Delta_V, params.opt_interp);
    
    % Reshape interpolated results to match input shape
    Interpolated_Delta_X = reshape(Interpolated_Delta_X, size(Query_X));
    Interpolated_Delta_V = reshape(Interpolated_Delta_V, size(Query_V));
    
    % Apply transformation: new_point = old_point + displacement
    current_x = Query_X + Interpolated_Delta_X;
    current_v = Query_V + Interpolated_Delta_V;
end

% Return evaluated coordinates
X_eval = current_x;
V_eval = current_v;

end
