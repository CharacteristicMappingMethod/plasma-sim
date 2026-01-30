function [v_refined, original_indices, dv_refined] = refine_velocity_grid(v, refine_factor, v_range)
% REFINE_VELOCITY_GRID Refines a velocity grid by inserting points between existing grid points
%
% Inputs:
%   v - Original velocity grid (1D array)
%   refine_factor - Refinement factor (integer >= 2)
%                   If 2: inserts 1 point between each pair
%                   If 3: inserts 2 points between each pair
%                   If n: inserts (n-1) points between each pair
%   v_range - Velocity range to refine [v_min, v_max] (optional, default: entire grid)
%
% Outputs:
%   v_refined - Refined velocity grid with inserted points
%   original_indices - Indices in v_refined that correspond to original grid points
%                      such that v_refined(original_indices) == v
%   dv_refined - Velocity grid spacing for refined grid (1D array)
%
% Example:
%   v = [-2, -1, 0, 1, 2];
%   [v_refined, indices] = refine_velocity_grid(v, 2);  % Refines entire grid by factor 2
%   [v_refined, indices] = refine_velocity_grid(v, 3, [-1, 1]);  % Refines only range [-1, 1] by factor 3

    % Input validation
    if nargin < 2
        error('At least 2 arguments required: v and refine_factor');
    end
    
    if refine_factor < 2
        error('Refinement factor must be >= 2');
    end
    
    if nargin < 3
        % Refine entire grid
        v_range = [min(v), max(v)];
    end
    
    % Validate velocity range
    if length(v_range) ~= 2 || v_range(1) >= v_range(2)
        error('v_range must be a 2-element array [v_min, v_max] with v_min < v_max');
    end
    
    % Find indices corresponding to the velocity range
    v1_idx = find(v >= v_range(1), 1, 'first');
    v2_idx = find(v <= v_range(2), 1, 'last');
    
    if isempty(v1_idx) || isempty(v2_idx) || v1_idx >= v2_idx
        error('No valid velocity range found for refinement');
    end
    
    % Initialize refined grid
    v_refined = v;
    
    % Process each segment in the specified range
    for i = v1_idx:(v2_idx-1)
        % Get the current segment endpoints
        v_start = v(i);
        v_end = v(i+1);
        
        % Create refined segment
        % For refine_factor = n, we want (n-1) intermediate points
        num_intermediate = refine_factor - 1;
        
        if num_intermediate > 0
            % Create intermediate points
            v_intermediate = linspace(v_start, v_end, refine_factor + 1);
            % Remove the endpoints since they already exist
            v_intermediate = v_intermediate(2:end-1);
            
            % Insert intermediate points into the refined grid
            % Find where to insert (after the current point)
            insert_pos = find(v_refined == v_start, 1, 'last') + 1;
            
            % Insert the intermediate points
            v_refined = [v_refined(1:insert_pos-1), v_intermediate, v_refined(insert_pos:end)];
        end
    end
    
    % Ensure the refined grid is sorted
    v_refined = sort(v_refined);
    
    % Remove any duplicate points (shouldn't happen, but safety check)
    v_refined = unique(v_refined);
    
    % Now find the original indices in the refined grid
    original_indices = [];
    for i = 1:length(v)
        % Find the index of the original point in the refined grid
        idx = find(abs(v_refined - v(i)) < 1e-10, 1);
        if ~isempty(idx)
            original_indices = [original_indices, idx];
        end
    end
    
    % Compute velocity grid spacing for refined grid
    dv_refined = zeros(size(v_refined));
    
    % For interior points, use average of adjacent spacings
    for i = 2:length(v_refined)-1
        dv_refined(i) = (v_refined(i+1) - v_refined(i-1)) / 2;
    end
    
    % For boundary points, use adjacent spacing
    if length(v_refined) > 1
        dv_refined(1) = v_refined(2) - v_refined(1);
        dv_refined(end) = v_refined(end) - v_refined(end-1);
    else
        dv_refined(1) = 1;  % Default spacing for single point
    end
    
end
