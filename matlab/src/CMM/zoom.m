function [f_zoom] = zoom(params, X, V, s)
% ZOOM - Apply composition mapping to evaluate distribution function
%
% Inputs: params (simulation parameters), X,V (coordinates), s (species, default=1)
% Output: f_zoom (distribution function values, 2D array)

% Handle optional species parameter
if nargin < 4
    % If no species specified, default to species 1
    s = 1;
end

% Initialize output array for distribution function values
fs = zeros(size(X,1), size(X,2));

% Check that input coordinates are within simulation domain for this species
if isfield(params.grids(s), 'dom') && length(params.grids(s).dom) >= 4
    dom = params.grids(s).dom;  % [x_min, v_min, x_max, v_max]
    x_min = dom(1); x_max = dom(3);
    v_min = dom(2); v_max = dom(4);
    
    if any(X(:) < x_min) || any(X(:) > x_max)
        warning('Some X coordinates are outside the simulation domain [%.3f, %.3f] for species %d', ...
                x_min, x_max, s);
    end
    
    if any(V(:) < v_min) || any(V(:) > v_max)
        warning('Some V coordinates are outside the simulation domain [%.3f, %.3f] for species %d', ...
                v_min, v_max, s);
    end
end

% Apply composition mapping to transform coordinates (X,V) -> (Xstar,Vstar)
% This uses the mapping stored in Map_stack for species s
%[Xstar, Vstar] = wrap_compose(params, params.grids(s),    squeeze(params.Map_stack(:,:,:,s,:)), X, V);
[X_star, V_star] = evaluate_map(squeeze(params.Map_stack(:,:,:,s,:)), params.grids(s), params, X, V);
% Evaluate initial condition at new positions
fini = params.fini{s};
fs = fini(X_star, V_star);

% Return the distribution function values
f_zoom = fs;
end