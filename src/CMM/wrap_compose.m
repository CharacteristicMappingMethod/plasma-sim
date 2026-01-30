function [Xstar, Vstar] = wrap_compose(params, grid, Map_stack, X, V)
% Compose maps by adding current map to stack and computing composition
%
% This function adds the current position and velocity maps to the existing
% map stack and computes the composition of all maps in the stack.
%
% Inputs:
%   params    - Simulation parameters structure
%   grid      - Grid structure
%   Map_stack - Array of existing maps [Nx, Nv, 2, Nmaps]
%   X, V      - Current position and velocity maps
%
% Outputs:
%   Xstar, Vstar - Composed position and velocity maps

    Nmaps = params.Nmaps;

    % Add current maps to the stack
    Map_stack(:,:,1,Nmaps+1) = X;
    Map_stack(:,:,2,Nmaps+1) = V;

    % Temporarily increment map count for composition
    params.Nmaps = Nmaps + 1;

    % Compute composition of all maps in stack
    composed_map = compose_maps_numerical(Map_stack, grid, params);

    % Restore original map count
    params.Nmaps = Nmaps;

    % Extract position and velocity components
    Xstar = composed_map(:,:,1);
    Vstar = composed_map(:,:,2);
end