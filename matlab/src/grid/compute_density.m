function [rho] = compute_density(fs, dv)
% COMPUTE_DENSITY Computes charge density from distribution function on uniform grid
%
% This function computes the charge density by integrating the distribution
% function over velocity space. For a uniform velocity grid, this is equivalent
% to a rectangular rule integration.
%
% INPUTS:
%   fs  - Distribution function (Nv x Nx array)
%        fs(i,j) = f(v_i, x_j) where v_i is velocity and x_j is position
%   dv  - Velocity grid spacing (scalar)
%        For uniform grids: dv = v(2) - v(1)
%
% OUTPUTS:
%   rho - Charge density (1 x Nx array)
%        rho(j) = sum_i(fs(i,j) * dv) for j = 1,...,Nx
%
% MATHEMATICAL BACKGROUND:
%   The charge density is defined as:
%   ρ(x) = ∫ f(x,v) dv
%
%   For a discrete uniform grid with spacing dv:
%   ρ(x_j) ≈ dv * Σ_i f(x_j, v_i)
%
%   This is the rectangular rule approximation to the integral.
%
% VECTORIZATION:
%   The function uses vectorized operations for efficiency:
%   - fs .* dv(:) multiplies each row of fs by the velocity spacing
%   - sum(..., 1) sums over velocity (rows) to get density at each x
%
% EXAMPLE:
%   % Create test data
%   Nv = 100; Nx = 50;
%   v = linspace(-3, 3, Nv);
%   x = linspace(0, 2*pi, Nx);
%   dv = v(2) - v(1);
%   fs = exp(-v'.^2/2) * ones(1,Nx);  % Maxwellian distribution
%   
%   % Compute density
%   rho = compute_density(fs, dv);
%
% See also: refine_velocity_grid

    rho = sum(fs .* dv(:), 1);

end
