function [fs, params] = NuFi(params,fs)
iT = params.it+1;
dt = params.dt;
for s = 1:params.Ns
    [X,V] = sympl_flow_Half(iT,dt,params.grids(s).X,params.grids(s).V,params.charge(s)/params.Mass(s)*params.Efield_list,params.grids(s),params);
    fini = params.fini{s};
    fs(:,:,s) = fini(X,V);
end
% Compute electric field
[Efield] =vPoisson(fs,params.grids,params.charge);
% Add external field
Efield = Efield + compute_external_Efield(params, params.grids(1).x, params.time + dt);
params.Efield = Efield;
params.Efield_list(:,iT) = Efield;

end






function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid, params)
% Symplectic flow for half time step in NuFi method
%
% This function applies a symplectic flow to advance particle positions
% and velocities for n-1 full steps plus one half step.
%
% Inputs:
%   n      - Number of steps
%   dt     - Time step size
%   X, V   - Position and velocity arrays
%   Efield - Electric field array
%   grid   - Grid structure
%   params - Parameters structure
%
% Outputs:
%   X, V   - Updated position and velocity arrays

    % Return early if n=1 (no flow needed)
    if n == 1
        return;
    end

    % Set up velocity field (NuFi method: velocity field is just V)
    Ux = @(X,V) V;
    
    % Set up acceleration field using direct interpolation
    Uv = @(X,V,E) -reshape(interp1d_periodic(X(:), params.grids(1).x, E(:), params.opt_interp), grid.size);

    % Apply symplectic flow
    while n > 2
        n = n - 1;
        X = X - dt * Ux(X,V);  % Position update (inverse signs for backward flow)
        V = V - dt * Uv(X,V,Efield(:,n));  % Velocity update
    end

    % Final half step
    X = X - dt * Ux(X,V);
    V = V - (dt / 2) * Uv(X,V,Efield(:,1));
end

