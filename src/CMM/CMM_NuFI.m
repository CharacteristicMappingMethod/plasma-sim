function [fs, params] = CMM_NuFI(params, fs)
% CMM-NuFI
%
% Inputs:
%   params - Simulation parameters structure
%   fs     - Distribution function array [Nv, Nx, Ns]
%
% Outputs:
%   fs     - Updated distribution function
%   params - Updated parameters structure

iT = params.it;
dt = params.dt;
N_remap = params.N_remap;

% Persistent counter for N_nufi that increments across function calls
persistent N_nufi;
if isempty(N_nufi)
    N_nufi = 1;
end
% Increment N_nufi counter
N_nufi = N_nufi + 1;

% Loop over all particle species
for s = 1:params.Ns

        % Use sample grid for flow and map grid for storage
        sample_grid = params.grids(s).sample;
        map_grid = params.grids(s).map;
        
        % Apply symplectic flow for half step on sample grid
        [Xrefined, Vrefined] = sympl_flow_Half(N_nufi, dt, sample_grid.X, sample_grid.V, ...
            params.charge(s)/params.Mass(s)*params.Efield_list, ...
            sample_grid, "CMM", params);

        % Compose with existing maps
        [Xstar, Vstar] = evaluate_map(squeeze(params.Map_stack(:,:,:,s,:)), map_grid, params, Xrefined, Vrefined);
        
        % Evaluate initial condition at new positions
        fini = params.fini{s};
        fs(:,:,s) = fini(Xstar, Vstar);
        
end


% Check if this is the final iteration
is_final_iteration = (iT == params.Nt_max) || (params.time + dt >= params.Tend);

% Perform remapping at specified intervals or on final iteration
if mod(iT, N_remap) == 0 || is_final_iteration
    % Add current maps to stack for all species
    Nmaps = params.Nmaps + 1;
    for s = 1:params.Ns
        params.Map_stack(:,:,1,s,Nmaps) = Xrefined(params.grids(s).idx_sample_to_map{1},params.grids(s).idx_sample_to_map{2});
        params.Map_stack(:,:,2,s,Nmaps) = Vrefined(params.grids(s).idx_sample_to_map{1},params.grids(s).idx_sample_to_map{2});
    end
    params.Nmaps = Nmaps;

    % Reset N_nufi counter after remapping
    N_nufi = 1;
end

% Update electric field for current step
[Efield] = vPoisson(fs, params.grids, params.charge);

% Add external field if defined
Efield = Efield + compute_external_Efield(params, params.grids(1).x, params.time + dt);


params.Efield = Efield;
params.Efield_list(:,N_nufi) = Efield;
end

function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid, method, params)
% Symplectic flow for half time step in CMM method
%
% This function applies a symplectic flow to advance particle positions
% and velocities for n-1 full steps plus one half step.
%
% Inputs:
%   n      - Number of steps (typically N_nufi)
%   dt     - Time step size
%   X, V   - Position and velocity arrays
%   Efield - Electric field array
%   grid   - Grid structure
%   method - Method type ("CMM" or "NuFi")
%   params - Parameters structure
%
% Outputs:
%   X, V   - Updated position and velocity arrays
opts = params.opt_interp;
% Return early if n=1 (no flow needed)
if n == 1
    return;
end

gridsize = size(V);
% Set up velocity field based on method
if method == "CMM"
    Vperiodic = grid.Vperiodic;
    Ux = @(X,V)interp2d_periodic(X, V+grid.Lv, grid.x, grid.v+grid.Lv, Vperiodic, opts); % all coordinates must be rescaled to [0,2pi]
else
    % NuFi method: velocity field is just V
    Ux = @(X,V) V;
end

% Set up acceleration field using direct interpolation
Uv = @(X,V,E) -reshape(interp1d_periodic(X(:), params.grids(1).x, E(:), params.opt_interp), gridsize);

% Apply symplectic flow based on number of maps
if params.Nmaps == 0
    % Vectorized implementation for Nmaps == 0
    if n > 2
        % Full steps: n-2 iterations
        for i = 1:(n-2)
            X = X - dt * Ux(X,V);  % Position update (inverse signs for backward flow)
            V = V - dt * Uv(X,V,Efield(:,n-i));  % Velocity update
        end
    end
    
    % Final half step
    X = X - dt * Ux(X,V);
    V = V - (dt / 2) * Uv(X,V,Efield(:,1));
else
    % Vectorized implementation for Nmaps > 0
    if n > 1
        % Full steps: n-1 iterations
        for i = 1:(n-1)
            X = X - dt * Ux(X,V);
            V = V - dt * Uv(X,V,Efield(:,n-i));
        end
    end
end
end