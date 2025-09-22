function [fs, params] = CMM(params, fs)
% CMM (Conservative Map Method) implementation for plasma simulation
%
% This function implements the Conservative Map Method for advancing
% the distribution function in time using symplectic flow and map composition.
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
% Initialize coordinate maps storage if not exists
if ~isfield(params, 'coordinate_maps')
    params.coordinate_maps.X = zeros(params.grids(1).size(1), params.grids(1).size(2), params.Ns);
    params.coordinate_maps.V = zeros(params.grids(1).size(1), params.grids(1).size(2), params.Ns);
    % Initialize storage for current maps
    params.coordinate_maps.submap_X = zeros(params.grids(1).size(1), params.grids(1).size(2), params.Ns);
    params.coordinate_maps.submap_V = zeros(params.grids(1).size(1), params.grids(1).size(2), params.Ns);
end
% Loop over all particle species
for s = 1:params.Ns

    % Apply symplectic flow for half step
    [X, V] = sympl_flow_Half(N_nufi, dt, params.grids(s).X, params.grids(s).V, ...
        params.charge(s)/params.Mass(s)*params.Efield_list, ...
        params.grids(s), "CMM", params);

    % Store current maps for potential remapping
    params.coordinate_maps.submap_X(:,:,s) = X;
    params.coordinate_maps.submap_V(:,:,s) = V;

    [detJ, ~, ~, ~, ~] = jacobian_determinant(X, V, params.grids(s));

    % we multiply with the weighting function to not count for
    % incompressibility errors on the periodified edges which are very
    % distored
    params.incomp_error(s)=max(abs(log(detJ(:))).*params.grids(s).Weights(:));
    % Compose with existing maps
    [Xstar, Vstar] = wrap_compose(params, params.grids(s), ...
        squeeze(params.Map_stack(:,:,:,s,:)), X, V);

    % Store coordinate maps for analysis
    params.coordinate_maps.X(:,:,s) = Xstar;
    params.coordinate_maps.V(:,:,s) = Vstar;

    % Evaluate initial condition at new positions
    fini = params.fini{s};
    fs(:,:,s) = fini(Xstar, Vstar);
end

params.max_incomp_error = max(params.incomp_error);

% Check if this is the final iteration
is_final_iteration = (iT == params.Nt_max) || (params.time + dt >= params.Tend);

% Perform remapping at specified intervals or on final iteration
if mod(iT, N_remap) == 0 || params.max_incomp_error > params.incomp_error_threshold || is_final_iteration
    % Add current maps to stack for all species
    Nmaps = params.Nmaps + 1;
    for s = 1:params.Ns
        params.Map_stack(:,:,1,s,Nmaps) = params.coordinate_maps.submap_X(:,:,s);
        params.Map_stack(:,:,2,s,Nmaps) = params.coordinate_maps.submap_V(:,:,s);
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

% Set up velocity field based on method
if method == "CMM"
    Vperiodic = grid.Vperiodic;
    Ux = @(X,V)interp2d_periodic(X, V+grid.Lv, grid.x, grid.v+grid.Lv, Vperiodic, opts); % all coordinates must be rescaled to [0,2pi]
else
    % NuFi method: velocity field is just V
    Ux = @(X,V) V;
end

% Set up acceleration field using direct interpolation
Uv = @(X,V,E) -reshape(interp1d_periodic(X(:), params.grids(1).x, E(:), params.opt_interp), grid.size);
% Apply symplectic flow based on number of maps

if params.Nmaps == 0

    while n > 2
        n = n - 1;
        X = X - dt * Ux(X,V);  % Position update (inverse signs for backward flow)
        V = V - dt * Uv(X,V,Efield(:,n));  % Velocity update
    end

    % Final half step
    X = X - dt * Ux(X,V);
    V = V - (dt / 2) * Uv(X,V,Efield(:,1));
else

    while n > 1
        n = n - 1;
        X = X - dt * Ux(X,V);
        V = V - dt * Uv(X,V,Efield(:,n));
    end
end
end