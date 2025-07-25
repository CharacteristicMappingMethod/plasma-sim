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

    iT = params.it + 1;
    dt = params.dt;
    N_remap = params.N_remap;

    % Persistent counter for N_nufi that increments across function calls
    persistent N_nufi;
    if isempty(N_nufi)
        N_nufi = 1;
    end

    % Loop over all particle species
    for s = 1:params.Ns
        % Increment N_nufi counter
        N_nufi = N_nufi + 1;

        % Apply symplectic flow for half step
        [X, V] = sympl_flow_Half(N_nufi, dt, params.grids(s).X, params.grids(s).V, ...
                                params.charge(s)/params.Mass(s)*params.Efield_list, ...
                                params.grids(s), "CMM");

        % Compose with existing maps
        [Xstar, Vstar] = wrap_compose(params, params.grids(s), ...
                                     squeeze(params.Map_stack(:,:,:,s,:)), X, V);

        % Evaluate initial condition at new positions
        fini = params.fini{s};
        fs(:,:,s) = fini(Xstar, Vstar);

        % Perform remapping at specified intervals
        if mod(iT, N_remap) == 0
            % Add current map to stack
            Nmaps = params.Nmaps + 1;
            params.Map_stack(:,:,1,s,Nmaps) = X;
            params.Map_stack(:,:,2,s,Nmaps) = V;
            params.Nmaps = Nmaps;

            % Recompute electric field
            [Efield] = vPoisson(fs, params.grids, params.charge);
            params.Efield_list(:,1) = Efield;

            % Reset N_nufi counter after remapping
            N_nufi = 1;
        end
    end

    % Update electric field for current step
    [Efield] = vPoisson(fs, params.grids, params.charge);
    params.Efield = Efield;
    params.Efield_list(:,N_nufi) = Efield;
end

function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid, method)
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
%
% Outputs:
%   X, V   - Updated position and velocity arrays

    mint = "lagrange";
    order = 3;

    % Return early if n=1 (no flow needed)
    if n == 1
        return;
    end

    % Periodic boundary function
    periodic = @(x) mod(x, grid.Lx - grid.dx);

    % Set up velocity field based on method
    if method == "CMM"
        Vperiodic = grid.Vperiodic;
        if mint == "lagrange"
            % Lagrange interpolation for velocity field
            Ux = @(X,V) reshape(lagrange2d_local_interp_periodic(X, V+grid.Lv, ...
                                                               grid.x, grid.v+grid.Lv, ...
                                                               Vperiodic, order), grid.size);
        else
            % Standard interpolation for velocity field
            Ux = @(X,V) interp2(grid.X, grid.V, Vperiodic, periodic(X), V, mint);
        end
    else
        % NuFi method: velocity field is just V
        Ux = @(X,V) V;
    end

    % Set up acceleration field
    if mint == "lagrange"
        % Lagrange interpolation for acceleration field
        Uv = @(X,V,E) -reshape(lagrange1d_local_interp_periodic(reshape(X,[],1), ...
                                                               grid.x, E(:), order), grid.size);
    else
        % Standard interpolation for acceleration field
        Uv = @(X,V,E) -reshape(interp1(grid.x, E, reshape(periodic(X),[],1), mint), grid.size);
    end

    % Apply symplectic flow: n-1 full steps
    while n > 2
        n = n - 1;
        X = X - dt * Ux(X,V);  % Position update (inverse signs for backward flow)
        V = V - dt * Uv(X,V,Efield(:,n));  % Velocity update
    end

    % Final half step
    X = X - dt * Ux(X,V);
    V = V - (dt / 2) * Uv(X,V,Efield(:,1));
end