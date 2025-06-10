function [params,data] = Sim(params)

% Initialize grids, distribution functions and output array (data)
[params, fs, data] = initialize_simulation(params);

% Plot initial conditions
params.Efield = vPoisson(fs,params.grids,params.charge);
params.Efield_list(:,1) = params.Efield;

% first plot the initial condition
plot_results(params, fs);

% Main loop over time
Nsamples = 0;
for iT = 1:params.Nt_max
    params.it = iT;
    % Perform a single time step
    [fs, params] = step(params, fs);
    
    % increase time
    time = params.dt * iT

    % Measurements
    [params]=measure(params, fs);

    % Plot results at each time step
    plot_results(params, fs);

    % save config at specific times
    if mod(iT,params.dit_save) == 0
        Nsamples = Nsamples + 1;
        data = save_config(params,data,fs,Nsamples);
    end

    % Check if simulation end time is reached
    if time >= params.Tend
        break;
    end
end

% save config to file
save_config(params,data,fs,Nsamples,1);

end
%% Helper Functions

function [fs, params] = step(params, fs)


    if params.method == "predcorr"
        [fs,params] = predictor_corrector(params,fs);
    elseif params.method == "predcorr_multi"
        [fs,params] = predictor_corrector_subcycling_electrons(params,fs);
    elseif params.method == "NuFi"
        [fs,params] = NuFi(params,fs);
    elseif params.method == "CMM"
        [fs,params] = CMM(params,fs);
    else
        display("error step")
    end


end

function [fs,params] = CMM(params,fs)
iT = params.it+1;
dt = params.dt;

for s = 1:params.Ns

    N_nufi = mod(iT-1,50) + 1;   % This ensures N_nufi is always between 1 and 100
    if N_nufi ==1
        N_nufi = N_nufi+1;
    end
  
    [X,V] = sympl_flow_Half(N_nufi,dt,params.grids(s).X,params.grids(s).V,params.charge(s)/params.Mass(s)*params.Efield_list,params.grids(s), "CMM");
    % 1.) do the map composition
%     Map_values = params.Map{s}(X,V);
    [Xstar,Vstar] = wrap_compose(params,params.grids(s),squeeze(params.Map_stack(:,:,:,s,:)),X,V); 
%     Xstar = Map_values(:,:,1);
%     Vstar = Map_values(:,:,2);
    % 2.) compose with inicond
    fini = params.fini{s};
    fs(:,:,s) = fini(Xstar,Vstar);

% 3.) Remapping ?
% e.g. after mod(iT,100) do remapping
    if mod(iT,50)==0
        % add mapstack
        Nmaps = params.Nmaps + 1;
        params.Map_stack(:,:,1,s,Nmaps) = X;
        params.Map_stack(:,:,2,s,Nmaps) = V;
        params.Nmaps = Nmaps;
%         composed_map = compose_maps_numerical(squeeze(params.Map_stack(:,:,:,s,:)), params.grids(s),params); % TODO
%         params.Map{s} = @(X,V) composed_map;
%         %reinitailize the E filed
        [Efield] =vPoisson(fs,params.grids,params.charge);
        params.Efield_list(:,1) = Efield;
    end
end


[Efield] =vPoisson(fs,params.grids,params.charge);
params.Efield = Efield;
params.Efield_list(:,N_nufi) = Efield;

end



function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid, method )
mint="spline";
order = 4; % lagrange interpolation order
if n == 1
    return;
end


periodic = @(x) mod(x,grid.Lx-grid.dx);

while n > 2
    n = n - 1;
    X = X - dt * V;  % Inverse signs; going backwards in time
    if mint=="lagrange"
        V = V + dt *reshape(lagrange_local_interp_periodic(reshape(X,[],1),grid.x,Efield(:,n),order),grid.size);
    else
        V = V + dt *reshape(interp1(grid.x,Efield(:,n),reshape(periodic(X),[],1),mint),grid.size);
    end

end

X = X - dt * V;
if mint=="lagrange"
    V = V + (dt / 2) *reshape(lagrange_local_interp_periodic(reshape(X,[],1),grid.x,Efield(:,1),order),grid.size);
else
    V = V + (dt / 2) *reshape(interp1(grid.x,Efield(:,1),reshape(periodic(X),[],1),mint),grid.size);
end

end



function plot_results(params, fs)
    % Plot results at each time step
    for s = 1:params.Ns
    subplot(params.Ns+1, 1, s);
    pcolor(params.grids(s).X, params.grids(s).V, fs(:, :, s)); shading flat;
    subtitle("$f_\mathrm{"+params.species_name(s)+"}$");
    colorbar();
    xlabel("$x$");
    ylabel("$v$");
    end
   
    subplot(params.Ns+1, 1, params.Ns+1);
    plot(params.grids(1).x, params.Efield); xlim([params.grids(1).x(1), params.grids(1).x(end)]);
    subtitle("$E$");
    colorbar();
    xlabel("$x$");

    pause(0.01); % Pause for visualization
end

function [Xstar, Vstar ] = wrap_compose(params, grid, Map_stack, X,V)
    Nmaps = params.Nmaps;    
    Map_stack(:,:,1,Nmaps+1) = X;
    Map_stack(:,:,2,Nmaps+1) = V;    
    params.Nmaps = Nmaps+1;
    composed_map = compose_maps_numerical(Map_stack,grid,params);  
    params.Nmaps = Nmaps;
    Xstar = composed_map(:,:,1);
    Vstar = composed_map(:,:,2);
end