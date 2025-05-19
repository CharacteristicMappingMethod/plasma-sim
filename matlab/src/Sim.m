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
    N_nufi = mod(iT,1000);
    [X,V] = sympl_flow_Half(N_nufi,dt,params.grids(s).X,params.grids(s).V,params.charge(s)/params.Mass(s)*params.Efield_list,params.grids(s), "CMM");
    % 1.) do the map composition
    Map_values = params.Map{s}(X,V);
    Xstar = Map_values(:,:,1);
    Vstar = Map_values(:,:,2);
    % 2.) compose with inicond
    fini = params.fini{s};
    fs(:,:,s) = fini(Xstar,Vstar);
end

% 3.) Remapping ?
% e.g. after mod(iT,100) do remapping
if mod(iT,1000)==0
        % add mapstack
        Nmaps = params.Nmaps + 1;
        params.Map_stack(:,:,1,s,Nmaps) = X;
        params.Map_stack(:,:,2,s,Nmaps) = V;
        params.Nmaps = Nmaps;

        params.Map = compose_maps(params.Map_stack); % TODO
end

% map1 = @(X,V) cat(3,X-10*V,V);
% V1 = grid.V;
% X1 = grid.X - 10*V1;

% map2 = @(X,V) cat(3,X-20*V,V);
% V2 = grid.V;
% X2 = grid.X - 20*V2;

% map3 = map1(map2(grid.X,grid.V))
% V3 = V1 \circ V2
% X3 = grid.X - 30 * grid.V;



[Efield] =vPoisson(fs,params.grids,params.charge);
params.Efield = Efield;
params.Efield_list(:,N_nufi) = Efield;


end

function [backward_map] = compose_maps(map_stack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets a stack of maps which are needed to 
% be composed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for map = map_stack
    map
end

backward_map = current_map;

end

% function [fs,params] = predictor_corrector_cmm(params,fs)
% 
%     iT = params.it+1;
% 
%     % Compute electric field
%     [Efield] = vPoisson(fs, params.grids, params.charge);
%     
%     % Advect distribution functions for half time step
%     for s = 1:params.Ns
%         f12(:, :, s) = Advect(fs(:, :, s), params.charge(s) / params.Mass(s) * Efield, params.grids(s), params.dt / 2);
%     end
% 
%     % Recompute electric field
%     [Efield] = vPoisson(f12, params.grids, params.charge);
% 
%     % Advect distribution functions for full time step
%     for s = 1:params.Ns
%         fs(:, :, s) = Advect(fs(:, :, s), params.charge(s) / params.Mass(s) * Efield, params.grids(s), params.dt);
%     end
% 
%     % save electric field;
%     params.Efield = Efield;
%     params.Efield_list(:,iT) = Efield;
% 
%     
% end




function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid, method )
mint="spline";
if n == 1
    return;
end

periodic = @(x) mod(x,grid.Lx-grid.dx);

% this is the periodic continuation of v given in the arxiv paper VP-CMM.
if method == "CMM"
    Vperiodic = grid.Vperiodic;
    Ux = @(X,V) interp2(grid.X, grid.V, Vperiodic, X, V, mint);
else % its the normal nufi method
    Ux =@(X,V) V;%interp2(grid.)
end

Uv = @(X,V,E) reshape(interp1(grid.x,E,reshape(periodic(X),[],1),mint),grid.size);


while n > 2
    n = n - 1;
    X = X - dt * Ux(X,V);  % Inverse signs; going backwards in time
    V = V + dt * Uv(X,V,Efield(:,n));
end

X = X - dt * Ux(X,V);
V = V + (dt / 2) *Uv(X,V,Efield(:,1));

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


