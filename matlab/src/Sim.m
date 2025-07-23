function [params,data] = Sim(params)

    % Initialize grids, distribution functions and output array (data)
    [params, fs, data] = initialize_simulation(params);
    
    % Plot initial conditions
    params.Efield = vPoisson(fs,params.grids,params.charge);
    params.Efield_list(:,1) = params.Efield;
    
    % Initialize CPU timing arrays
    params.tcpu = zeros(params.Nt_max, 1);
    params.time_array = zeros(params.Nt_max, 1);
    
    % first plot the initial condition
    plot_results(params, fs);
    
    % Main loop over time
    Nsamples = 0;
    time = 0;
    for iT = 1:params.Nt_max
        tic_iter = tic();  % Start timing for this iteration
        
        params.it = iT;
        % Perform a single time step
        [fs, params] = step(params, fs);
        
        % increase time
        time = time + params.dt;
        params.time_array(iT) = time;
    
        % Measurements
        [params]=measure(params, fs);
    
        % Plot results at each time step
        plot_results(params, fs);
        
        % Record CPU time for this iteration
        params.tcpu(iT) = toc(tic_iter);
    
        % save config at specific times
        if mod(iT,params.dit_save) == 0
            Nsamples = Nsamples + 1;
            data = save_config(params,data,fs,Nsamples);
        end
        
        fprintf("iter: %d, time: %.1f, dt: %.2f\n",iT, time, params.dt)
        % Check if simulation end time is reached
        if time >= params.Tend
            % Trim arrays to actual size
            params.tcpu = params.tcpu(1:iT);
            params.time_array = params.time_array(1:iT);
            break;
        end
    end
    
    % save config to file
    save_config(params,data,fs,Nsamples,1);
    
    end
    %% Helper Functions
    
    function [fs, params] = step(params, fs)
    
    
        if params.method == "predcorr"
            if isfield(params,"dt_adapt_tolerance")
                [fs,params] = predictor_corrector_dt_adaptive(params,fs);
            else
                [fs,params] = predictor_corrector(params,fs);
            end
        elseif params.method == "predcorr_multi"
            [fs,params] = predictor_corrector_subcycling_electrons(params,fs);
        elseif params.method == "predcorr_hybrid"
            [fs,params] = predictor_corrector2(params,fs);
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
N_remap = params.N_remap;
for s = 1:params.Ns

    N_nufi = mod(iT-1,N_remap) + 1;   % This ensures N_nufi is always between 1 and 100
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
    if mod(iT,N_remap)==0
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
    

    
function [X,V] = sympl_flow_Half(n, dt, X,V, Efield, grid,method)
mint="lagrange";
order = 3;
if n == 1
    return;
end

periodic = @(x) mod(x,grid.Lx-grid.dx);

% this is the periodic continuation of v given in the arxiv paper VP-CMM.
if method == "CMM"
    Vperiodic = grid.Vperiodic;
    if mint=="lagrange"
        %%% attention grid values for the langrange interpolation should be
        %%% always between [0,L-dx] not negative!
    Ux = @(X,V) reshape(lagrange2d_local_interp_periodic(X, V+grid.Lv, grid.x, grid.v+grid.Lv, Vperiodic, order),grid.size);
    else
    Ux = @(X,V) interp2(grid.X, grid.V, Vperiodic, periodic(X), V, mint); 
    end
else % its the normal nufi method
    Ux =@(X,V) V;%interp2(grid.)
end

    if mint=="lagrange"
        % currently doesnt work??
        Uv = @(X,V,E) -reshape(lagrange1d_local_interp_periodic(reshape(X,[],1),grid.x,E(:),order),grid.size);
    else
        Uv = @(X,V,E) -reshape(interp1(grid.x,E,reshape(periodic(X),[],1),mint),grid.size);
    end


while n > 2
    n = n - 1;
    X = X - dt * Ux(X,V);  % Inverse signs; going backwards in time
    V = V - dt *Uv(X,V,Efield(:,n));
end

X = X - dt * Ux(X,V);
V = V - (dt / 2) *Uv(X,V,Efield(:,1));

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
    
    
