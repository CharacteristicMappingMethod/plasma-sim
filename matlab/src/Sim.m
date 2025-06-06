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
    
        N_nufi = mod(iT-1,10) + 1;   % This ensures N_nufi is always between 1 and 10
        if N_nufi ==1
            N_nufi = N_nufi+1;
        end  
        [X,V] = sympl_flow_Half(N_nufi,dt,params.grids(s).X,params.grids(s).V,params.charge(s)/params.Mass(s)*params.Efield_list,params.grids(s), "CMM");
        fini = params.fini{s};
        fs(:,:,s) = fini(X,V);
    
    % 3.) Remapping ?
    % e.g. after mod(iT,10) do remapping
        if mod(iT,10)==0
            % add mapstack
            Nmaps = params.Nmaps + 1;
            params.Map_stack(:,:,1,s,Nmaps) = X;
            params.Map_stack(:,:,2,s,Nmaps) = V;
            params.Nmaps = Nmaps;

            composed_map = compose_maps_numerical(squeeze(params.Map_stack(:,:,:,s,:)), params.grids(s),params); 
   
            Xstar = composed_map(:,:,1);
            Vstar = composed_map(:,:,2);
            %reinitailize the E filed and inital condition
            
            
            fini = params.fini{s};
            fs(:,:,s) = fini(Xstar,Vstar);
            [Efield] =vPoisson(fs,params.grids,params.charge);
        end
    end
    
 
[Efield] =vPoisson(fs,params.grids,params.charge);
params.Efield = Efield;
params.Efield_list(:,N_nufi) = Efield;
    
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
        Ux = @(X,V) interp2(grid.X, grid.V, Vperiodic, X, V, mint); % TODO: check if this is correct, because V is not periodic 
    else % its the normal nufi method
        Ux =@(X,V) V;%interp2(grid.)
    end
    
    Uv = @(X,V,E) reshape(interp1(grid.x,E,reshape(periodic(X),[],1),mint),grid.size);
    
    % For CMM method, we need to go backwards in time using the stored fields
    field_index = 1;
    while n > 2
        n = n - 1;
        field_index = field_index + 1;
        X = X - dt * Ux(X,V);  % Inverse signs; going backwards in time
        if field_index <= size(Efield, 2)
            V = V + dt * Uv(X,V,Efield(:,field_index));
        else
            % If we don't have enough history, use the last available field
            V = V + dt * Uv(X,V,Efield(:,end));
        end
    end
    
    X = X - dt * Ux(X,V);
    V = V + (dt / 2) * Uv(X,V,Efield(:,1));
    
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
    
    
    