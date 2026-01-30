function [params,data] = Sim(params)
% Main simulation function for plasma physics simulation
%
% This function runs the main simulation loop, handling initialization,
% time stepping, measurements, plotting, and data saving.
%
% Inputs:
%   params - Simulation parameters structure
%
% Outputs:
%   params - Updated parameters structure
%   data   - Simulation data structure

    % Initialize grids, distribution functions and output array (data)
    [params, fs, data] = initialize_simulation(params);

    % Plot initial conditions
    params.Efield = vPoisson(fs,params.grids,params.charge);
    params.Efield_list(:,1) = params.Efield;
    params.time = 0;

    % Initialize CPU timing arrays
    params.tcpu = zeros(params.Nt_max, 1);
    params.time_array = zeros(params.Nt_max, 1);

    % First plot the initial condition
    plot_results(params, fs);

    % Main loop over time
    Nsamples = 0;
    time = 0;
    for iT = 1:params.Nt_max
        tic_iter = tic();  % Start timing for this iteration

        params.it = iT;
        % Perform a single time step
        [fs, params] = step(params, fs);

        % Increase time
        time = time + params.dt;
        params.time = time;
        params.time_array(iT) = time;

        % Measurements
        if mod(iT,params.measure_freq) ==0
            [params]=measure(params, fs);
        end

        % Plot results at each time step
        if mod(iT,params.plot_freq) == 0
            plot_results(params, fs);
        end
        % Record CPU time for this iteration
        params.tcpu(iT) = toc(tic_iter);

        % Save config at specific times
        if mod(iT,params.dit_save) == 0
            Nsamples = Nsamples + 1;
            data = save_config(params,data,fs,Nsamples);
        end
        
        if params.method == "CMM" || params.method == "CMM_vargrid"
            fprintf("iter: %d, time: %.1f, dt: %.2f, incomp_error: %.2e Nmaps: %d, cpu_time: %.3f s\n",iT, time, params.dt, params.max_incomp_error, params.Nmaps, params.tcpu(iT))
        else
            fprintf("iter: %d, time: %.1f, dt: %.2f, cpu_time: %.3f s\n",iT, time, params.dt, params.tcpu(iT))
        end
        
        % Check if simulation end time is reached
        if time >= params.Tend
            % Trim arrays to actual size
            params.tcpu = params.tcpu(1:iT);
            params.time_array = params.time_array(1:iT);
            break;
        end
    end

    % Save config to file
    save_config(params,data,fs,Nsamples,1);
    
    % Print final timing summary
    total_cpu_time = sum(params.tcpu);
    avg_cpu_time = mean(params.tcpu);
    fprintf("\n=== Simulation Complete ===\n");
    fprintf("Total iterations: %d\n", iT);
    fprintf("Total CPU time: %.3f s (%.2f min)\n", total_cpu_time, total_cpu_time/60);
    fprintf("Average time per iteration: %.3f s\n", avg_cpu_time);
    fprintf("Last iteration CPU time: %.3f s\n", params.tcpu(end));
end

%% Helper Functions

function [fs, params] = step(params, fs)
% Perform a single time step using the specified method
%
% This function dispatches to the appropriate time stepping method
% based on the params.method field.
%
% Inputs:
%   params - Simulation parameters structure
%   fs     - Distribution function array
%
% Outputs:
%   fs     - Updated distribution function
%   params - Updated parameters structure

    if params.method == "predcorr"
        [fs,params] = predictor_corrector(params,fs);    
    elseif params.method == "NuFi"
        [fs,params] = NuFi(params,fs);
    elseif params.method == "CMM-NuFI"
        [fs,params] = CMM_NuFI(params,fs);
    else
        error("Unknown method: %s", params.method);
    end
end

function plot_results(params, fs)
% Plot results at each time step
%
% This function creates subplots showing the distribution function
% for each species and the electric field.
%
% Inputs:
%   params - Simulation parameters structure
%   fs     - Distribution function array

    % Plot distribution function for each species
    for s = 1:params.Ns
        subplot(params.Ns+1, 1, s);
        pcolor(params.grids(s).Xsample_grid, params.grids(s).Vsample_grid, fs(:, :, s));
        shading flat;
        subtitle("$f_\mathrm{"+params.species_name(s)+"}$");
        colorbar();
        xlabel("$x$");
        ylabel("$v$");
    end

    % Get electric field with external field added
    Efield = params.Efield;

    % Plot electric field
    subplot(params.Ns+1, 1, params.Ns+1);
    plot(params.grids(1).x, Efield);
    xlim([params.grids(1).x(1), params.grids(1).x(end)]);
    subtitle("$E$");
    colorbar();
    xlabel("$x$");

    pause(0.01); % Pause for visualization
end
