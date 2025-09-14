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
        params.time_array(iT) = time;

        % Measurements
        [params]=measure(params, fs);

        % Plot results at each time step
        plot_results(params, fs);

        % Record CPU time for this iteration
        params.tcpu(iT) = toc(tic_iter);

        % Save config at specific times
        if mod(iT,params.dit_save) == 0
            Nsamples = Nsamples + 1;
            data = save_config(params,data,fs,Nsamples);
        end

        fprintf("iter: %d, time: %.1f, dt: %.2f \n",iT, time, params.dt)
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
        pcolor(params.grids(s).X, params.grids(s).V, fs(:, :, s));
        shading flat;
        subtitle("$f_\mathrm{"+params.species_name(s)+"}$");
        colorbar();
        xlabel("$x$");
        ylabel("$v$");
    end

    % Plot electric field
    subplot(params.Ns+1, 1, params.Ns+1);
    plot(params.grids(1).x, params.Efield);
    xlim([params.grids(1).x(1), params.grids(1).x(end)]);
    subtitle("$E$");
    colorbar();
    xlabel("$x$");

    pause(0.01); % Pause for visualization
end
