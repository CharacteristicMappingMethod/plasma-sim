function [params, fs, data] = initialize_simulation(params)
% Initialize grids and distribution functions
for s = 1:params.Ns
    if params.Ns ==1
        Lv = params.Lv;
    else
        Lv = params.Lv_s(s);
    end
    [grid] = make_periodic_grid(params.Lx,Lv,params.Nx,params.Nv);
    grid.method = "spline";
    params.grids(s) = grid;
end

if params.method=="CMM_vargrid" % works only for one species
    if params.Ns>1
        error("CMM_vargrid currently does not support Ns>1");
    end
    [v_refined_partial, original_indices, dv_refined] = refine_velocity_grid(params.grids(1).v, params.refine_factor, params.v_range);
    params.grids(1).v_refined = v_refined_partial;
    params.grids(1).original_indices = original_indices;
    params.grids(1).dv = dv_refined;
    [Xref,Vref] = meshgrid(params.grids(1).x,v_refined_partial);
    params.grids(1).Xsample_grid = Xref;
    params.grids(1).Vsample_grid = Vref;

    params.grids(1).size_sample_grid=size(Xref);
end


% Maximal Iteration number:
% check if maximal time iteration number Nt_max fits final time Tend
% this is necessary since size of t = 0he preallocated arrays may depend on
% Nt_max.
if params.Nt_max> params.Tend/params.dt
    params.Nt_max = ceil(params.Tend/params.dt);
end


% Initialize distribution functions
for s = 1:params.Ns
    fini = params.fini{s};
    fs(:, :, s) = fini(params.grids(s).Xsample_grid,params.grids(s).Vsample_grid);
    % Method specific initialization:
    if params.method == "CMM" || params.method == "CMM_vargrid"
        params.max_incomp_error = 0;
        params.Nmaps=0;
        
        % Calculate maximum number of maps needed based on remapping frequency
        % Maps are stored every N_remap iterations, plus one for final iteration
        % Also account for potential remapping due to incompressibility threshold
        params.Nmaps_max = ceil(params.Nt_max / params.N_remap) + 1;
        
        % Add safety margin for potential threshold-based remapping (20% extra)
        params.Nmaps_max = ceil(params.Nmaps_max);
        
        % Validate Nmaps_max is reasonable (prevent excessive memory allocation)
        if params.Nmaps_max > 10000
            warning('Nmaps_max is very large (%d). Consider increasing N_remap or reducing Nt_max.', params.Nmaps_max);
        end
        
        % Ensure minimum size for edge cases
        if params.Nmaps_max < 2
            params.Nmaps_max = 2;
        end
        
        % Initialize map stack with calculated maximum size
        params.Map_stack = zeros(params.grids(s).size(1), params.grids(s).size(2), 2, params.Ns, params.Nmaps_max);
        
        % Debug information
        fprintf('CMM initialization: Nt_max=%d, N_remap=%d, Nmaps_max=%d\n', ...
                params.Nt_max, params.N_remap, params.Nmaps_max);
    
    end
end

% output data storage allocation
if isfield(params, 'dt_save')
    dt_save = params.dt_save; % how often do we want to save solution?
    dit_save = dt_save/params.dt;
    params.dit_save = dit_save;
    if rem(dit_save,1)~=0 || dt_save < params.dt
        % dit_save is not an integer
        % therefore dt is not an divisor of dt_save
        error('Huston we have a problem: dt_save is not correct.');
    end
    % allocate
    Nsamples = params.Nt_max/dit_save;
    Nsize = [params.grids(1).size_sample_grid(:)',Nsamples,params.Ns];
    data.fs=zeros(Nsize);
    data.Efield = zeros([params.grids(1).Nx,Nsamples]);
    data.time = dt_save*[1:Nsamples];

else
    % never save anything
    params.dit_save = params.Nt_max + 2;
    data = [];
end



% Directories
if ~isfield(params, 'data_dir')
    if ~isfield(params, 'root_dir')
        root = "./";
    else
        root = params.root_dir;
    end
    params.data_dir= root + "/data/" + params.mycase+"_Tend"+num2str(params.Tend)+"_"+params.method+"/";
end

% default interpolation parameters
if ~isfield(params, 'opt_interp')
    params.opt_interp = struct('scheme', 'lagrange-bary', 'order', 3, 'use_mex', true);
end

% default plotting frequency
if ~isfield(params, 'plot_freq')
    params.plot_freq = 1;
else
    if params.plot_freq == 0
        params.plot_freq = params.Nt_max;
    end
end

% default plotting frequency
if ~isfield(params, 'measure_freq')
    params.measure_freq = 1;
else
    if params.measure_freq == 0
        params.measure_freq = params.Nt_max;
    end
end

if ~isfield(params, 'incomp_error_threshold')
    params.incomp_error_threshold = 1e+10;
else
    params.N_remap = params.Nt_max; % if incomp_error_threshold is set, remapping parameter is set to Nt_max to avoid "manual" remapping
end


if ~exist(params.data_dir, 'dir')
    mkdir(params.data_dir);
end

end