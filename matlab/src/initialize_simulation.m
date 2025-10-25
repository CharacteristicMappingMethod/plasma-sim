function [params, fs, data] = initialize_simulation(params)
% Initialize grids and distribution functions
for s = 1:params.Ns
    if params.Ns ==1
        Lv = params.Lv;
    else
        Lv = params.Lv_s(s);
    end
    
    % Check if Nsample and Nmap are defined, otherwise use Nx, Nv
    if isfield(params, 'Nsample') && isfield(params, 'Nmap')
        % Use separate sample and map grids
        Nsample = params.Nsample;
        Nmap = params.Nmap;
        
        % Expand scalar values to arrays
        if isscalar(Nsample)
            Nsample = [Nsample, Nsample];
        end
        if isscalar(Nmap)
            Nmap = [Nmap, Nmap];
        end
        
        % Validate that Nsample is a multiple of Nmap
        if any(mod(Nsample, Nmap) ~= 0)
            error('Nsample must be a multiple of Nmap for each dimension');
        end
        
        % Validate that Nmap is not bigger than Nsample
        if any(Nmap > Nsample)
            error('Nmap must not be bigger than Nsample');
        end
        
        % Automatically switch method based on grid configuration
        if ~isequal(Nsample, Nmap)
            % Use CMM_vargrid for dual-grid setup
            if params.method == "CMM"
                fprintf('Switching from CMM to CMM_vargrid due to dual-grid setup (Nsample=%s, Nmap=%s)\n', ...
                        mat2str(Nsample), mat2str(Nmap));
                params.method = "CMM_vargrid";
            end
        else
            % Use CMM for single-grid setup
            if params.method == "CMM_vargrid"
                fprintf('Switching from CMM_vargrid to CMM due to single-grid setup (Nsample=%s, Nmap=%s)\n', ...
                        mat2str(Nsample), mat2str(Nmap));
                params.method = "CMM";
            end
        end
        
        % Create sample grid
        [grid_sample] = make_periodic_grid(params.Lx, Lv, Nsample(1), Nsample(2));
        grid_sample.method = "spline";
        
        % Create map grid
        [grid_map] = make_periodic_grid(params.Lx, Lv, Nmap(1), Nmap(2));
        grid_map.method = "spline";
        
        % Store both grids
        params.grids(s).sample = grid_sample;
        params.grids(s).map = grid_map;
        
        % Create index mapping from sample grid to map grid
        % Since sample grid is finer, we need to map each sample point to nearest map point
        ratio_x = Nsample(1) / Nmap(1);
        ratio_v = Nsample(2) / Nmap(2);
        
        idx_sample_to_map = {1:ratio_x:Nsample(1), 1:ratio_v:Nsample(2)};
        params.grids(s).idx_sample_to_map = idx_sample_to_map;
        
        % Test that the index mapping is correct
        sample_X_mapped = grid_sample.X(idx_sample_to_map{1}, idx_sample_to_map{2});
        sample_V_mapped = grid_sample.V(idx_sample_to_map{1}, idx_sample_to_map{2});
        
        if ~isequal(sample_X_mapped, grid_map.X)
            error('Index mapping test failed: sample grid X coordinates do not match map grid X coordinates');
        end
        
        if ~isequal(sample_V_mapped, grid_map.V)
            error('Index mapping test failed: sample grid V coordinates do not match map grid V coordinates');
        end
    
        
        % Keep original grid structure for compatibility
        params.grids(s).x = grid_sample.x;
        params.grids(s).v = grid_sample.v;
        params.grids(s).X = grid_sample.X;
        params.grids(s).V = grid_sample.V;
        params.grids(s).Xsample_grid = grid_sample.X;
        params.grids(s).Vsample_grid = grid_sample.V;
        params.grids(s).size = grid_sample.size;
        params.grids(s).size_sample_grid = grid_sample.size_sample_grid;
        params.grids(s).dom = grid_sample.dom;
        params.grids(s).dx = grid_sample.dx;
        params.grids(s).dv = grid_sample.dv;
        params.grids(s).Lx = grid_sample.Lx;
        params.grids(s).Lv = grid_sample.Lv;
        params.grids(s).Nx = grid_sample.Nx;
        params.grids(s).Nv = grid_sample.Nv;
        params.grids(s).Dx = grid_sample.Dx;
        params.grids(s).Dv = grid_sample.Dv;
        params.grids(s).kx = grid_sample.kx;
        params.grids(s).kx2 = grid_sample.kx2;
        params.grids(s).kv = grid_sample.kv;
        params.grids(s).kv2 = grid_sample.kv2;
        params.grids(s).Vperiodic = grid_sample.Vperiodic;
        params.grids(s).Weights = grid_sample.Weights;
        
    else
        % Use traditional single grid
        [grid] = make_periodic_grid(params.Lx, Lv, params.Nx, params.Nv);
        grid.method = "spline";
        params.grids(s) = grid;
    end
end

if params.method=="CMM_vargrid" && isfield(params, 'refine_factor') && params.Ns==1 % works only for one species
    [v_refined_partial, original_indices, dv_refined] = refine_velocity_grid(params.grids(1).v, params.refine_factor, params.v_range);
    params.grids(1).v_refined = v_refined_partial;
    params.grids(1).idx_sample_to_map = {1:params.Nx, original_indices};
    params.grids(1).dv = dv_refined;
    [Xref,Vref] = meshgrid(params.grids(1).x,v_refined_partial);
    params.grids(1).Xsample_grid = Xref;
    params.grids(1).Vsample_grid = Vref;

    params.grids(1).size_sample_grid=size(Xref);
    
    % Test that the index mapping is correct for legacy case
    sample_X_mapped = Xref(params.grids(1).idx_sample_to_map{1}, params.grids(1).idx_sample_to_map{2});
    sample_V_mapped = Vref(params.grids(1).idx_sample_to_map{1}, params.grids(1).idx_sample_to_map{2});
    
    if ~isequal(sample_X_mapped, params.grids(1).X)
        error('Legacy index mapping test failed: sample grid X coordinates do not match original grid X coordinates');
    end
    
    if ~isequal(sample_V_mapped, params.grids(1).V)
        error('Legacy index mapping test failed: sample grid V coordinates do not match original grid V coordinates');
    end

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
        if params.method == "CMM_vargrid"
            params.Map_stack = zeros(params.grids(s).map.size(1), params.grids(s).map.size(2), 2, params.Ns, params.Nmaps_max);
        else
            params.Map_stack = zeros(params.grids(s).size(1), params.grids(s).size(2), 2, params.Ns, params.Nmaps_max);
        end
        
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
    data.time = dt_save*(1:Nsamples);

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
    if params.incomp_error_threshold<=0
        params.incomp_error_threshold=1e+10;
    else
        params.N_remap = params.Nt_max; % if incomp_error_threshold is set, remapping parameter is set to Nt_max to avoid "manual" remapping
    end
end


if ~exist(params.data_dir, 'dir')
    mkdir(params.data_dir);
end

end