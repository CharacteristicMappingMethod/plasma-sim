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

% Initialize distribution functions
for s = 1:params.Ns
    fini = params.fini{s};
    fs(:, :, s) = fini(params.grids(s).X,params.grids(s).V);
    % Method specific initialization:
    if params.method == "CMM"
       
        params.Nmaps=0;
   
        params.Map_stack(:,:,:,s,:) = zeros(params.grids(s).size(1), params.grids(s).size(2), params.Ns, 1); 
    
    end
end

% Maximal Iteration number:
% check if maximal time iteration number Nt_max fits final time Tend
% this is necessary since size of t = 0he preallocated arrays may depend on
% Nt_max.
if params.Nt_max> params.Tend/params.dt
    params.Nt_max = ceil(params.Tend/params.dt);
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
    Nsize = [params.grids(1).size(:)',Nsamples,params.Ns];
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

if ~isfield(params, 'incomp_error_threshold')
    params.incomp_error_threshold = 1e+10;
else
    params.N_remap = params.Nt_max; % if incomp_error_threshold is set, remapping parameter is set to Nt_max to avoid "manual" remapping
end


if ~exist(params.data_dir, 'dir')
    mkdir(params.data_dir);
end

end