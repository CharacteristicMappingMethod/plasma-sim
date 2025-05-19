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
        params.Map{s} = @(X,V) cat(3,X,V);
    end
end

% Maximal Iteration number:
% check if maximal time iteration number Nt_max fits final time Tend
% this is necessary since size of the preallocated arrays may depend on
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

if ~exist(params.data_dir, 'dir')
    mkdir(params.data_dir);
end

end