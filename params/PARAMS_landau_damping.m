% Initialize parameter struct
params = struct();

%% General Simulation Parameters
params.mycase       = "landau_damping";   % Options: "landau_damping", "two_stream"
params.method       = "predcorr";         % Time integration method
params.Nt_max       = 4000000;               % Max number of time steps
params.dt           = 0.1;                % Time step size
params.Tend         = 100;                 % Simulation end time
params.dt_save      = 10;                 % Save data every dt_save time units
params.N_remap = 20;
%% Grid Parameters
params.Nsample = [2^8, 2^8];                % number of grid points in the sample grid
params.Nmap = [2^5, 2^5];                   % number of grid points saved in the map grid
params.k            = 0.5;                % Wave number
params.Lx           = 2 * pi / params.k;  % Length of spatial domain
params.Lv           = 12;                 % Velocity domain length (default)
params.Lv_s         = [12, 0.1 * pi];     % Velocity domain for each species
params.plot_freq = 10;                      % Nr. of iterations between plotting (if 0 no plotting)
params.measure_freq = 1;                    % Nr. of iterations between measurements (if 0 no measurements)
%% Species Parameters
params.Ns           = 1;                  % Number of species
params.species_name = ["electrons"];      % Names of species
params.Mr           = 1;                  % Mass ratio (if multiple species)
params.Mass         = [1];                % Mass of each species
params.charge       = [1];                % Charge of each species

%% Initial Condition Parameters
params.alpha        = 0.01;                % Perturbation amplitude
params.f0 = @(x,v) (1 + params.alpha * cos(x * params.k)) ...
                ./ sqrt(2 * pi) .* exp(-(v).^2 / 2);
params.fini = {params.f0};                % Initial distribution function

% Interpolation parameters
opts.scheme = 'lagrange-bary';
%opts.scheme = 'bspline';
opts.order = 3;
opts.use_mex = true;
params.opt_interp = opts;