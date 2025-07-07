% Initialize parameter struct
params = struct();

%% General Simulation Parameters
params.mycase       = "landau_damping";   % Options: "landau_damping", "two_stream"
params.method       = "predcorr";         % Time integration method
params.Nt_max       = 4000;               % Max number of time steps
params.dt           = 0.5;                % Time step size
params.Tend         = 80;                 % Simulation end time
params.dt_save      = 10;                 % Save data every dt_save time units

%% Grid Parameters
params.Nx           = 2^8;                % Number of spatial grid points
params.Nv           = 2^8;                % Number of velocity grid points
params.k            = 0.5;                % Wave number
params.Lx           = 2 * pi / params.k;  % Length of spatial domain
params.Lv           = 12;                 % Velocity domain length (default)
params.Lv_s         = [12, 0.1 * pi];     % Velocity domain for each species

%% Species Parameters
params.Ns           = 1;                  % Number of species
params.species_name = ["electrons"];      % Names of species
params.Mr           = 1;                  % Mass ratio (if multiple species)
params.Mass         = [1];                % Mass of each species
params.charge       = [1];                % Charge of each species

%% Initial Condition Parameters
params.alpha        = 0.05;                % Perturbation amplitude
params.f0 = @(x,v) (1 + params.alpha * cos(x * params.k)) ...
                ./ sqrt(2 * pi) .* exp(-(v).^2 / 2);
params.fini = {params.f0};                % Initial distribution function
