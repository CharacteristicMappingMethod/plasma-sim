params = struct();
params.mycase = "two_stream";          % "two_stream"
params.Nx = 2^8;                            % Number of spatial grid points
params.Nv = 2^8;                            % Number of velocity grid points
params.Ns = 1;                              % Number of species (electrons and ions)
params.method="predcorr";
%params.dt_adapt_tolerance = 0.1;            % realtive error allowed during time integration
%params.dt_interv = [1/100, 1];
params.species_name = ["ions"]; % name of the different species
params.Mr = 1;                           % Mass ratio for ions
params.Mass = [1];               % Mass of species
params.charge = [-1];                    % Charge of species
params.Nt_max = 4000;                       % Maximum number of time steps
params.dt = 1/4;                            % Time step size
params.dt_save = 10;                        % Save after dt_save time
params.Tend = 40;                           % End time of simulation

% Initial condition parameters
params.k = 0.2;                             % Wave number
params.eps = 5e-2;                         % Perturbation amplitude

% Electrons
params.v0 = 3;                   % Electron drift velocity
% Grid parameters
params.Lx = 2 * pi / params.k;    % Spatial domain length
params.Lv = 5 * pi; % Velocity domain lengths for each species

params.f0 = @(x, v) (1 + params.eps*cos(params.k*x)) ./(2*sqrt(2*pi)).*(exp(-(v-params.v0).^2/2)+exp(-(v+params.v0).^2/2));
params.fini = {params.f0};

