params = struct();
params.mycase = "keen_waves_weak";          % "keen_waves_weak"
params.Nx = 2^7;                            % Number of spatial grid points
params.Nv = 2^7;                            % Number of velocity grid points
params.Ns = 1;                              % Number of species (electrons only)
params.refine_factor=30;
params.v_range=[1.2, 1.6]
params.N_remap = 40;
params.method="CMM_vargrid";
%params.dt_adapt_tolerance = 0.1;            % realtive error allowed during time integration
%params.dt_interv = [1/100, 1];
params.species_name = ["electrons"]; % name of the different species
params.Mr = 1;                           % Mass ratio
params.Mass = [1];               % Mass of species
params.charge = [-1];                    % Charge of species
params.Nt_max = 200000;                    % Maximum number of time steps
params.dt = 0.5;                          % Time step size
params.dt_save = 200;                      % Save after dt_save time
params.Tend = 5000;                        % End time of simulation
params.plot_freq = 10;                      % Nr. of iterations between plotting (if 0 no plotting)
params.measure_freq = 1;                    % Nr. of iterations between measurements (if 0 no measurements)
% Keen waves parameters from paper (weak drive)
params.kDr = 0.26;                          % Drive wave number
params.wDr = 0.37;                          % Drive frequency
params.aDr = 0.00625;                       % Drive amplitude (weak)
params.TDr = 200;                           % Drive duration

% Drive timing parameters
params.t0 = 0;                              % Start time
params.tL = 69;                             % Left ramp start
params.twL = 20;                            % Left ramp width
params.twR = 20;                            % Right ramp width
params.tR = 207 + params.TDr;               % Right ramp start

% Grid parameters
params.Lx = 2 * pi / params.kDr;            % Spatial domain length
params.Lv = 6;                              % Velocity domain length

% Initial condition: spatially uniform Maxwellian
params.f0 = @(x, v) exp(-v.^2/2) / sqrt(2*pi);
params.fini = {params.f0};

% External drive field: E_Pond(x,t) = a_Dr * k_Dr * a(t) * sin(k_Dr * x - w_Dr * t)
params.g = @(t) 0.5 * (tanh((t - params.tL) / params.twL) - tanh((t - params.tR) / params.twR));
params.a = @(t) (params.g(t) - params.g(params.t0)) / (1 - params.g(params.t0));
params.E_ext = @(x, t) params.aDr * params.kDr * params.a(t) .* sin(params.kDr * x - params.wDr * t);

% Interpolation parameters
opts.scheme = 'lagrange-bary';
opts.order = 3;
opt.use_mex = true;
params.opt_interp = opts;
