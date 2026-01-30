function [fs,params] = predictor_corrector_dt_adaptive(params,fs)

dt = params.dt;
dt_max = params.dt_interv(2);
dt_min = params.dt_interv(1);
Error_tolerance = params.dt_adapt_tolerance;
alpha = 0.9; % conservative factor for increasing the time step
iT = params.it + 1;

% Species indices
s_elec = find(params.species_name == "electrons");
s_ion  = find(params.species_name == "ions");

% Grid and physical parameters
charge = params.charge;
mass   = params.Mass;
grids  = params.grids;
iT = params.it+1;

% save the current distribution in case time step is rejected
fs_prev = fs;

% Compute electric field
[Efield] = vPoisson(fs, params.grids, params.charge);

% Advect distribution functions for half time step
for s = 1:params.Ns
    f12(:, :, s) = Advect(fs(:, :, s), params.charge(s) / params.Mass(s) * Efield, params.grids(s), params.dt / 2);
end


% Advect distribution functions for half time step
for s = 1:params.Ns
    fpred(:, :, s) = Advect(f12(:, :, s), params.charge(s) / params.Mass(s) * Efield, params.grids(s), params.dt / 2);
end

% Recompute electric field
[Efield] = vPoisson(f12, params.grids, params.charge);

% Advect distribution functions for full time step
for s = 1:params.Ns
    fs(:, :, s) = Advect(fs(:, :, s), params.charge(s) / params.Mass(s) * Efield, params.grids(s), params.dt);
end

Error = norm(fs(:) - fpred(:))/norm(fs(:));

if Error < Error_tolerance
    params.dt = min(dt_max, alpha*dt * sqrt((Error_tolerance/Error)));
else
    params.dt = max(dt_min, alpha*dt * sqrt((Error_tolerance/Error)));
    fs = predictor_corrector_dt_adaptive(params,fs_prev);
end

% save electric field;
params.Efield = Efield;
params.Efield_list(:,iT) = Efield;


end

