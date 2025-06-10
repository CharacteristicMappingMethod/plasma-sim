function [fs,params] = predictor_corrector_hybrid(params,fs)
% Time step
dt = params.dt;
iT = params.it + 1;

% Species indices
s_elec = find(params.species_name == "electrons");
s_ion  = find(params.species_name == "ions");

% Grid and physical parameters
charge = params.charge;
mass   = params.Mass;
grids  = params.grids;
T_e    = 100;  % Electron temperature (assumed constant)

% Spatial and velocity grids for electrons
X = grids(s_elec).X;
v = grids(s_elec).v;

%% --- Step 1: Compute initial charge density and electric field ---
% Ion density (charge-weighted)
n_ions = charge(s_ion) * sum(fs(:,:,s_ion), 2) * grids(s_ion).dv;

% Electron density (charge-weighted)
n_elec = charge(s_elec) * sum(fs(:,:,s_elec), 2) * grids(s_elec).dv;

% Total charge density
rho = n_ions + n_elec;

% Solve Poisson equation to get E and φ at t^n
[Efield, phi_old] = Poisson(rho, grids);

%% --- Step 2: Advect ions (first half step) ---
accel_ion = charge(s_ion) / mass(s_ion) * Efield;
f12(:,:,s_ion) = Advect(fs(:,:,s_ion), accel_ion, grids(s_ion), dt/2);

%% --- Step 3: Recompute E field after half-step ion push ---
n_ions = charge(s_ion) * sum(f12(:,:,s_ion), 2) * grids(s_ion).dv;
rho = n_ions + n_elec;
[Efield, phi] = Poisson(rho, grids);

%% --- Step 4: Advect ions full step with updated E field ---
fs(:,:,s_ion) = Advect(fs(:,:,s_ion), charge(s_ion)/mass(s_ion)*Efield, grids(s_ion), dt);

%% --- Step 5: Electron fluid–kinetic update ---
% Electron fluid density from Boltzmann response
[Phi, V] = meshgrid(phi, v);
n0 = 1;  % Reference density (can be computed from initial conditions)
f_elec_fluid = @(X, V, Phi_field) params.fe0(X, V) .* exp(charge(s_elec) * Phi_field / T_e);

% Compute kinetic part: δg = f - f_fluid
f_fluid = f_elec_fluid(X, V, Phi);
dg = fs(:,:,s_elec) - f_fluid;

% First half source step
S_half_1 = source_term(params, dg, charge(s_elec) * f_fluid / T_e, phi, phi_old, -Efield);
dg = dg + 0.5 * dt * S_half_1;

% Advection step (semi-Lagrangian)
accel_elec = charge(s_elec) / mass(s_elec) * Efield;
dg = Advect(dg, accel_elec, grids(s_elec), dt);

% Recompute φ and E for next half source step
f_half = f_elec_fluid(X, V, Phi) + dg;
n_elec = charge(s_elec) * sum(f_half, 2) * grids(s_elec).dv;
n_ions = charge(s_ion)  * sum(fs(:,:,s_ion), 2) * grids(s_ion).dv;
rho = n_ions + n_elec;

[Efield, phi] = Poisson(rho, grids);
[Phi, V] = meshgrid(phi, v);  % Update meshgrid

% Second half source step
f_fluid = f_elec_fluid(X, V, Phi);  % Update f_fluid with new φ
S_half_2 = source_term(params, dg, charge(s_elec) * f_fluid / T_e, phi, phi_old, -Efield);
dg = dg + 0.5 * dt * S_half_2;

% Final electron distribution: fluid + kinetic
fs(:,:,s_elec) = f_fluid + dg;

%% update e field
% save electric field;
params.Efield = Efield;
params.Efield_list(:,iT) = Efield;


end

function [Efield,phi] = Poisson(rho, grids)


Nx = grids(1).Nx;
Lx = grids(1).Lx;
kx = (2*pi/Lx) * [-Nx/2:Nx/2-1];
kx = fftshift(kx);

% laplacian is division -|k|^2
K2 = kx.^2;
K2(1) = 1; % avoid devision by 0params.fe0(X,V)
% to avoid a division by zero, we set the zeroth wavenumber to one.
% this leaves it's respective Fourier coefficient unaltered, so the
% zero mode of Sk is conserved.dphi_dx_h = 1i*phi_fft.*kx(1,:); This way, Sk's zero mode implicitly
% defined the zero mode of the result
% Note that the zero mode is NOT uniquely defined: in a periodic
% setting, the solution of Laplace's (or Poisson's) equation is only
% defined up to a constant! You can freely overwrite the zero mode,
% therefore.
b = fft(rho);
phi_fft =  -b ./(K2); % solves second equation of vlassov poisson
phi_fft(1) = 0; % set mean to zero
dphi_dx_h = 1i*phi_fft.*kx;
Efield = -reshape(ifft(dphi_dx_h, "symmetric"), 1, []);
phi = reshape(ifft(phi_fft, "symmetric"), 1, []);
end

function S = source_term(params,dg, fe0, phi, phi_old, dphi_dx)
    % Approximate d(phi)/dt
    dphi_dt = (phi - phi_old) / params.dt;

    

    % Source term
    S = -fe0 .* (dphi_dt + params.grids(1).V .* dphi_dx);
end