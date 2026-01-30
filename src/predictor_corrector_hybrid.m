function [fs,params] = predictor_corrector_hybrid(params,fs)

iT = params.it+1;
dt = params.dt;
iT = params.it + 1;

% Species indices
s_elec = find(params.species_name == "electrons");
s_ion  = find(params.species_name == "ions");

% Grid and physical parameters
charge = params.charge;
mass   = params.Mass;
grids  = params.grids;
Nv = params.Nv;
T_e = 1;

v_th = sqrt(T_e/mass(s_elec));
maxw = @(v) 1/sqrt(2*pi*v_th.^2) * exp(-v.^2/(2*v_th^2));
fadiabatic =@(Phi_field) maxw(grids(s_elec).V).*exp(charge(s_elec)*Phi_field/T_e);

% Compute electric field
[Efield,phi] = vPoisson(fs, params.grids, params.charge);
Phi = repmat(phi(:), 1, Nv);

% compute effective potential
ne = sum(fs(:,:,s_elec))*params.grids(s_elec).dv;
n0 = sum(maxw(grids(s_elec).V)*grids(s_elec).dv);
phi_eff = log(ne./n0)*T_e/charge(s_elec);
Phi_eff = repmat(phi_eff(:), 1, Nv);

% Advect distribution functions for half time step
f12(:, :, s_ion) = Advect(fs(:, :, s_ion), params.charge(s_ion) / params.Mass(s_ion) * Efield, params.grids(s_ion), params.dt / 2);

df = fs(:, : ,s_elec) - fadiabatic(Phi_eff);
df12 = Advect(df, params.charge(s_elec) / params.Mass(s_elec) * Efield, params.grids(s_elec), params.dt / 2);

f12(:, : ,s_elec) = fadiabatic(Phi);


% Recompute electric field
[Efield,phi] = vPoisson(f12, params.grids, params.charge);
Phi = repmat(phi(:), 1, Nv);

% Advect distribution functions for full time step
fs(:, :, s_ion) = Advect(fs(:, :, s_ion), params.charge(s_ion) / params.Mass(s_ion) * Efield, params.grids(s_ion), params.dt);
fs(:, : ,s_elec) = fadiabatic(Phi);

% save electric field;
params.Efield = Efield;
params.Efield_list(:,iT) = Efield;


end


function [Efield,phi] = vPoisson(fs, grids, charge)

% compute charge density
rho = 0;
Ns = length(grids);
for s = 1:Ns
    rho = rho + charge(s)*sum(fs(:,:,s))*grids(s).dv;
end

Nx = grids(1).Nx;
Lx = grids(1).Lx;
kx = (2*pi/Lx) * [-Nx/2:Nx/2-1];
kx = fftshift(kx);

% laplacian is division -|k|^2
K2 = kx.^2;
K2(1) = 1; % avoid devision by 0
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