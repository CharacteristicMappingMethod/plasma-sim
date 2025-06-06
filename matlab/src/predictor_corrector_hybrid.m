function [fs,params] = predictor_corrector_hybrid(params,fs)
dt = params.dt;
s_elec = find(params.species_name=="electrons");
s_ion =  find(params.species_name=="ions");

charge = params.charge;

% ion density
n_ions = charge(s_ion)*sum(fs(:,:,s_ion))*params.grids(s_ion).dv;

% electron density
n_elec = charge(s_elec)*sum(fs(:,:,s_elec))*params.grids(s_elec).dv;

iT = params.it+1;

% Compute electric field
rho = n_ions + n_elec;
[Efield,phi_old] = Poisson(rho, params.grids);

%% advect ions
f12(:, :, s_ion) = Advect(fs(:, :, s_ion), params.charge(s_ion) / params.Mass(s_ion) * Efield, params.grids(s_ion), params.dt / 2);

% Recompute electric field
n_ions = charge(s_ion)*sum(fs(:,:,s_ion))*params.grids(s_ion).dv;
rho = n_ions + n_elec;
[Efield,phi] = Poisson(rho, params.grids);

% Advect distribution functions for full time step
fs(:, :, s_ion) = Advect(fs(:, :, s_ion), params.charge(s_ion) / params.Mass(s_ion) * Efield, params.grids(s_ion), params.dt);

%% advect electrons
% fluid part:

T_e = 10;
X = params.grids(s_elec).X;
v = params.grids(s_elec).v;
[Phi,V] = meshgrid(phi,v);
n_elec_fluid = (1) .*exp(charge(s_elec)*Phi/T_e);

f_elec_fluid =@(X,V,Phi_field) params.fe0(X,V).*exp(charge(s_elec)*Phi_field/T_e);

% kinetic part
dg = fs(:,:,s_elec) - f_elec_fluid(X,V,Phi);
dg = dg + 0.5 * dt * source_term(params,dg,charge(s_elec)*f_elec_fluid(X,V,Phi)/T_e, phi, phi_old, -Efield);
% left hand side:
dg = Advect(dg, params.charge(s_elec) / params.Mass(s_elec) * Efield, params.grids(s_elec), params.dt);
% right hand side:
fhalf =  f_elec_fluid(X,V,Phi)+ dg;
n_ions = charge(s_ion)*sum(fs(:,:,s_ion))*params.grids(s_ion).dv;
n_elec = charge(s_elec)*sum(fhalf)*params.grids(s_elec).dv;
rho = n_ions + n_elec;

[Efield,phi] = Poisson(rho, params.grids);
[Phi,V] = meshgrid(phi,v);


dg = Advect(dg, params.charge(s_elec) / params.Mass(s_elec) * Efield, params.grids(s_elec), params.dt);
dg = dg + 0.5 * dt * source_term(params,dg,charge(s_elec)*  f_elec_fluid(X,V,Phi)/T_e, phi, phi_old, -Efield);

%n_elec_kinetic =
fs(:,:,s_elec) = f_elec_fluid(X,V,Phi) + dg;

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