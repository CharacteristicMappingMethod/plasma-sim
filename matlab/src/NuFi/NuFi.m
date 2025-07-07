function [fs, params] = NuFi(params,fs)
iT = params.it+1;
dt = params.dt;
for s = 1:params.Ns
    [X,V] = sympl_flow_Half(iT,dt,params.grids(s).X,params.grids(s).V,params.charge(s)/params.Mass(s)*params.Efield_list,params.grids(s));
    fini = params.fini{s};
    fs(:,:,s) = fini(X,V);
end
[Efield] =vPoisson(fs,params.grids,params.charge);
params.Efield = Efield;
params.Efield_list(:,iT) = Efield;

end






function [X, V] = sympl_flow_Half(n, dt, X, V, Efield, grid)
mint="lagrange";
order = 4; % lagrange interpolation order
if n == 1
    return;
end


periodic = @(x) mod(x,grid.Lx-grid.dx);

while n > 2
    n = n - 1;
    X = X - dt * V;  % Inverse signs; going backwards in time
    if mint=="lagrange"
        V = V + dt *reshape(lagrange1d_local_interp_periodic(reshape(X,[],1),grid.x,Efield(:,n),order),grid.size);
    else
        V = V + dt *reshape(interp1(grid.x,Efield(:,n),reshape(periodic(X),[],1),mint),grid.size);
    end

end

X = X - dt * V;
if mint=="lagrange"
    V = V + (dt / 2) *reshape(lagrange1d_local_interp_periodic(reshape(X,[],1),grid.x,Efield(:,1),order),grid.size);
else
    V = V + (dt / 2) *reshape(interp1(grid.x,Efield(:,1),reshape(periodic(X),[],1),mint),grid.size);
end

end

