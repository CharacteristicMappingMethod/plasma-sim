%

function [X,V] = sympl_flow(n, dt, X,V, ZEfield, grid)
mint="spline";
if n == 1
    return;
end

size_grid = size(V);

periodic = @(x) mod(x,grid.Lx-grid.dx)

% Omit the following line for Psi_tilda:
V = V + (dt / 2) * reshape(interp1(grid.x,ZEfield(:,n),reshape(periodic(X),[],1),mint),size_grid);

while n > 2
    n = n - 1;
    X = X - dt * V;  % Inverse signs; going backwards in time
    V = V + dt *reshape(interp1(grid.x,ZEfield(:,n),reshape(periodic(X),[],1),mint),size_grid);
end

X = X - dt * V;
V = V + (dt / 2) *reshape(interp1(grid.x,ZEfield(:,1),reshape(periodic(X),[],1),mint),size_grid);

end


