%

function [X,V] = sympl_flow(params, n, dt, X,V, ZEfield, grid)
if n == 1
    return;
end

data_size = size(V);
Uv = @(X,V,E) reshape(interp1d_periodic(X(:), grid.x, E(:), params.opt_interp), data_size);
% Omit the following line for Psi_tilda:
V = V + (dt / 2) * Uv(X,V,ZEfield(:,n));

while n > 2
    n = n - 1;
    X = X - dt * V;  % Inverse signs; going backwards in time
    V = V + dt *Uv(X,V,ZEfield(:,n));
end

X = X - dt * V;
V = V + (dt / 2) *Uv(X,V,ZEfield(:,1));

end


