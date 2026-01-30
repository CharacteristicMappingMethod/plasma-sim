
function F_interp = lagrange2d_local_interp_periodic_original(x_target, y_target, xgrid, ygrid, Fgrid, order, method)
% Wrapper to select between 'normal' and 'barycentric' 2D interpolation

    if nargin < 7
        method = 'normal';
    end

    switch lower(method)
        case 'normal'
            F_interp = lagrange2d_local_interp_periodic_normal(x_target, y_target, xgrid, ygrid, Fgrid, order);
        case 'barycentric'
            F_interp = lagrange2d_local_interp_barycentric(x_target, y_target, xgrid, ygrid, Fgrid, order);
        otherwise
            error('Unknown method "%s". Choose "normal" or "barycentric".', method);
    end
end



function F_interp = lagrange2d_local_interp_barycentric(x_target, y_target, xgrid, ygrid, Fgrid, order)

    Nx = length(xgrid); Ny = length(ygrid);
    dx = xgrid(2) - xgrid(1); dy = ygrid(2) - ygrid(1);

    x_target = x_target(:);
    y_target = y_target(:);
    Nt = length(x_target);
    F_interp = zeros(Nt, 1);

    half = floor(order/2);
    x_local = -half+(0:order);
    y_local = x_local;

    wx = barycentric_weights(x_local);
    wy = barycentric_weights(y_local);

    for i = 1:Nt
        x = x_target(i); y = y_target(i);
        idx_x = x/dx; jx = floor(idx_x); delta_x = idx_x - jx;
        idxs_x = mod(jx - half + (0:order), Nx) + 1;

        idx_y = y/dy; jy = floor(idx_y); delta_y = idx_y - jy;
        idxs_y = mod(jy - half + (0:order), Ny) + 1;

        fx = zeros(order+1, 1);
        for jj = 1:(order+1)
            row = Fgrid(idxs_y(jj), idxs_x);
            fx(jj) = barycentric_interp(delta_x, x_local, row, wx);
        end

        F_interp(i) = barycentric_interp(delta_y, y_local, fx', wy);
    end
end

function w = barycentric_weights(x_nodes)
    n = length(x_nodes);
    w = ones(1, n);
    for j = 1:n
        for k = [1:j-1, j+1:n]
            w(j) = w(j) / (x_nodes(j) - x_nodes(k));
        end
    end
end

function p = barycentric_interp(x, x_nodes, y_nodes, w)
    diff = x - x_nodes;
    zero_idx = find(abs(diff) < 1e-14, 1);
    if ~isempty(zero_idx)
        p = y_nodes(zero_idx);
        return;
    end
    num = sum(w .* y_nodes ./ diff);
    denom = sum(w ./ diff);
    p = num / denom;
end




function F_interp = lagrange2d_local_interp_periodic_normal(x_target, y_target, xgrid, ygrid, Fgrid, order)
% 2D local Lagrange interpolation with periodic boundaries.
% xgrid, ygrid: 1D periodic grids
% Fgrid: function values on tensor product grid (size: Ny x Nx)
% x_target, y_target: vectors of target points in [0,1)
% order: interpolation order (e.g. 4 for cubic -> 5 stencil points)

    Nx = length(xgrid);
    Ny = length(ygrid);
    dx = xgrid(2) - xgrid(1);
    dy = ygrid(2) - ygrid(1);

    % Ensure column vectors
    x_target = x_target(:);
    y_target = y_target(:);
    Nt = length(x_target);
    F_interp = zeros(Nt, 1);

    half = floor(order/2);
    local_nodes = -half+(0:order); % Relative stencil positions

    for i = 1:Nt
        % Target point
        x = x_target(i);
        y = y_target(i);

        % Index and local x
        idx_x = x/dx;
        jx = floor(idx_x);
        delta_x = idx_x - jx;
        idxs_x = mod(jx - half + (0:order), Nx) + 1;
        x_local = local_nodes;

        % Index and local y
        idx_y = y/dy;
        jy = floor(idx_y);
        delta_y = idx_y - jy;
        idxs_y = mod(jy - half + (0:order), Ny) + 1;
        y_local = local_nodes;

        % Interpolate in x for each row in y
        fx = zeros(order+1, 1);
        for jj = 1:(order+1)
            row = Fgrid(idxs_y(jj), idxs_x);
            fx(jj) = 0;
            for j = 1:(order+1)
                Ljx = lagrange_basis(delta_x, x_local, j);
                fx(jj) = fx(jj) + row(j) * Ljx;
            end
        end

        % Interpolate result in y
        fy = 0;
        for j = 1:(order+1)
            Ljy = lagrange_basis(delta_y, y_local, j);
            fy = fy + fx(j) * Ljy;
        end

        F_interp(i) = fy;
    end
end

function Lj = lagrange_basis(x, x_nodes, j)
% Computes the j-th Lagrange basis polynomial evaluated at x
    Lj = 1;
    for m = 1:length(x_nodes)
        if m ~= j
            Lj = Lj * (x - x_nodes(m)) / (x_nodes(j) - x_nodes(m) + 1e-32);
        end
    end
end
