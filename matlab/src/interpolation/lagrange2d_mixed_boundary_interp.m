function F_interp = lagrange2d_mixed_boundary_interp(x_target, y_target, xgrid, ygrid, Fgrid, order)
% 2D local Lagrange interpolation with PERIODIC boundaries in X 
% and CLAMPED (non-periodic) boundaries in Y.

    Nx = length(xgrid);
    Ny = length(ygrid);
    dx = xgrid(2) - xgrid(1);
    dy = ygrid(2) - ygrid(1);
    
    x_target = x_target(:);
    y_target = y_target(:);
    Nt = length(x_target);
    F_interp = zeros(Nt, 1);
    
    half = floor(order/2);
    local_nodes = -half+(0:order);

    for i = 1:Nt
        x = x_target(i);
        y = y_target(i);
        
        % --- Periodic X dimension (unchanged) ---
        idx_x = x/dx;
        jx = floor(idx_x);
        delta_x = idx_x - jx;
        idxs_x = mod(jx - half + (0:order), Nx) + 1;
        
        % --- Non-Periodic Y dimension (CORRECTED) ---
        % Find the nearest grid index for y
        jy = round((y - ygrid(1)) / dy) + 1; 
        
        % Get stencil indices, but clamp them to stay within the grid
        stencil_y_indices = jy - half + (0:order);
        idxs_y = max(1, min(Ny, stencil_y_indices)); % Clamp to [1, Ny]
        
        % The target y relative to the *local* stencil's grid points
        delta_y = (y - ygrid(idxs_y(half+1))) / dy;

        % --- Interpolation (modified for new delta_y) ---
        fx = zeros(order+1, 1);
        y_local = local_nodes;
        
        % Interpolate in x for each row in y-stencil
        for jj = 1:(order+1)
            row = Fgrid(idxs_y(jj), idxs_x);
            fx(jj) = lagrange_basis_eval(delta_x, local_nodes, row);
        end
        
        % Interpolate result in y
        fy = lagrange_basis_eval(delta_y, y_local, fx');
        
        F_interp(i) = fy;
    end
end

% Helper function to make the main loop cleaner
function val = lagrange_basis_eval(x, nodes, values)
    val = 0;
    for j = 1:length(nodes)
        Lj = 1;
        for m = 1:length(nodes)
            if m ~= j
                Lj = Lj * (x - nodes(m)) / (nodes(j) - nodes(m) + 1e-32);
            end
        end
        val = val + values(j) * Lj;
    end
end