function f_interp = lagrange_local_interp_periodic(x_target, xgrid, ygrid, order, method)
% Wrapper to choose interpolation method: 'normal' or 'barycentric'
    if nargin < 5
        method = 'normal';
    end

    switch lower(method)
        case 'normal'
            f_interp = lagrange_local_interp_periodic_standard(x_target, xgrid, ygrid, order);
        case 'barycentric'
            f_interp = lagrange_local_interp_barycentric(x_target, xgrid, ygrid, order);
        otherwise
            error('Unknown method "%s". Choose "normal" or "barycentric".', method);
    end
end





% -------------------------------------------------------------------------
% STANDARD LAGRANGE INTERPOLATION (slow and instable)
% -------------------------------------------------------------------------
function f_interp = lagrange_local_interp_periodic_standard(x_target, xgrid, ygrid, order)
% Interpolates y = f(x) at x_target using local Lagrange interpolation.
% xgrid: equispaced periodic grid (assumed 1D)
% ygrid: function values at xgrid
% order: number of points used in local stencil (e.g., 4 for cubic)
% x_target: vector of target points in [0,1)

    N = length(xgrid);
    delta_x = xgrid(2) - xgrid(1);  % uniform spacing assumed
    f_interp = zeros(size(x_target));
    half = floor(order/2);
    x_local = -half+(0:order);
    for i = 1:length(x_target)
        x = x_target(i);
        
        
        % Index of closest node
        idx0 = x/delta_x;
        j0 = floor(idx0);
        
        delta_idx = idx0 - j0;

        % Local grid and values
        idx_list = mod(j0 - half + (0:order), N) + 1;        
        y_local = ygrid(idx_list);

        % Evaluate local Lagrange interpolant
        p = 0;
        for j = 1:order+1
            Lj = lagrange_basis(delta_idx, x_local, j);
            p = p + y_local(j) * Lj;
        end

        f_interp(i) = p;
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





% -------------------------------------------------------------------------
% STANDARD LAGRANGE INTERPOLATION (slow and instable)
% -------------------------------------------------------------------------
function f_interp = lagrange_local_interp_barycentric(x_target, xgrid, ygrid, order)
    N = length(xgrid);
    delta_x = xgrid(2) - xgrid(1);
    f_interp = zeros(size(x_target));
    half = floor(order/2);
    x_local = -half+(0:order);
    w = barycentric_weights(x_local);

    for i = 1:length(x_target)
        x = x_target(i);
        idx0 = x/delta_x;
        j0 = floor(idx0);
        delta_idx = idx0 - j0;
        idx_list = mod(j0 - half + (0:order), N) + 1;
        y_local = ygrid(idx_list);
        f_interp(i) = barycentric_interp(delta_idx, x_local, y_local, w);
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