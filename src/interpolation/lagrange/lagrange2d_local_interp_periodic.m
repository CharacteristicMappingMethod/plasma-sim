

function F_interp = lagrange2d_local_interp_periodic(x_target, y_target, xgrid, ygrid, Fgrid, order, method)
% Wrapper to dispatch to 2D periodic Lagrange interpolation implementation
% method: 'normal' or 'barycentric'

    if nargin < 7
        method = 'barycentric';
    end

    switch lower(method)
        case 'normal'
            F_interp = lagrange2d_local_interp_periodic_normal(x_target, y_target, xgrid, ygrid, Fgrid, order);
        case 'barycentric'
            F_interp = lagrange2d_local_interp_periodic_barycentric(x_target, y_target, xgrid, ygrid, Fgrid, order);
        otherwise
            error('Unknown method "%s". Choose "normal" or "barycentric".', method);
    end
end