

function F_interp = lagrange1d_local_interp_periodic(x_target, xgrid, Fgrid, order, method)
% Wrapper to dispatch to 2D periodic Lagrange interpolation implementation
% method: 'normal' or 'barycentric'

    if nargin < 7
        method = 'barycentric';
    end

    switch lower(method)
        case 'normal'
            F_interp = lagrange1d_local_interp_periodic_normal(x_target, xgrid, Fgrid, order);
        case 'barycentric'
            F_interp = lagrange1d_local_interp_periodic_barycentric(x_target,  xgrid,  Fgrid, order)';
        case 'matlab_barycentric'
            F_interp = lagrange_local_interp_periodic(x_target, xgrid, Fgrid, order, "baycentric");
        case 'matlab_normal'
            F_interp = lagrange_local_interp_periodic(x_target,  xgrid, Fgrid, order, "normal");
        otherwise
            error('Unknown method "%s". Choose "normal" or "barycentric".', method);
    end
end