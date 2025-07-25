function F_interp = interp2d_periodic(xq, yq, xgrid, ygrid, Fgrid, opts)
% Unified 2D periodic interpolation wrapper
%
% F_interp = interp2d_periodic(xq, yq, xgrid, ygrid, Fgrid, opts)
%
% Inputs:
%   xq, yq   - query points (vectors or meshgrid arrays)
%   xgrid, ygrid - grid points (vectors)
%   Fgrid    - function values on grid (size: [length(ygrid), length(xgrid)])
%   opts     - (optional) struct with fields:
%       .scheme   - 'bspline' (default), 'lagrange-bary', 'lagrange-normal'
%       .order    - (for lagrange) interpolation order (default: 4)
%       .degree   - (for bspline) spline degree (default: 3)
%       .use_mex  - use MEX if available (default: true)
%       .cache    - use cache for bspline (default: true)
%
% Output:
%   F_interp - interpolated values at (xq, yq), same shape as xq/yq

    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'scheme'),    opts.scheme = 'lagrange-bary'; end
    if ~isfield(opts, 'order'),     opts.order = 4; end
    if ~isfield(opts, 'degree'),    opts.degree = 3; end
    if ~isfield(opts, 'use_mex'),   opts.use_mex = true; end
    if ~isfield(opts, 'cache'),     opts.cache = true; end

    % Flatten query points for evaluation
    sz = size(xq);
    xqv = xq(:); yqv = yq(:);
    Nx = length(xgrid); Ny = length(ygrid);

    switch lower(opts.scheme)
        case 'bspline'
            % Use B-spline periodic interpolation
            if opts.use_mex && exist('bspline_periodic_eval_2d_mex_cpp', 'file') == 3
                interp = periodic_bspline_interpolant_2d(Fgrid, opts.degree, opts.degree, 'mex', opts.cache);
            else
                interp = periodic_bspline_interpolant_2d(Fgrid, opts.degree, opts.degree, 'matlab', opts.cache);
            end
            F_interp = interp.evaluate(xqv, yqv);
        case 'lagrange-bary'
            % Use Lagrange barycentric periodic interpolation
            if opts.use_mex && (exist('lagrange2d_local_interp_periodic_barycentric.mexa64', 'file') == 2 || ...
                               exist('lagrange2d_local_interp_periodic_barycentric.mexmaci64', 'file') == 2)
                F_interp = lagrange2d_local_interp_periodic_barycentric(xqv, yqv, xgrid, ygrid, Fgrid, opts.order);
            else
                F_interp = lagrange2d_local_interp_periodic(xqv, yqv, xgrid, ygrid, Fgrid, opts.order, 'barycentric');
            end
        case 'lagrange-normal'
            % Use Lagrange normal periodic interpolation
            if opts.use_mex && (exist('lagrange2d_local_interp_periodic_normal.mexa64', 'file') == 2 || ...
                               exist('lagrange2d_local_interp_periodic_normal.mexmaci64', 'file') == 2)
                F_interp = lagrange2d_local_interp_periodic_normal(xqv, yqv, xgrid, ygrid, Fgrid, opts.order);
            else
                F_interp = lagrange2d_local_interp_periodic(xqv, yqv, xgrid, ygrid, Fgrid, opts.order, 'normal');
            end
    end
    % Reshape to match input
    F_interp = reshape(F_interp, sz);
end