function F_interp = interp1d_periodic(xq, xgrid, Fgrid, opts)
% Unified 1D periodic interpolation wrapper
%
% F_interp = interp1d_periodic(xq, xgrid, Fgrid, opts)
%
% Inputs:
%   xq       - query points (vector)
%   xgrid    - grid points (vector)
%   Fgrid    - function values on grid (vector)
%   opts     - (optional) struct with fields:
%       .scheme   - 'bspline' (default), 'lagrange-bary', 'lagrange-normal'
%       .order    - (for lagrange) interpolation order (default: 4)
%       .degree   - (for bspline) spline degree (default: 3)
%       .use_mex  - use MEX if available (default: true)
%       .cache    - use cache for bspline (default: true)
%
% Output:
%   F_interp - interpolated values at xq, same shape as xq

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'scheme'),    opts.scheme = 'lagrange-bary'; end
    if ~isfield(opts, 'order'),     opts.order = 3; end
    if ~isfield(opts, 'degree'),    opts.degree = 3; end
    if ~isfield(opts, 'use_mex'),   opts.use_mex = true; end
    if ~isfield(opts, 'cache'),     opts.cache = true; end

    % Flatten query points for evaluation
    sz = size(xq);
    xqv = xq(:);

    switch lower(opts.scheme)
        case 'bspline'
            % Use B-spline periodic interpolation
            % For now, use MATLAB's built-in spline with periodic extension
            % This provides proper 4th order convergence for cubic B-splines
            
            % Create periodic data by extending the grid and function values
            dx = xgrid(2)-xgrid(1);
            xgrid_extended = [xgrid(:); xgrid(end)+dx]; % Extend grid
            Fgrid_extended = [Fgrid(:); Fgrid(1)]; % Periodic extension
            
            % Use MATLAB's spline function
            F_interp = interp1(xgrid_extended, Fgrid_extended, mod(xqv,xgrid_extended(end)), 'spline');
        case 'lagrange-bary'
            % Use Lagrange barycentric periodic interpolation
            if opts.use_mex && (exist('lagrange1d_local_interp_periodic_barycentric.mexa64', 'file') == 2 || ...
                               exist('lagrange1d_local_interp_periodic_barycentric.mexmaci64', 'file') == 2)
                F_interp = lagrange1d_local_interp_periodic_barycentric(xqv, xgrid, Fgrid, opts.order);
            else
                F_interp = lagrange1d_local_interp_periodic(xqv, xgrid, Fgrid, opts.order, 'barycentric');
            end
        case 'lagrange-normal'
            % Use Lagrange normal periodic interpolation
            if opts.use_mex && (exist('lagrange1d_local_interp_periodic_normal.mexa64', 'file') == 2 || ...
                               exist('lagrange1d_local_interp_periodic_normal.mexmaci64', 'file') == 2)
                F_interp = lagrange1d_local_interp_periodic_normal(xqv, xgrid, Fgrid, opts.order);
            else
                F_interp = lagrange1d_local_interp_periodic(xqv, xgrid, Fgrid, opts.order, 'normal');
            end
    end
    % Reshape to match input
    F_interp = reshape(F_interp, sz);
end
