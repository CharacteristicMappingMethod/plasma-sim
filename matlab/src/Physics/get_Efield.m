function E = get_Efield(params, Efield, x, optional_time)
% Interpolate Efield array and add external ponderomotive Efield
%
% Inputs:
%   params  - simulation parameters structure
%   Efield  - electric field array on the grid
%   x       - query points (can be array)
%
% Output:
%   E       - interpolated electric field at query points
%
% The function uses interp1d_periodic for interpolation and adds
% an external ponderomotive field if params.E_ext is defined

    % Get current time
    if nargin < 4
        t = params.time_array(end);
    else
        t = optional_time;
    end
    
    % Interpolate the main Efield using interp1d_periodic
    E = interp1d_periodic(x, params.grids(1).x, Efield, params.opt_interp);
    
    % Add external ponderomotive Efield if specified
    E_tot = E + compute_external_Efield(params, x, t);
    
    % Ensure output has the same shape as input x
    E = reshape(E, size(x));
end
