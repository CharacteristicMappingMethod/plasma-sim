function [detJ, dXdx, dXdv, dVdx, dVdv] = jacobian_determinant(X, V, grid, method)
% JACOBIAN_DETERMINANT Compute Jacobian determinant of map [X,V] using Fourier transforms
%
% This function computes the Jacobian determinant of the coordinate transformation
% (x,v) -> (X,V) using Fourier spectral differentiation for high accuracy.
%
% Inputs:
%   X, V  - Transformed coordinate arrays [Nv, Nx] (from CMM.m)
%   grid  - Grid structure containing coordinate information
%   method - Method to use for differentiation ("fourier" or "FD")
%
% Outputs:
%   detJ  - Jacobian determinant array [Nv, Nx]
%   dXdx  - Partial derivative ∂X/∂x [Nv, Nx]
%   dXdv  - Partial derivative ∂X/∂v [Nv, Nx] 
%   dVdx  - Partial derivative ∂V/∂x [Nv, Nx]
%   dVdv  - Partial derivative ∂V/∂v [Nv, Nx]
%
% The Jacobian determinant is computed as:
%   detJ = (∂X/∂x)(∂V/∂v) - (∂X/∂v)(∂V/∂x)

if (nargin < 4 || isempty(method))
    method = "fourier";
end

    % Get grid dimensions
    Nx = grid.Nx;
    Nv = grid.Nv;
    
    % Get coordinate arrays
    x = grid.x;  % [1, Nx] - spatial coordinates
    v = grid.v;  % [1, Nv] - velocity coordinates
    
    % Ensure X and V are the same size as the grid
    if ~isequal(size(X), [Nv, Nx]) || ~isequal(size(V), [Nv, Nx])
        error('X and V must have size [Nv, Nx] = [%d, %d]', Nv, Nx);
    end
    
    % Compute partial derivatives using Fourier spectral differentiation
    if method == "fourier"
        % 1. ∂X/∂x: differentiate X with respect to x (along columns)
        dXdx = fourier_derivative_1d(X-grid.X, grid.Lx, 2);
        
        % 2. ∂X/∂v: differentiate X with respect to v (along rows)  
        dXdv = fourier_derivative_1d(X-grid.X, 2*grid.Lv, 1);
        
        % 3. ∂V/∂x: differentiate V with respect to x (along columns)
        dVdx = fourier_derivative_1d(V-grid.V, grid.Lx, 2);
        
        % 4. ∂V/∂v: differentiate V with respect to v (along rows)
        dVdv = fourier_derivative_1d(V-grid.V, 2*grid.Lv, 1);
        
        % Compute Jacobian determinant
        detJ = (dXdx+1) .* (dVdv+1) - dXdv .* dVdx;
        fprintf("max(abs(detJ(:)-1)) = %3.2e\n", max(abs(detJ(:)-1)))
    elseif method == "FD"
        % 1. ∂X/∂x: differentiate X with respect to x (along columns)
        dXdx = grid.Dx*X;
                
        % 2. ∂X/∂v: differentiate X with respect to v (along rows)  
        dXdv = X*grid.Dv;
        
        % 3. ∂V/∂x: differentiate V with respect to x (along columns)
        dVdx = grid.Dx*V;
        
        % 4. ∂V/∂v: differentiate V with respect to v (along rows)
        dVdv = V*grid.Dv;
        
        % Compute Jacobian determinant
        detJ = (dXdx) .* (dVdv) - dXdv .* dVdx;
        detJ = detJ(2:end-1,2:end-1);
        fprintf("max(abs(detJ2(:)-1)) = %3.2e\n", max(abs(detJ(:)-1)))
    end
    
    % Optional: Check for negative determinants (non-physical)
    if any(detJ(:) < 0)
     %   warning('Negative Jacobian determinants detected! Map may not be orientation-preserving.');
    end
    
    % Optional: Check for very small determinants (near-singular map)
    min_det = min(detJ(:));
    if min_det < 1e-10
      %  warning('Very small Jacobian determinants detected (min = %.2e). Map may be near-singular.', min_det);
    end
    
end

function df = finite_difference(f, L, dim)

    % Get size of f
    sz = size(f);
    
    if dim == 1
        % Differentiate along first dimension (rows)
        df = Dx*f;

    elseif dim == 2
        % Differentiate along second dimension (columns)
        N = sz(2);
        % Create wavenumber array for this dimension
        k = (2*pi/L) * [-N/2:N/2-1];
        k = fftshift(k);
        % Compute derivative
        f_hat = fft(f, [], 2);
        df_hat = f_hat .* (1i * k);
        df = ifft(df_hat, [], 2, 'symmetric');
        
    else
        error('dim must be 1 or 2');
    end
    
end



function df = fourier_derivative_1d(f, L, dim)
% FOURIER_DERIVATIVE_1D Compute 1D derivative using Fourier spectral method
%
% Inputs:
%   f     - Function values [Nv, Nx] or [Nx, Nv]
%   coord - Coordinate array (x or v)
%   dim   - Dimension along which to differentiate (1 or 2)
%
% Output:
%   df    - Derivative    % Get coordinate spacing

    
    % Get size of f
    sz = size(f);
    
    if dim == 1
        % Differentiate along first dimension (rows)
        N = sz(1);
        % Create wavenumber array for this dimension
        k = (2*pi/L) * [-N/2:N/2-1];
        k = fftshift(k);
        % Compute derivative
        f_hat = fft(f, [], 1);
        df_hat = 1i * k(:) .* f_hat;
        df = ifft(df_hat, [], 1, 'symmetric');
        
    elseif dim == 2
        % Differentiate along second dimension (columns)
        N = sz(2);
        % Create wavenumber array for this dimension
        k = (2*pi/L) * [-N/2:N/2-1];
        k = fftshift(k);
        % Compute derivative
        f_hat = fft(f, [], 2);
        df_hat = f_hat .* (1i * k);
        df = ifft(df_hat, [], 2, 'symmetric');
        
    else
        error('dim must be 1 or 2');
    end
    
end

