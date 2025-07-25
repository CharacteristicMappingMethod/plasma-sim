function interp = periodic_bspline_interpolant_2d(u, degree_x, degree_y, method, use_cache)
% u: 2D data values (size: [nx, ny])
% degree_x: degree of the B-spline in x
% degree_y: degree of the B-spline in y
% method: 'matlab', 'mex', or 'auto' (default: 'matlab')
% use_cache: (optional, default true) if true, use persistent cache for divisor FFT
% Returns a struct with fields:
%   coeffs, degree_x, degree_y, nx, ny, evaluate (function handle), method

    if nargin < 4 || isempty(method)
        method = 'matlab';
    end
    method = lower(method);
    if strcmp(method, 'auto')
        if exist('bspline_periodic_eval_2d_mex_cpp', 'file') == 3
            method = 'mex';
        else
            method = 'matlab';
        end
    end
    if nargin < 5
        use_cache = true;
    end

    [nx, ny] = size(u);
    dx = 2*pi / nx;
    dy = 2*pi / ny;
    x = linspace(0, 2*pi - dx, nx);
    y = linspace(0, 2*pi - dy, ny);
    bspline_divisor = zeros(nx, ny);
    bspline_divisor(1,1) = 1;
    [Xg, Yg] = ndgrid(x, y); % Use ndgrid for 'ij' indexing

    persistent divisor_fft_cache;
    if isempty(divisor_fft_cache)
        divisor_fft_cache = containers.Map();
    end
    cache_key = sprintf('nx%d_ny%d_dx%d_dy%d_degx%d_degy%d', nx, ny, degree_x, degree_y, degree_x, degree_y);
    if use_cache && isKey(divisor_fft_cache, cache_key)
        bspline_divisor_eval_fft = divisor_fft_cache(cache_key);
    else
        bspline_divisor_eval = bspline_periodic_eval_2d(bspline_divisor, degree_x, degree_y, Xg, Yg);
        bspline_divisor_eval_fft = fft2(bspline_divisor_eval);
        if use_cache
            divisor_fft_cache(cache_key) = bspline_divisor_eval_fft;
        end
    end
    coeffs = real(ifft2(fft2(u) ./ bspline_divisor_eval_fft));

    interp.coeffs = coeffs;
    interp.degree_x = degree_x;
    interp.degree_y = degree_y;
    interp.nx = nx;
    interp.ny = ny;
    interp.method = method;

    if strcmp(method, 'mex')
        interp.evaluate = @(YQ, XQ) bspline_periodic_eval_2d_mex(coeffs, degree_x, degree_y, XQ, YQ);
    else
        interp.evaluate = @(YQ, XQ) bspline_periodic_eval_2d(coeffs, degree_x, degree_y, XQ, YQ);
    end
end

function vq = bspline_periodic_eval_2d(coeffs, degree_x, degree_y, XQ, YQ)
    [nx, ny] = size(coeffs);
    hx = 2*pi / nx;
    hy = 2*pi / ny;
    vq = zeros(size(XQ));
    for idx = 1:numel(XQ)
        x = mod(XQ(idx), 2*pi);
        y = mod(YQ(idx), 2*pi);
        ix = floor(x / hx);
        iy = floor(y / hy);
        xi = (x - ix*hx) / hx;
        yi = (y - iy*hy) / hy;
        % For cubic, spread is always 2 (covers 4 points)
        spread_x = 2;
        spread_y = 2;
        local_coeffs_x = zeros(1, 2*spread_x+1);
        temp = zeros(1, 2*spread_y+1);
        for j = -spread_y:spread_y
            for i = -spread_x:spread_x
                idx_x = mod(ix + i, nx) + 1; % MATLAB 1-based
                idx_y = mod(iy + j, ny) + 1;
                local_coeffs_x(i+spread_x+1) = coeffs(idx_x, idx_y);
            end
            temp(j+spread_y+1) = deboor_eval(degree_x, xi, local_coeffs_x);
        end
        vq(idx) = deboor_eval(degree_y, yi, temp);
    end
end

function val = deboor_eval(p, x, c)
    % Only works for p == 3 (cubic B-spline)
    assert(p == 3, 'Only cubic B-splines (p=3) are supported.');
    % The knot vector and padding are chosen to match the Python code
    t = [-4, -4, -4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 4, 4, 4];
    k = 8; % index of knot interval containing x (MATLAB 1-based)
    c_pad = [0, 0, 0, c, 0, 0, 0];
    d = c_pad(k-3:k); % 4 elements for cubic
    for r = 1:3
        for j = 3:-1:r
            if t(j+1+k-r) ~= t(j+k-3)
                alpha = (x - t(j+k-3)) / (t(j+1+k-r) - t(j+k-3));
                d(j+1) = (1.0 - alpha) * d(j) + alpha * d(j+1);
            end
        end
    end
    val = d(4);
end


function vq = bspline_periodic_eval_2d_mex(coeffs, degree_x, degree_y, XQ, YQ)
    % MEX wrapper for 2D periodic B-spline evaluation
    % coeffs: 2D coefficient array [nx, ny]
    % degree_x, degree_y: B-spline degrees (should be 3 for cubic)
    % XQ, YQ: query points (can be arrays of any size)
    % Returns: interpolated values at query points

        % Ensure XQ and YQ are column vectors
        XQ = XQ(:);
        YQ = YQ(:);

        % Call the MEX function
        vq = bspline_periodic_eval_2d_mex_cpp(coeffs, degree_x, degree_y, YQ, XQ);

        % Reshape output to match input shape if XQ/YQ were not column vectors
        if ~iscolumn(XQ) || ~iscolumn(YQ)
            vq = reshape(vq, size(XQ));
        end
    end