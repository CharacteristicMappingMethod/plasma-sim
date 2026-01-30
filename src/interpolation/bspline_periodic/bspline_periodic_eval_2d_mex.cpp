#include "mex.h"
#include <vector>
#include <cmath>
#include <stdexcept>

// De Boor algorithm for cubic B-spline evaluation (p == 3)
double deboor_eval(int p, double x, const std::vector<double>& c) {
    // Only works for p == 3 (cubic B-spline)
    if (p != 3) {
        throw std::invalid_argument("Only cubic B-splines (p=3) are supported.");
    }

    // The knot vector and padding are chosen to match the MATLAB code
    std::vector<double> t = {-4, -4, -4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 4, 4, 4};
    int k = 7; // index of knot interval containing x (MATLAB 1-based)

    // Pad coefficients with zeros
    std::vector<double> c_pad = {0, 0, 0};
    c_pad.insert(c_pad.end(), c.begin(), c.end());
    c_pad.insert(c_pad.end(), {0, 0, 0});

    std::vector<double> d = {c_pad[k-3], c_pad[k-2], c_pad[k-1], c_pad[k]}; // 4 elements for cubic

    for (int r = 1; r <= 3; ++r) {
        for (int j = 3; j >= r; --j) {
            if (t[j+1+k-r] != t[j+k-3]) {
                double alpha = (x - t[j+k-3]) / (t[j+1+k-r] - t[j+k-3]);
                d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j];
            }
        }
    }

    return d[3]; // Return the final value
}

// 2D periodic B-spline evaluation
std::vector<double> bspline_periodic_eval_2d(
    const std::vector<double>& coeffs,
    int degree_x, int degree_y,
    const std::vector<double>& XQ,
    const std::vector<double>& YQ,
    int nx, int ny
) {
    double hx = 2*M_PI / nx;
    double hy = 2*M_PI / ny;
    int n_points = XQ.size();
    std::vector<double> vq(n_points, 0.0);

    // For cubic, spread is always 2 (covers 4 points)
    int spread_x = 2;
    int spread_y = 2;

    for (int idx = 0; idx < n_points; ++idx) {
        double x = fmod(XQ[idx], 2*M_PI);
        if (x < 0) x += 2*M_PI;
        double y = fmod(YQ[idx], 2*M_PI);
        if (y < 0) y += 2*M_PI;

        int ix = static_cast<int>(floor(x / hx));
        int iy = static_cast<int>(floor(y / hy));
        double xi = (x - ix*hx) / hx;
        double yi = (y - iy*hy) / hy;

        std::vector<double> local_coeffs_x(2*spread_x+1, 0.0);
        std::vector<double> temp(2*spread_y+1, 0.0);

        for (int j = -spread_y; j <= spread_y; ++j) {
            for (int i = -spread_x; i <= spread_x; ++i) {
                int idx_x = ((ix + i) % nx + nx) % nx; // Handle negative indices
                int idx_y = ((iy + j) % ny + ny) % ny;
                local_coeffs_x[i+spread_x] = coeffs[idx_y + idx_x * ny]; // Column-major indexing
            }
            temp[j+spread_y] = deboor_eval(degree_x, xi, local_coeffs_x);
        }

        vq[idx] = deboor_eval(degree_y, yi, temp);
    }

    return vq;
}

// MEX gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("bspline_periodic_eval_2d:nrhs",
                         "Five inputs required: coeffs, degree_x, degree_y, XQ, YQ");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("bspline_periodic_eval_2d:nlhs", "One output allowed.");
    }

    // Get input arguments
    double* coeffs_ptr = mxGetPr(prhs[0]);
    int degree_x = static_cast<int>(mxGetScalar(prhs[1]));
    int degree_y = static_cast<int>(mxGetScalar(prhs[2]));
    double* XQ_ptr = mxGetPr(prhs[3]);
    double* YQ_ptr = mxGetPr(prhs[4]);

    // Get dimensions
    const mwSize* coeffs_dims = mxGetDimensions(prhs[0]);
    int nx = static_cast<int>(coeffs_dims[0]);
    int ny = static_cast<int>(coeffs_dims[1]);
    int n_points = static_cast<int>(mxGetNumberOfElements(prhs[3]));

    // Check that XQ and YQ have the same size
    if (mxGetNumberOfElements(prhs[4]) != n_points) {
        mexErrMsgIdAndTxt("bspline_periodic_eval_2d:input",
                         "XQ and YQ must have the same number of elements.");
    }

    // Convert to std::vector
    std::vector<double> coeffs(coeffs_ptr, coeffs_ptr + nx * ny);
    std::vector<double> XQ(XQ_ptr, XQ_ptr + n_points);
    std::vector<double> YQ(YQ_ptr, YQ_ptr + n_points);

    // Call the evaluation function
    std::vector<double> vq;
    try {
        vq = bspline_periodic_eval_2d(coeffs, degree_x, degree_y, XQ, YQ, nx, ny);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("bspline_periodic_eval_2d:exception", e.what());
    }

    // Set output
    plhs[0] = mxCreateDoubleMatrix(n_points, 1, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);
    for (int i = 0; i < n_points; ++i) {
        out_ptr[i] = vq[i];
    }
}