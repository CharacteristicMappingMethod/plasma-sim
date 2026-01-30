#include "mex.h"
#include <vector>
#include <cmath>
#include <stdexcept>

// Helper: computes the j-th Lagrange basis polynomial at x
double lagrange_basis(double x, const std::vector<double>& x_nodes, int j) {
    double Lj = 1.0;
    int n = x_nodes.size();
    for (int m = 0; m < n; ++m) {
        if (m != j) {
            Lj *= (x - x_nodes[m]) / (x_nodes[j] - x_nodes[m] + 1e-32);
        }
    }
    return Lj;
}

// Main interpolation function
std::vector<double> lagrange_local_interp_periodic(
    const std::vector<double>& x_target,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid,
    int order
) {
    int N = xgrid.size();
    if (N < order + 1) throw std::invalid_argument("xgrid too small for given order");
    double delta_x = xgrid[1] - xgrid[0]; // uniform spacing assumed
    std::vector<double> f_interp(x_target.size(), 0.0);
    int half = order / 2;

    // Build local x stencil offsets
    std::vector<double> x_local(order + 1);
    for (int k = 0; k <= order; ++k)
        x_local[k] = -half + k;

    for (size_t i = 0; i < x_target.size(); ++i) {
        double x = x_target[i];
        double idx0 = x / delta_x;
        int j0 = static_cast<int>(std::floor(idx0));
        double delta_idx = idx0 - j0;

        // Build periodic index list
        std::vector<double> y_local(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (j0 - half + k) % N;
            if (idx < 0) idx += N;
            y_local[k] = ygrid[idx];
        }

        // Local Lagrange interpolation
        double p = 0.0;
        for (int j = 0; j <= order; ++j) {
            double Lj = lagrange_basis(delta_idx, x_local, j);
            p += y_local[j] * Lj;
        }
        f_interp[i] = p;
    }
    return f_interp;
}

// MEX gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check number of inputs
    if (nrhs != 4)
        mexErrMsgIdAndTxt("lagrange1d_local_interp_periodic:nrhs", "Four inputs required: x_target, xgrid, ygrid, order");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("lagrange1d_local_interp_periodic:nlhs", "One output allowed.");

    // Get input arguments
    size_t n_target = mxGetNumberOfElements(prhs[0]);
    size_t n_grid = mxGetNumberOfElements(prhs[1]);
    size_t n_ygrid = mxGetNumberOfElements(prhs[2]);
    if (n_grid != n_ygrid)
        mexErrMsgIdAndTxt("lagrange1d_local_interp_periodic:input", "xgrid and ygrid must have the same length.");

    // Convert x_target
    double* x_target_ptr = mxGetPr(prhs[0]);
    std::vector<double> x_target(x_target_ptr, x_target_ptr + n_target);

    // Convert xgrid
    double* xgrid_ptr = mxGetPr(prhs[1]);
    std::vector<double> xgrid(xgrid_ptr, xgrid_ptr + n_grid);

    // Convert ygrid
    double* ygrid_ptr = mxGetPr(prhs[2]);
    std::vector<double> ygrid(ygrid_ptr, ygrid_ptr + n_ygrid);

    // Convert order
    int order = static_cast<int>(mxGetScalar(prhs[3]));

    // Call the interpolation function
    std::vector<double> f_interp;
    try {
        f_interp = lagrange_local_interp_periodic(x_target, xgrid, ygrid, order);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("lagrange1d_local_interp_periodic:exception", e.what());
    }

    // Set output
    plhs[0] = mxCreateDoubleMatrix(n_target, 1, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);
    for (size_t i = 0; i < n_target; ++i)
        out_ptr[i] = f_interp[i];
}
