#include "mex.h"
#include <vector>
#include <cmath>
#include <stdexcept>

// Compute barycentric weights for 1D nodes
std::vector<double> barycentric_weights(const std::vector<double>& x_nodes) {
    int n = x_nodes.size();
    std::vector<double> w(n, 1.0);
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            if (k != j)
                w[j] /= (x_nodes[j] - x_nodes[k]);
        }
    }
    return w;
}

// Evaluate barycentric interpolant
double barycentric_interp(double x, const std::vector<double>& x_nodes,
                          const std::vector<double>& y_nodes,
                          const std::vector<double>& w) {
    double num = 0.0, denom = 0.0;
    for (size_t j = 0; j < x_nodes.size(); ++j) {
        double dx = x - x_nodes[j];
        if (std::abs(dx) < 1e-14)
            return y_nodes[j];
        double t = w[j] / dx;
        num += t * y_nodes[j];
        denom += t;
    }
    return num / denom;
}

// Main 1D periodic local barycentric interpolation
std::vector<double> lagrange1d_barycentric_interp(
    const std::vector<double>& x_target,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid,
    int order
) {
    int Nt = x_target.size();
    int N = xgrid.size();
    double dx = xgrid[1] - xgrid[0];
    int half = order / 2;

    std::vector<double> f_interp(Nt, 0.0);
    std::vector<double> local_nodes(order + 1);
    for (int k = 0; k <= order; ++k)
        local_nodes[k] = -half + k;

    std::vector<double> w = barycentric_weights(local_nodes);

    for (int i = 0; i < Nt; ++i) {
        double x = x_target[i];
        double idx_f = x / dx;
        int j0 = static_cast<int>(std::floor(idx_f));
        double delta = idx_f - j0;

        std::vector<int> idxs(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (j0 - half + k) % N;
            if (idx < 0) idx += N;
            idxs[k] = idx;
        }

        std::vector<double> y_local(order + 1);
        for (int k = 0; k <= order; ++k)
            y_local[k] = ygrid[idxs[k]];

        f_interp[i] = barycentric_interp(delta, local_nodes, y_local, w);
    }

    return f_interp;
}

// MEX entry point
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    if (nrhs != 4)
        mexErrMsgTxt("Usage: f_interp = lagrange1d_local_interp_barycentric(xq, xgrid, fgrid, order)");

    // Inputs
    double* xq_ptr = mxGetPr(prhs[0]);
    double* xgrid_ptr = mxGetPr(prhs[1]);
    double* fgrid_ptr = mxGetPr(prhs[2]);
    int order = static_cast<int>(mxGetScalar(prhs[3]));

    int Nt = static_cast<int>(mxGetNumberOfElements(prhs[0]));
    int N = static_cast<int>(mxGetNumberOfElements(prhs[1]));

    // Copy to std::vector
    std::vector<double> xq(xq_ptr, xq_ptr + Nt);
    std::vector<double> xgrid(xgrid_ptr, xgrid_ptr + N);
    std::vector<double> fgrid(fgrid_ptr, fgrid_ptr + N);

    // Call interpolation
    std::vector<double> f_interp = lagrange1d_barycentric_interp(xq, xgrid, fgrid, order);

    // Output
    plhs[0] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
    std::copy(f_interp.begin(), f_interp.end(), mxGetPr(plhs[0]));
}