#include "mex.h"
#include <vector>
#include <cmath>
#include <stdexcept>

// Helper: Barycentric weights
std::vector<double> barycentric_weights(const std::vector<double>& x_nodes) {
    int n = x_nodes.size();
    std::vector<double> w(n, 1.0);
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            if (k != j) {
                w[j] /= (x_nodes[j] - x_nodes[k]);
            }
        }
    }
    return w;
}

// Helper: Barycentric interpolation
double barycentric_interp(double x, const std::vector<double>& x_nodes,
                          const std::vector<double>& y_nodes,
                          const std::vector<double>& w) {
    double num = 0.0, denom = 0.0;
    for (size_t j = 0; j < x_nodes.size(); ++j) {
        double dx = x - x_nodes[j];
        if (std::abs(dx) < 1e-14)
            return y_nodes[j];  // exact match
        double term = w[j] / dx;
        num += term * y_nodes[j];
        denom += term;
    }
    return num / denom;
}

// Main 2D interpolation
std::vector<double> lagrange2d_barycentric_interp(
    const std::vector<double>& x_target,
    const std::vector<double>& y_target,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid,
    const std::vector<double>& Fgrid,
    int Nx, int Ny,
    int order
) {
    int Nt = x_target.size();
    double dx = xgrid[1] - xgrid[0];
    double dy = ygrid[1] - ygrid[0];
    int half = order / 2;

    std::vector<double> F_interp(Nt, 0.0);
    std::vector<double> local_nodes(order + 1);
    for (int k = 0; k <= order; ++k)
        local_nodes[k] = -half + k;

    std::vector<double> w = barycentric_weights(local_nodes);

    for (int i = 0; i < Nt; ++i) {
        double x = x_target[i];
        double y = y_target[i];

        int jx = static_cast<int>(std::floor(x / dx));
        double delta_x = x / dx - jx;
        std::vector<int> idxs_x(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (jx - half + k) % Nx;
            if (idx < 0) idx += Nx;
            idxs_x[k] = idx;
        }

        int jy = static_cast<int>(std::floor(y / dy));
        double delta_y = y / dy - jy;
        std::vector<int> idxs_y(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (jy - half + k) % Ny;
            if (idx < 0) idx += Ny;
            idxs_y[k] = idx;
        }

        std::vector<double> fx(order + 1);
        for (int jj = 0; jj <= order; ++jj) {
            std::vector<double> row(order + 1);
            for (int j = 0; j <= order; ++j) {
                row[j] = Fgrid[idxs_y[jj] + idxs_x[j] * Ny];
            }
            fx[jj] = barycentric_interp(delta_x, local_nodes, row, w);
        }

        F_interp[i] = barycentric_interp(delta_y, local_nodes, fx, w);
    }

    return F_interp;
}

// MEX entry point
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {
    if (nrhs != 6)
        mexErrMsgTxt("Usage: F_interp = lagrange2d_local_interp_barycentric(xq, yq, xgrid, ygrid, Fgrid, order)");

    // Inputs
    double* xq_ptr = mxGetPr(prhs[0]);
    double* yq_ptr = mxGetPr(prhs[1]);
    double* xgrid_ptr = mxGetPr(prhs[2]);
    double* ygrid_ptr = mxGetPr(prhs[3]);
    double* F_ptr = mxGetPr(prhs[4]);
    int order = static_cast<int>(mxGetScalar(prhs[5]));

    int Nt = static_cast<int>(mxGetNumberOfElements(prhs[0]));
    int Nx = static_cast<int>(mxGetNumberOfElements(prhs[2]));
    int Ny = static_cast<int>(mxGetNumberOfElements(prhs[3]));

    // Copy to std::vector
    std::vector<double> xq(xq_ptr, xq_ptr + Nt);
    std::vector<double> yq(yq_ptr, yq_ptr + Nt);
    std::vector<double> xgrid(xgrid_ptr, xgrid_ptr + Nx);
    std::vector<double> ygrid(ygrid_ptr, ygrid_ptr + Ny);
    std::vector<double> F(F_ptr, F_ptr + Nx * Ny);

    // Run interpolation
    std::vector<double> F_interp = lagrange2d_barycentric_interp(xq, yq, xgrid, ygrid, F, Nx, Ny, order);

    // Output
    plhs[0] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
    std::copy(F_interp.begin(), F_interp.end(), mxGetPr(plhs[0]));
}
