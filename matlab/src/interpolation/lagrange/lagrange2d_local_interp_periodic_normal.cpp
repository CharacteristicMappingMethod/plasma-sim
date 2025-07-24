#include "mex.h"
#include <vector>
#include <cmath>
#include <stdexcept>

// Helper: computes the j-th Lagrange basis polynomial evaluated at x
inline double lagrange_basis(double x, const std::vector<double>& x_nodes, int j) {
    double Lj = 1.0;
    int n = x_nodes.size();
    for (int m = 0; m < n; ++m) {
        if (m != j) {
            Lj *= (x - x_nodes[m]) / (x_nodes[j] - x_nodes[m] + 1e-32);
        }
    }
    return Lj;
}

// 2D local Lagrange interpolation with periodic boundaries
std::vector<double> lagrange2d_local_interp_periodic(
    const std::vector<double>& x_target,
    const std::vector<double>& y_target,
    const std::vector<double>& xgrid,
    const std::vector<double>& ygrid,
    const std::vector<double>& Fgrid, // flattened column-major (MATLAB style)
    int Nx, int Ny,
    int order
) {
    int Nt = x_target.size();
    double dx = xgrid[1] - xgrid[0];
    double dy = ygrid[1] - ygrid[0];
    int half = order / 2;
    std::vector<double> local_nodes(order + 1);
    for (int k = 0; k <= order; ++k) local_nodes[k] = -half + k;
    std::vector<double> F_interp(Nt, 0.0);

    for (int i = 0; i < Nt; ++i) {
        double x = x_target[i];
        double y = y_target[i];
        // X indices
        double idx_x = x / dx;
        int jx = static_cast<int>(std::floor(idx_x));
        double delta_x = idx_x - jx;
        std::vector<int> idxs_x(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (jx - half + k) % Nx;
            if (idx < 0) idx += Nx;
            idxs_x[k] = idx;
        }
        // Y indices
        double idx_y = y / dy;
        int jy = static_cast<int>(std::floor(idx_y));
        double delta_y = idx_y - jy;
        std::vector<int> idxs_y(order + 1);
        for (int k = 0; k <= order; ++k) {
            int idx = (jy - half + k) % Ny;
            if (idx < 0) idx += Ny;
            idxs_y[k] = idx;
        }
        // Interpolate in x for each row in y
        std::vector<double> fx(order + 1, 0.0);
        for (int jj = 0; jj <= order; ++jj) {
            // row = Fgrid(idxs_y[jj], idxs_x)
            std::vector<double> row(order + 1);
            for (int j = 0; j <= order; ++j) {
                // MATLAB: Fgrid(row, col) = Fgrid(y, x), column-major
                row[j] = Fgrid[idxs_y[jj] + idxs_x[j] * Ny];
            }
            fx[jj] = 0.0;
            for (int j = 0; j <= order; ++j) {
                double Ljx = lagrange_basis(delta_x, local_nodes, j);
                fx[jj] += row[j] * Ljx;
            }
        }
        // Interpolate result in y
        double fy = 0.0;
        for (int j = 0; j <= order; ++j) {
            double Ljy = lagrange_basis(delta_y, local_nodes, j);
            fy += fx[j] * Ljy;
        }
        F_interp[i] = fy;
    }
    return F_interp;
}

// MEX gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 6)
        mexErrMsgIdAndTxt("lagrange2d_local_interp_periodic:nrhs", "Six inputs required: x_target, y_target, xgrid, ygrid, Fgrid, order");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("lagrange2d_local_interp_periodic:nlhs", "One output allowed.");

    // Get input arguments
    size_t n_target = mxGetNumberOfElements(prhs[0]);
    size_t n_target_y = mxGetNumberOfElements(prhs[1]);
    if (n_target != n_target_y)
        mexErrMsgIdAndTxt("lagrange2d_local_interp_periodic:input", "x_target and y_target must have the same length.");
    size_t Nx = mxGetNumberOfElements(prhs[2]);
    size_t Ny = mxGetNumberOfElements(prhs[3]);
    size_t nF = mxGetNumberOfElements(prhs[4]);
    if (nF != Nx * Ny)
        mexErrMsgIdAndTxt("lagrange2d_local_interp_periodic:input", "Fgrid size does not match xgrid and ygrid.");

    // Convert x_target, y_target
    double* x_target_ptr = mxGetPr(prhs[0]);
    double* y_target_ptr = mxGetPr(prhs[1]);
    std::vector<double> x_target(x_target_ptr, x_target_ptr + n_target);
    std::vector<double> y_target(y_target_ptr, y_target_ptr + n_target);
    // Convert xgrid, ygrid
    double* xgrid_ptr = mxGetPr(prhs[2]);
    double* ygrid_ptr = mxGetPr(prhs[3]);
    std::vector<double> xgrid(xgrid_ptr, xgrid_ptr + Nx);
    std::vector<double> ygrid(ygrid_ptr, ygrid_ptr + Ny);
    // Convert Fgrid (column-major)
    double* Fgrid_ptr = mxGetPr(prhs[4]);
    std::vector<double> Fgrid(Fgrid_ptr, Fgrid_ptr + nF);
    // Convert order
    int order = static_cast<int>(mxGetScalar(prhs[5]));

    // Call the interpolation function
    std::vector<double> F_interp;
    try {
        F_interp = lagrange2d_local_interp_periodic(x_target, y_target, xgrid, ygrid, Fgrid, Nx, Ny, order);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("lagrange2d_local_interp_periodic:exception", e.what());
    }

    // Set output
    plhs[0] = mxCreateDoubleMatrix(n_target, 1, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);
    for (size_t i = 0; i < n_target; ++i)
        out_ptr[i] = F_interp[i];
}