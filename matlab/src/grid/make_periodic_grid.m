function [grid] = make_periodic_grid(Lx,Lv,Nx,Nv)

    x = [0:Nx-1] * Lx / Nx; % Spatial grid
    v = [0:Nv-1] * 2*Lv / Nv-Lv; % Spatial grid
    %v = linspace(-Lv,Lv,Nv);

    dx = x(2) - x(1); % Spatial grid spacing
    dv = v(2) - v(1); % Velocity grid spacing

    % Create meshgrid for spatial and velocity coordinates
    [X,V] = meshgrid(x, v);

    grid.x  = x; grid.v  = v;
    grid.X  = X; grid.V  = V;
    grid.Xsample_grid = X; grid.Vsample_grid = V; % sample grid is the grid at which f is evaluated.
    grid.size_sample_grid = size(X);
    grid.dx = dx;grid.dv = dv;
    grid.Lx = Lx;grid.Lv = Lv;
    grid.Nx = Nx;grid.Nv = Nv;


    % Create finite difference operators
    % Create 2nd order central difference operators
    grid.Dx = spdiags(ones(Nx,1)*[-1,0,1], [-1,0,1], Nx, Nx);
    grid.Dv = spdiags(ones(Nv,1)*[-1,0,1], [-1,0,1], Nv, Nv);

    % Apply periodic boundary conditions
    grid.Dx(1,end) = -1;  % Connect first and last points
    grid.Dx(end,1) = 1;
    grid.Dv(1,end) = -1;
    grid.Dv(end,1) = 1;

    % Scale by grid spacing
    grid.Dx = grid.Dx / (2*dx);
    grid.Dv = grid.Dv / (2*dv);
    
    grid.size = size(X);
    grid.dom = [0, -Lv, Lx-dx, Lv-dv];
    
    % Fourier
    kx = (2*pi/Lx) * [-Nx/2:Nx/2-1];
    grid.kx = fftshift(kx);
    grid.kx2 = grid.kx.^2;
    grid.kx2(1) = 1;
    kv = (pi/Lv) * [-Nv/2:Nv/2-1];
    grid.kv = fftshift(kv);
    grid.kv2 = grid.kv.^2;
    grid.kv2(1) = 1;

    % periodic flow velocity
    [vper,sigma] = velocity_periodicfication(grid);
    [~,Vperiodic] = meshgrid(x,vper);
    [~,Weights] = meshgrid(x,abs(sigma)<1e-12);
    grid.Vperiodic = Vperiodic;
    grid.Weights = Weights;
    
end

