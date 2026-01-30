
function [params]=measure(params, fs)


iT = params.it;
time = params.dt * iT;
Efield = params.Efield;

% Loop over species
for s = 1:params.Ns
    % Get species-specific grid and distribution function
    grid = params.grids(s);
    f = fs(:,:,s);
    Vgrid = grid.Vsample_grid;
    dv = grid.dv;
    % Calculate diagnostics
    density = compute_density(f,grid.dv);
    Mass = sum(density, "all") * grid.dx;
    Momentum = sum((f .* Vgrid) .* dv(:), "all") * grid.dx;
    Epot = 0.5 * sum(Efield.^2) * grid.dx;
    Ekin = 0.5 * sum(f .* (Vgrid.^2).*dv(:), "all") * grid.dx;
    Etot = Epot + Ekin;
    L2norm = sum((abs(f).^2).*dv(:), "all") * grid.dx;
    
    if isfield(params, 'incomp_error')
        incomp_error = params.incomp_error(s);
    else
        incomp_error = 0;
    end
    % Create a filename for the species
    species_name = params.species_name(s); % Species name
    filename = fullfile(params.data_dir, species_name+".csv");

    rho_modes = fourier_modes(density, 5);
    % Create a table row for the current measurement
    new_row = table(iT, time, Mass, Momentum, Epot, Ekin, Etot, L2norm, incomp_error, rho_modes(1), rho_modes(2), rho_modes(3), rho_modes(4), rho_modes(5),...
        'VariableNames', {'it','time','Mass', 'Momentum', 'Epot', 'Ekin', 'Etot', 'L2norm', 'incomp_error', 'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

    % Check if the file already exists
    if exist(filename, 'file') && iT>1
        % Load the existing table
        existing_table = readtable(filename);
        % Append the new row
        updated_table = [existing_table; new_row];
    else
        % Create a new table with the current row
        updated_table = new_row;
        params.diagnostic_filename(s) = filename;
    end

    % Write the updated table to the file
    writetable(updated_table, filename);
end
end


function [mode_abs] = fourier_modes(field1D, nmodes)

field1D_fft = fft(field1D);
mode_abs = abs(field1D_fft(2:nmodes+1));

end