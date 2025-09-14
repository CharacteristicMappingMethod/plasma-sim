
function [params]=measure(params, fs)


iT = params.it;
time = params.dt * iT;
Efield = params.Efield;

% Loop over species
for s = 1:params.Ns
    % Get species-specific grid and distribution function
    grid = params.grids(s);
    f = fs(:,:,s);
    Vgrid = grid.V;

    % Calculate diagnostics
    Mass = sum(f, "all") * grid.dx * grid.dv;
    Momentum = sum(f .* Vgrid, "all") * grid.dx * grid.dv;
    Epot = 0.5 * sum(Efield.^2) * grid.dx;
    Ekin = 0.5 * sum(f .* (Vgrid.^2), "all") * grid.dx * grid.dv;
    Etot = Epot + Ekin;
    L2norm = sum(abs(f).^2, "all") * grid.dx * grid.dv;
    
    if isfield(params, 'incomp_error')
        incomp_error = params.incomp_error(s,iT);
    else
        incomp_error = 0;
    end
incomp_error
    % Create a filename for the species
    species_name = params.species_name(s); % Species name
    filename = fullfile(params.data_dir, species_name+".csv");

    % Create a table row for the current measurement
    new_row = table(iT, time, Mass, Momentum, Epot, Ekin, Etot, L2norm, incomp_error, ...
        'VariableNames', {'it','time','Mass', 'Momentum', 'Epot', 'Ekin', 'Etot', 'L2norm', 'incomp_error'});

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


function [Emode_abs] = fourier_modes(Efield, k_list, grid)
% this function implements (63-65) of https://arxiv.org/pdf/1009.3046.pdf
% for odd fourier numbers (k=0.5 etc)
ik = 1;
Ek = fft(Efield);
Emode_abs = abs(Ek(2:5))/grid.size(1);
%     for k = k_list
%         Esin_k = sum(Efield.*sin(2*pi*k*grid.x/grid.L(1)))*grid.dx;
%         Ecos_k = sum(Efield.*cos(2*pi*k*grid.x/grid.L(1)))*grid.dx;
%         Emode_abs(ik) = sqrt(abs(Esin_k)^2+abs(Ecos_k)^2)/grid.L(1);
%         ik = ik + 1;
%     end

end