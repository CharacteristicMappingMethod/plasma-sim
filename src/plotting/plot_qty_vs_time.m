function h = plot_qty_vs_time(params, qty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a quantity (qty) from a given simulation
% Input:    qty    ... Quantity as a string (e.g. "Epot", "Etot")
%           params ... structure of parameters including params.data_dir
% Output:   plot of the quantity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


species_colors = lines(params.Ns); % Assign different colors for each species

for s = 1:params.Ns
    species_name = params.species_name(s); % Species name
    filename = fullfile(params.data_dir, species_name+ '.csv');

    if exist(filename, 'file')
        filename
        % Read the diagnostic data
        data = readtable(filename);
        time = data.time; % Assuming each row is a new iteration
        data_spec{s} = data;
        if qty == "Etot"
            h(s) = semilogy(time, data.Etot, 'Color', species_colors(s,:), 'DisplayName', species_name);
            title('Total Energy'); xlabel('time'); ylabel('Etot'); grid on; hold on
        elseif qty == "Epot"
            h(s) = semilogy(time, data.Epot);
            title('Potential Energy'); xlabel('time'); ylabel('$E_\mathrm{pot}$'); grid on; hold on
            return
        end

    else
        warning('File %s does not exist. Skipping.', filename);
    end

end


end