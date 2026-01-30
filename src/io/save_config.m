function data = save_config(params,data,fs,Nsamples,save_to_file)

if nargin < 5  % Check if save_to_file is provided
    save_to_file = 0;
end


for species = 1:params.Ns
    data.fs(:,:,Nsamples,species) = fs(:,:,species);
    data.Efield(:,Nsamples) = params.Efield;
end

% Save data to file if requested
if save_to_file == 1
    if isfield(params, 'data_dir')
        filename = fullfile(params.data_dir, 'config_data.mat');
        save(filename, 'data',"params");
        fprintf('Data saved to %s\n', filename);
    else
        error('params.data_dir is not defined or is not a valid string.');
    end
end
end