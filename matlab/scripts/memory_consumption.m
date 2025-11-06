clear all;
clc
close all;
% Add paths
addpath(genpath('../src/'),genpath('../params/'),"../");
DEFAULTS
% Memory consumption comparison for NuFi and CMM-NuFi methods
%pic_dir = pic_dir;
% Parameters
N_sample = 1024;
Nnufi_max = 500;

% N range (total iterations)
N_vec = linspace(1,Nnufi_max);

% Fixed parameter
d = 2;

% Parameter arrays to loop over (combinations of N_remap and N_mapgrid)
N_remap_values = [5, 20, 100];
N_mapgrid_values = [64, 64, 64];

% Plot
figure(4);

% Calculate number of combinations
num_combinations = length(N_remap_values) * length(N_mapgrid_values);
colors = lines(num_combinations);
idx = 1;

% Calculate NuFi memory (only once since d is constant)
memory_nufi = (N_vec + 1) * N_sample^d * 8;

% Plot NuFi
plot(N_vec, memory_nufi / (1024)^2, 'k-', 'LineWidth', 2, ...
     'DisplayName', sprintf('NuFi $N_f=%d$ ($d=%d$)',N_sample, d));
hold on;

% Loop over all combinations of N_remap and N_mapgrid
for i = 1:length(N_remap_values)
    N_remap = N_remap_values(i);
        N_mapgrid = N_mapgrid_values(i);
        
        % Calculate N_maps
        N_maps = floor(N_vec / N_remap);
        
        % Calculate CMM-NuFi memory
        memory_cmm_nufi = (N_maps * 2*d * N_mapgrid^(2*d) + ...
                          (N_vec - N_maps*N_remap + 1) * N_sample^d) * 8;
        
        % Plot CMM-NuFi
        plot(N_vec, memory_cmm_nufi / (1024)^2, '--', 'LineWidth', 1.5, ...
             'Color', colors(idx,:), 'DisplayName', ...
             sprintf('CMM-NuFi ($N_\\mathrm{remap}=%d$, $N_{\\boldsymbol{\\chi}}=%d$)', ...
                     N_remap, N_mapgrid));
        idx = idx + 1;
end

ylim([0,2*10^4])
xlabel('$N$ (total iterations)');
ylabel('Memory [MB]');
title(sprintf('Memory Consumption $d=%d$',d));
legend('Location', 'northeast');
grid on;


save_fig_tikz(pic_dir+ "/memory_d"+num2str(d));
