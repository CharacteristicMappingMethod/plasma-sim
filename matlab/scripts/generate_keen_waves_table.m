% Generate LaTeX table comparing Keen wave parameters
clear; clc; close all;

% Add paths
addpath(genpath('../src/'), genpath('../params/'), "../");
DEFAULTS

% Load parameters for both drives
PARAMS_keen_waves_weak;
params_weak = params;

PARAMS_keen_waves_canonical;
params_canonical = params;

% Create LaTeX table
fprintf('Generating LaTeX table for Keen wave parameters...\n');

% Open file for writing
fid = fopen('../images/keen_waves_parameters_table.tex', 'w');

% Write LaTeX table header
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Keen Wave Parameters Comparison}\n');
fprintf(fid, '\\label{tabl:keen_waves_params}\n');
fprintf(fid, '\\begin{tabular}{lcc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Parameter & Weak Drive & Canonical Drive \\\\\n');
fprintf(fid, '\\midrule\n');

% Grid parameters
fprintf(fid, '\\multicolumn{3}{l}{\\textbf{Grid Parameters}} \\\\\n');
fprintf(fid, '$N_x \\times N_v$ & $%d\\times%d$ & $%d\\times%d$ \\\\\n', params_weak.Nx,params_weak.Nv, params_canonical.Nx,params_canonical.Nv);
fprintf(fid, '$L_x$ & $[%.1f,%.1f]$ & $[%.1f,%.1f]$ \\\\\n', 0,params_weak.Lx, 0, params_canonical.Lx);
fprintf(fid, '$L_v$ & $[%.1f,%.1f]$ & $[%.1f,%.1f]$ \\\\\n', -params_weak.Lv,params_weak.Lv, -params_canonical.Lv, params_canonical.Lv);
fprintf(fid, '\\addlinespace\n');

% Time parameters
fprintf(fid, '\\multicolumn{3}{l}{\\textbf{Time Parameters}} \\\\\n');
fprintf(fid, '$T_\\mathrm{end}$ & %.0f & %.0f \\\\\n', params_weak.Tend, params_canonical.Tend);
fprintf(fid, '$\\Delta t$ & %.1f & %.1f \\\\\n', params_weak.dt, params_canonical.dt);
fprintf(fid, '$\\Delta t_\\mathrm{save}$ & %.0f & %.0f \\\\\n', params_weak.dt_save, params_canonical.dt_save);
fprintf(fid, '$N_\\mathrm{remap}$ & %d & %d \\\\\n', params_weak.N_remap, params_canonical.N_remap);
fprintf(fid, '\\addlinespace\n');

% Drive parameters
fprintf(fid, '\\multicolumn{3}{l}{\\textbf{Drive Parameters}} \\\\\n');
fprintf(fid, '$k_\\mathrm{Dr}$ & %.3f & %.3f \\\\\n', params_weak.kDr, params_canonical.kDr);
fprintf(fid, '$\\omega_\\mathrm{Dr}$ & %.3f & %.3f \\\\\n', params_weak.wDr, params_canonical.wDr);
fprintf(fid, '$a_\\mathrm{Dr}$ & %.3f & %.3f \\\\\n', params_weak.aDr, params_canonical.aDr);
fprintf(fid, '$T_\\mathrm{Dr}$ & %.0f & %.0f \\\\\n', params_weak.TDr, params_canonical.TDr);
fprintf(fid, '$t_\\mathrm{L}$ & %.0f & %.0f \\\\\n', params_weak.tL, params_canonical.tL);
fprintf(fid, '$t_\\text{w}$ & %.0f & %.0f \\\\\n', params_weak.twL, params_canonical.twL);
fprintf(fid, '$t_\\mathrm{w}$ & %.0f & %.0f \\\\\n', params_weak.twR, params_canonical.twR);
fprintf(fid, '$t_\\text{R}$ & %.0f & %.0f \\\\\n', params_weak.tR, params_canonical.tR);
fprintf(fid, '\\addlinespace\n');



% Table footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

% Close file
fclose(fid);

fprintf('LaTeX table saved to: ../images/keen_waves_parameters_table.tex\n');

% Also display the table in console
fprintf('\nGenerated LaTeX table:\n');
fprintf('=====================\n\n');

% Read and display the file
fid = fopen('../images/keen_waves_parameters_table.tex', 'r');
while ~feof(fid)
    line = fgetl(fid);
    fprintf('%s\n', line);
end
fclose(fid);
