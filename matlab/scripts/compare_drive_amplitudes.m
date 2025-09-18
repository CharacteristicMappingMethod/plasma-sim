% Compare external electric field amplitudes for weak and canonical drives
clear; clc; close all;

% Add paths
addpath(genpath('../src/'), genpath('../params/'), "../");
DEFAULTS

% Load parameters for both drives
PARAMS_keen_waves_weak;
params_weak = params;

PARAMS_keen_waves_canonical;
params_canonical = params;

% Time vector for plotting
t = linspace(0, 500, 1000);

% Compute a(t) for both drives
a_weak = params_weak.a(t);
a_canonical = params_canonical.a(t);

% Create figure
figure('Position', [100, 100, 1000, 600]);

% Plot both amplitude functions
plot(t, a_weak, 'b-', 'LineWidth', 2, 'DisplayName', 'Weak Drive');
hold on;
plot(t, a_canonical, 'r-', 'LineWidth', 2, 'DisplayName', 'Canonical Drive');

% Add vertical lines for key timing parameters
% Weak drive timing
xline(params_weak.tL, 'b--', 'LineWidth', 1, 'DisplayName', 'Weak: Left ramp start');
xline(params_weak.tR, 'b--', 'LineWidth', 1, 'DisplayName', 'Weak: Right ramp start');

% Canonical drive timing
xline(params_canonical.tL, 'r--', 'LineWidth', 1, 'DisplayName', 'Canonical: Left ramp start');
xline(params_canonical.tR, 'r--', 'LineWidth', 1, 'DisplayName', 'Canonical: Right ramp start');

% Formatting
xlabel('Time $t$', 'Interpreter', 'latex');
ylabel('Amplitude $a(t)$', 'Interpreter', 'latex');
title('External Electric Field Amplitude Comparison', 'Interpreter', 'latex');
legend("Weak drive", "Canonical Drive",'Location', 'best');
grid on;
xlim([0, 500]);

% Save figure
save_fig_tikz("../images/drive_amplitude_comparison");


