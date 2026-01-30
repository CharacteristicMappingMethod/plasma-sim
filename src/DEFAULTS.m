%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab using latex
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'DefaultAxesTitle','latex');
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load library and parameters
%addpath(genpath('../../lib/'),genpath('../../params/'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directories
pic_dir="./images/";
display("[Defaults] images will be stored in: "+ pic_dir)

dat_dir="./data/";
display("[Defaults] images will be stored in: "+ dat_dir)


% Set default axes properties
set(0, 'DefaultAxesFontSize', 14);        % Larger font size for axes
set(0, 'DefaultAxesLineWidth', 1.5);      % Thicker axes lines
set(0, 'DefaultAxesXMinorTick', 'on');    % Enable minor ticks
set(0, 'DefaultAxesYMinorTick', 'on');

% Set default line properties for plots
set(0, 'DefaultLineLineWidth', 2);        % Thicker plot lines
set(0, 'DefaultLineMarkerSize', 8);       % Larger markers

% Set default legend properties
set(0, 'DefaultLegendFontSize', 14);      % Larger font for legends