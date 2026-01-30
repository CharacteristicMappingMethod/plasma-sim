%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main script for CMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('./src/'),genpath('./params/'))
DEFAULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select case:
PARAMS_two_stream;
%PARAMS_landau_damping;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate
tic()
params.method = "CMM-NuFI";
[params1, data] = Sim(params);
t(1) = toc()

tic()
params.method = "NuFi";
[params2, data] = Sim(params);
t(2) = toc()

fprintf("tcpu CMM-Nufi: %2.2f\n", t(1))
fprintf("tcpu NuFI: %2.2f \n", t(2))



%% some plotting
figure(44)
h = plot_qty_vs_time(params1,"Epot");
%figure(45)
hold on
h = plot_qty_vs_time(params2,"Epot");
h.Color="r"
h.LineStyle="--"
legend(params1.method,params2.method)