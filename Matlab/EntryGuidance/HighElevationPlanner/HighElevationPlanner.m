%High Elevation Planner
clear; clc; close all

% p = [15;130;160]; %Works well for 780
p = [5,51,250]; % 12km
opt = optimset('tolX',1e-3,'tolFun',1e-8);
mars = Mars();
vm = VehicleModel();
[p,fval,flag,output] = fminsearch(@(p) HighElevationCostFunction(p,mars,vm), p, opt);
[cost,t,x] = HighElevationCostFunction(p, mars, vm);

EntryPlots(t,x)
close all