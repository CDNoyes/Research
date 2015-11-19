%High Elevation Planner
clear; clc; close all

% p = [15;130;160]; %Works well for 780

%Initial Guess
p = [5,51,250]; % 12km

%Setup the optimization tolerances
opt = optimset('tolX',1e-3,'tolFun',1e-8);

%Create the standard models to be used
mars = Mars();
vm = VehicleModel();

%Define the target location
DR = 780;
CR = 0;

%Call the optimization routine
[p,fval,flag,output] = fminsearch(@(p) HighElevationCostFunction(p,mars,vm,DR,CR), p, opt);
%Compute the state trajectory
[cost,t,x] = HighElevationCostFunction(p, mars, vm,DR,CR);

EntryAnalysis(t,x)
% EntryPlots(t,x)
% close all