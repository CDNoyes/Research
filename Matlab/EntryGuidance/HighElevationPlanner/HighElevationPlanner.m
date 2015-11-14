%High Elevation Planner
clear
p = [15;130;160]; %Works well for 780
opt = optimset('tolX',1e-3,'tolFun',1e-8);
[p,fval,flag,output] = fminsearch(@HighElevationCostFunction, p, opt);
[cost,t,x] = HighElevationCostFunction(p);

% EntryPlots(t,x)