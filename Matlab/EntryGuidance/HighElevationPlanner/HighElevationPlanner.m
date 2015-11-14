%High Elevation Planner
clear
p = [15;130;160];
opt = optimset('tolX',1e-6,'tolFun',1e-6);
[p,fval,flag,output] = fminsearch(@HighElevationCostFunction, p, opt);
[cost,t,x] = HighElevationCostFunction(p);

% EntryPlots(t,x)