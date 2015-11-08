%Rosenbrock "banana" function
clear; clc;

f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% [x,fopt,iter] = TrustRegion(f,[0;0]);

[x,fval,flag,output] = fminsearch(f,[0;0]);

c = @(x) [x(2); -x(2)];
[x_constrained, f_constrained, f_evals, penalty] = Penalty(f,c,[1;1]);