%PENALTY Solves a constrained optimization problem.
%   PENALTY(F,C,X0) minimizes the the objective function F subject to the
%   constraints C starting from an initial guess X0. The optimization uses
%   MATLAB's fminsearch routine, a derivative-free method that uses the
%   Nelder-Mead downhill simplex algorithm. Constraints are handled via a
%   penalty function that iteratively enforces the constraints. Inequality
%   constraints should be written in the form c(x) >= 0. Equality
%   constraints can be included by using two inequalities.
%
%   Current implementation uses a very basic heuristic for determining the
%   initial penalty, and each iteration increases the penalty by a factor
%   of 10.
%
%   Outputs:
%   X, the optimal point
%   FVAL, the value of the objective function at X
%   FUNCCOUNT, the total number of function evaluations over all iterations
%   S, a structure with the initial and final penalty values

function [x,fval,funcCount,S] = Penalty(f,c,x0)
if isempty(c)
    c = @(x) 0;
    c0 = 1;
    warning('In the absence of constraints, Penalty is equivalent to using fminsearch.')
    iterMax = 1;
else
    c0 = norm(min(0,c(x0)).^2);
    iterMax = 10;
end
if ~c0
    c0 = 1;
end
f0 = f(x0);
if ~f0
    f0 = 1;
end
sigma = 4*f0/c0;
S.initial = sigma;
x = x0;
tol = 1e-6;
xold = x0+10*tol;
funcCount = 0;
iter = 0;

while norm(x-xold) > tol && iter < iterMax
    phi = @(x) f(x) + sum(sigma*min(0,c(x)).^2);
    xold = x;
    [x,fval,~,output] = fminsearch(phi,xold);
    sigma = sigma*10;
    funcCount = output.funcCount + funcCount;
    iter = iter + 1;
end
S.final = sigma/10; %Since it gets multiplied by 10 before actually being used again
end